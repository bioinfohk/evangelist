import requests, sys
from ftplib import FTP
from dna_puller.sliding_parser import SlidingParser
import os, shutil
import gzip
import json, time


class DnaPuller:
    def __init__(self, species, jsons = True, remove_data = True, types = ['dna', 'cdna', 'cds', 'dna_sm'], parse = True, window_size = 1000, jsons_dir = 'jsons'):
        self.species = species
        self.types = types
        self.data = {}
        self.jsons = jsons
        self.remove_data = remove_data
        self.parse = parse
        self.jsons_dir = jsons_dir
        self.parser = SlidingParser(window_size)

    def download_and_parse_data(self):
        for speciess in self.species:
            
            if not os.path.exists(speciess):
              os.mkdir(speciess)
            
            self.data[speciess] = {}
            for type in self.types:
                self.actual_type = type
                ftp = FTP('ftp.ensembl.org')
                ftp.login()
                print('getting ' + type + ' for ' + speciess)
                self.data[speciess][type] = {}
                
                dirs = self.ftp_cwd(speciess.lower(), type)
                
                ftp.cwd(self.ftp_cwd(speciess.lower(), type))
                
                # Get filenames from Ensembl FTP server and validate them if they're correct filenames for our  
                # analysis (in this point linkage groups) based on rules in validate_name method
                self.filenames = []
                ftp.retrlines('NLST', callback=self.validate_name)

                # If no filename is valid, try with less restrictive rules.
                if len(self.filenames) == 0:
                      ftp.retrlines('NLST', callback=self.validate_name_toplevel)

                # Filtration of files by dna type - In Ensembl, there is 'dna_sm', 'dna_rm' and 'dna' in dna folder.
                # That means that previous step gives me all linkage groups three times - for 'dna_sm', 'dna_rm' and 'dna'
                # This code filtrates requested dna type
                files_to_download = []
                for file in self.filenames:
                    if '.' + type + '.' in file:
                        files_to_download.append(file)

                os.mkdir(os.path.join(speciess, type))
                for filename in files_to_download:
                    print('download file: ' + filename)
                    file_path = os.path.join(speciess, type, filename)
                    fasta_path = os.path.join(speciess, type, filename[0:-3])
                    # Download file from Ensembl to local storage by FTP
                    with open(file_path, 'wb') as file:
                        ftp.retrbinary('RETR ' + filename, file.write, 102400)
                    # unziping
                    handle = gzip.open(file_path)
                    with open(fasta_path, 'wb') as out:
                        for line in handle:
                            out.write(line)
                    # parsing / analysis by sliding window method        
                    if self.parse:
                        self.data[speciess][type].update(self.parser.parse_file(fasta_path))
                ftp.quit()
            if self.jsons:
                if not os.path.exists(self.jsons_dir):
                  os.mkdir(self.jsons_dir)
                json_path = os.path.join(self.jsons_dir, speciess + '.json')
                print('creating ' + json_path)
                with open(json_path, 'w') as fp:
                    json.dump(self.data[speciess], fp)
            if self.remove_data:    
                print('removing ' + os.path.join(speciess))  
                shutil.rmtree(os.path.join(speciess))
            time.sleep(5)        
              
    # creating path for Ensembl FTP server based on animal kind and type of requested dna 
    def ftp_cwd(self, species, type):
        if 'dna' in type:
              return '/pub/release-101/fasta/' + species + '/' + 'dna' + '/'
        else:
              return '/pub/release-101/fasta/' + species + '/' + type + '/'
    
    def validate_name(self, name):
        for type in self.types:
            if ('.' + type + '.' in name) and ('.toplevel.' not in name) and ('.nonchromosomal.' not in name) and ('.MT.' not in name) and ('.abinitio.' not in name) and ('.alt.' not in name):
                self.filenames.append(name)
              
    def validate_name_toplevel(self, name):
        for type in self.types:
            if ('.' + type + '.' in name) and ('.toplevel.' in name) and ('.nonchromosomal.' not in name) and ('.MT.' not in name) and ('.abinitio.' not in name)and ('.alt.' not in name):
                self.filenames.append(name)          

    # Method for getting list of species for given group ie. Vertebrate
    # See: https://rest.ensembl.org/documentation/info/info_genomes_taxonomy
    # Example: dna_puller.get_species('vertebrate')
    # Result: JSON with response
    @staticmethod
    def get_species(klass_species):
        url_names = []
        server = "https://rest.ensembl.org"
        ext = "/info/genomes/taxonomy/"+ klass_species +"?"

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        for json_data in r.json():
            url_names.append(json_data['url_name'])

        return url_names
      

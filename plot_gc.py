from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import pandas as pd


def uppercase_ratio(sequence):
	"""
	Returns ration of upper-case letters (incl. N) in the sequence.

	>>> uppercase_ratio("AGCTNNN")
	1.0
	>>> uppercase_ratio("N")
	1.0
	>>> uppercase_ratio([])
	0.0
	>>> uppercase_ratio("agcN")
	0.25
	"""
	seq_len = len(sequence)
	if seq_len == 0:
		return 0.0

	upper = sum(1 for letter in sequence if letter.isupper())
	return upper / seq_len


class GCPlot:
	"""
	Example usage:

	>>> from glob import glob
	>>> gcplot = GCPlot()
	>>> gcplot.prepare_data(glob('Esox_lucius/dna_sm/*.fa'))
	>>> fig = gcplot.plot_all_in_one()
	>>> fig.savefig('Esox_lucius/gc.png')
	"""
	def __init__(self, window_size=1000, stride=1000):
		self.stride = stride
		self.window_size = window_size
		self.data = None

	def prepare_data(self, fasta_files):
		data = (rec for f in fasta_files for rec in SeqIO.parse(f, 'fasta'))
		data = ({'df': self._compute_features(lg.seq, lg.name), 'len': len(lg.seq), 'name': lg.name} for lg in data)
		data = sorted(data, key = lambda d: d['len'], reverse = True)
		self.data = data

	def _compute_features(self, sequence, debug_name):
		index = []
		current_gc = []
		current_uppercase_ratio = []
		seq_len = len(sequence)

		print('Computing', debug_name, seq_len)

		for start in range(0, seq_len, self.stride):
			end = min(seq_len, start + self.window_size)
			window = sequence[start:end]

			index.append(start)
			current_gc.append(GC(window))
			current_uppercase_ratio.append(uppercase_ratio(window) * 100)

		return pd.DataFrame({'gc': current_gc, 'upper': current_uppercase_ratio}, index = index)

	def plot_all_in_one(self, alpha=0.2):
		"""
		Plots GC % and uppercase letter % for all sequences in fasta files.

		:return: Figure
		"""

		width = self.data[0]['len']
		fig, axes_arr = plt.subplots(nrows=len(self.data), figsize=(20, 6 * len(self.data)), dpi=200)

		for d, ax in zip(self.data, axes_arr):
			df = d['df']
			print('Drawing', d['name'])
			scatter = ax.scatter(df.index, df['gc'], c=df['upper'], vmin=0, vmax=1, marker='.', cmap='jet', alpha=alpha)

			colorbar = fig.colorbar(scatter, ax=ax)
			colorbar.set_label('Uppercase letters %')
			colorbar.set_alpha(1)
			colorbar.draw_all()

			ax.set_title(d['name'])
			ax.set_xlabel('Postion, bp')
			ax.set_ylabel('GC %')
			ax.set_ylim(0, 100)
			ax.set_xlim(0, width)

		return fig

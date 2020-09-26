class GnuplotGenerator:

    def __init__(self, title, xrange):
        self.title = title
        self.xrange = xrange
        self.plots = []

    def add_palette(self, palette):
        self.palette = palette

    def add_plot(self, file, pattern, title):
        plot_attrs = {}

        plot_attrs['file'] = file
        plot_attrs['pattern'] = pattern
        plot_attrs['title'] = title

        self.plots.append(plot_attrs)

    def set_term(self, type, output_file, width, ratio):
        self.output_type = type
        self.output_file = output_file
        self.ratio = ratio
        self.width = width

    def prepare_definition(self):
        lines = []

        if self.palette:
            lines.append('set palette defined ' + self.palette)

        lines.append('set xrange ' + self.xrange)

        height = (self.width * self.ratio) * len(self.plots)

        lines.append('set term ' + self.output_type + ' size ' + str(self.width) + ', ' + str(int(height)))
        lines.append('set size ratio ' + str(self.ratio))
        lines.append('set output "' + self.output_file + '"')

        lines.append('set multiplot layout ' + str(len(self.plots)) + ',1 title "' + self.title + '" font ",26"')

        lines.append('set tmargin 0.2')
        lines.append('set bmargin 0.2')

        lines.append('set colorbox')
        lines.append('unset colorbox')

        for plot in self.plots:
            lines.append('set title "' + plot['title'] + '" font ",20"')
            lines.append("plot '" + plot['file'] + "' using " + plot['pattern'] + " w p pt 7 palette z lw 3")

        return lines





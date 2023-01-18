import os
import sys

from subprocess import call

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS, _TAGGERS

def write_line(file, line):
    file.write(line+'\n')

def escape_string_tex_style(string):
    return string.replace('_', '\_')

years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]
regions = [
'Pass',
'Fail',
]
channels = [
'muo',
'ele',
]
prepostfits = [
'prefitRaw',
# 'postfitCombine',
]
taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
'ak8_t_btagDCSV__tau',
'hotvr_t__tau',
'ak8_w__partnet',
'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}


class SlideCreator():

    def __init__(self):

        self.years = years
        self.regions = regions
        self.channels = channels

        self.tagger = taggers['hotvr_t__tau']
        # self.tagger = taggers['ak8_t__tau']
        # self.tagger = taggers['ak8_w__partnet']
        self.pt_bins = self.tagger.var_intervals
        # self.mscSplitting = 'mscTop3'
        self.mscSplitting = 'mscW3'

        self.tex_file_path = 'test_slides_{}.tex'.format(self.tagger.name)

    def create_tex_file(self):

        with open(self.tex_file_path, 'w') as tex_file:
            #__________________________________________________
            # latex header
            write_line(tex_file, r'''\documentclass[aspectratio=43,9pt]{beamer}''')
            write_line(tex_file, r'''\usetheme{CambridgeUS}''')
            write_line(tex_file, r'''\usecolortheme{seahorse}''')
            write_line(tex_file, r'''\usefonttheme{serif}''')
            write_line(tex_file, r'''\usepackage{makecell}''')
            write_line(tex_file, r'''\hypersetup{colorlinks,linkcolor=black,urlcolor=blue}''')
            write_line(tex_file, r'''\author{Christopher Matthies}''')
            write_line(tex_file, r'''\institute[U. Hamburg]{Institut f{\"u}r Experimentalphysik, Universit{\"a}t Hamburg}''')
            write_line(tex_file, r'''\title[Run 2 UL t jet tagging calibration]{Run 2 UL t jet tagging calibration}''')
            write_line(tex_file, r'''\date{Created \today}''')
            #__________________________________________________
            # begin document
            write_line(tex_file, r'''\begin{document}''')
            write_line(tex_file, r'''\setlength{\tabcolsep}{0.5em}''')

            # write_line(tex_file, r'''\begin{frame}{asdf}''')
            # write_line(tex_file, r'''jkl;''')
            # write_line(tex_file, r'''\end{frame}''')
            #__________________________________________________
            # iterate over slides for fit templates
            sorted_pt_bins = []
            for pt_bin in self.pt_bins.values():
                if not pt_bin.fit == True: continue
                sorted_pt_bins.append(pt_bin)
            # sorted_pt_bins = [x for x in self.pt_bins.values()]
            sorted_pt_bins = sorted(sorted_pt_bins, key=lambda x: x.var_min, reverse=False)
            for wp in self.tagger.get_wp():
                for year in self.years:
                    for prepostfit in prepostfits:
                        write_line(tex_file, r'''\begin{frame}[fragile]{\small{}'''+escape_string_tex_style(self.tagger.name)+r', WP: '+escape_string_tex_style(wp.name)+r''', \textcolor{black}{UL '''+_YEARS.get(year).get('year')+r''' \bgroup\small(MiniAODv2/NanoAODv9)\egroup}}''')
                        write_line(tex_file, r'''\begin{center}''')
                        write_line(tex_file, r'''\scalebox{0.7}{%''')
                        write_line(tex_file, r'''\begin{tabular}{c'''+r'c'*len(sorted_pt_bins)+r'''}''')
                        write_line(tex_file, r'''\textcolor{blue}{\bfseries{}'''+(r'pre-fit' if prepostfit.lower().startswith('prefit') else r'post-fit')+r'''} & \multicolumn{'''+str(len(sorted_pt_bins))+r'''}{c}{\textcolor{black}{\textit{lower} \quad{} $\longleftarrow$ \quad{} Probe jet $p_\text{T}$ [GeV] \quad{} $\longrightarrow$ \quad{} \textit{higher}}}\\[.2cm]''')
                        write_line(tex_file, r'''& '''+r'''& '''.join([pt_bin.get_interval_tex() for pt_bin in sorted_pt_bins])+r'\\[.2cm]')
                        for region in self.regions:
                            for channel in self.channels:
                                write_line(tex_file, r'''\makecell[b]{\textcolor{'''+(r'green' if region=='Pass' else r'red')+r'''}{'''+region.upper()+r'''} \\ '''+(r'\textmu{}' if channel=='muo' else r'e')+r'''+jets \\ ~\\} %''')
                                for pt_bin in sorted_pt_bins:
                                    # if not pt_bin.fit == True: continue
                                    # pdf_file_path = 'not_defined'
                                    substring = '-'.join([self.tagger.name, wp.name, pt_bin.name, year, self.tagger.fit_variable])
                                    region_string = '_'.join(['Main', region, channel])
                                    pdf_file_name = '-'.join(['Plot', prepostfit, substring, region_string, self.mscSplitting, 'withLeg'])+'.pdf'
                                    pdf_file_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', self.tagger.name, substring, 'plots', pdf_file_name)
                                    write_line(tex_file, r'''& \includegraphics[height=1.8cm]{{'''+pdf_file_path+'''}}''')
                                write_line(tex_file, r'''\\''')
                        write_line(tex_file, r'''\end{tabular}%''')
                        write_line(tex_file, r'''}''')
                        write_line(tex_file, r'''\end{center}''')
                        write_line(tex_file, r'''\end{frame}''')
            #__________________________________________________
            # end document
            write_line(tex_file, r'''\end{document}''')

    def compile_tex(self):

        # need to compile twice
        for i in range(0, 2):
            call("pdflatex {}".format(self.tex_file_path), shell=True)

if __name__=='__main__':

    sc = SlideCreator()
    sc.create_tex_file()
    sc.compile_tex()

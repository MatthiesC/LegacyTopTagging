from __future__ import print_function

import os
import sys

# from subprocess import call
import subprocess

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS, _TAGGERS, _NULL_WP

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
'prefitCombine',
'postfitCombine',
]
taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
# 'ak8_t_btagDCSV__tau',
# 'hotvr_t__tau',
# 'ak8_w__partnet',
'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}


class SlideCreator():

    def __init__(self, tagger_name, mscSplitting=None):

        self.years = years
        self.regions = regions
        self.channels = channels

        self.tagger = taggers[tagger_name]
        self.pt_bins = self.tagger.var_intervals
        self.pt_bin_total_range = None
        for pt_bin in self.pt_bins.values():
            if pt_bin.total_range == True:
                self.pt_bin_total_range = pt_bin
                break
        self.null_wp = _NULL_WP

        self.other_variables = []
        if self.tagger.name.startswith('ak8_t'):
            self.mscSplitting = 'mscTop3'
            self.other_variables = [
                'output_probejet_AK8_pt',
                'output_probejet_AK8_mass',
                'output_probejet_AK8_mSD',
                'output_probejet_AK8_tau32',
                'output_probejet_AK8_maxDeepJet',
                'output_probejet_AK8_maxDeepCSV',
                # 'output_probejet_AK8_MDDeepAK8_TvsQCD',
            ]
        elif self.tagger.name.startswith('ak8_w'):
            self.mscSplitting = 'mscW3'
            self.other_variables = [
                'output_probejet_AK8_pt',
                'output_probejet_AK8_mass',
                'output_probejet_AK8_mSD',
                'output_probejet_AK8_ParticleNet_WvsQCD',
            ]
        elif self.tagger.name.startswith('hotvr_t'):
            self.mscSplitting = 'mscTop3'
            self.other_variables = [
                'output_probejet_HOTVR_pt',
                'output_probejet_HOTVR_mass',
                'output_probejet_HOTVR_mpair',
                'output_probejet_HOTVR_fpt1',
                'output_probejet_HOTVR_nsub',
                'output_probejet_HOTVR_tau32',
            ]
        if 'MDdeepak8' in self.tagger.name:
            self.other_variables.append('output_probejet_AK8_MDDeepAK8_TvsQCD')
        if mscSplitting is not None:
            self.mscSplitting = mscSplitting

        self.tex_file_dir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'plots')
        os.system('mkdir -p '+self.tex_file_dir)
        self.tex_file_name = 'slides-{}-{}-{}.tex'.format(self.tagger.name, self.mscSplitting, '_'.join(prepostfits))
        self.tex_file_path = os.path.join(self.tex_file_dir, self.tex_file_name)

    def create_tex_file(self):

        with open(self.tex_file_path, 'w') as tex_file:
            #__________________________________________________
            # latex header
            write_line(tex_file, r'''\documentclass[aspectratio=43,9pt]{beamer}''')
            write_line(tex_file, r'''\usepackage{pgfpages}''')
            write_line(tex_file, r'''\pgfpagesuselayout{resize to}[physical paper height=288mm,physical paper width=384mm]''') # resize the slides to 3x higher zoom
            write_line(tex_file, r'''\usetheme{CambridgeUS}''')
            write_line(tex_file, r'''\usecolortheme{seahorse}''')
            write_line(tex_file, r'''\usefonttheme{serif}''')
            write_line(tex_file, r'''\usepackage{makecell}''')
            write_line(tex_file, r'''\usepackage{array}''')
            write_line(tex_file, r'''\hypersetup{colorlinks,linkcolor=black,urlcolor=blue}''')
            write_line(tex_file, r'''\author{Christopher Matthies}''')
            write_line(tex_file, r'''\institute[U. Hamburg]{Institut f{\"u}r Experimentalphysik, Universit{\"a}t Hamburg}''')
            write_line(tex_file, r'''\title[UHH2/LegacyTopTagging framework]{Run 2 UL jet tagging calibration\\ with the UHH2/LegacyTopTagging framework}''')
            write_line(tex_file, r'''\def\mydate{\leavevmode\hbox{\the\year-\twodigits\month-\twodigits\day}}''')
            write_line(tex_file, r'''\def\twodigits#1{\ifnum#1<10 0\fi\the#1}''')
            write_line(tex_file, r'''\date{Created with \ensuremath\heartsuit{} on \mydate}''')
            #__________________________________________________
            # begin document
            write_line(tex_file, r'''\begin{document}''')
            write_line(tex_file, r'''\setlength{\tabcolsep}{0.1em}''')
            #__________________________________________________
            # slide with substructure and other important variables
            if self.other_variables is not None and len(self.other_variables) > 0:
                for channel in self.channels:
                    write_line(tex_file, r'''\begin{frame}[fragile]{\small{}\textbf{Tagger: '''+escape_string_tex_style(self.tagger.name)+r'}, no WP applied \\ \textcolor{orange}{'''+(r'''\textmu{}''' if channel=='muo' else r'e')+r'''+jets}, \bgroup\small(Summer20 MiniAODv2/NanoAODv9)\egroup}''')
                    write_line(tex_file, r'''\begin{center}''')
                    write_line(tex_file, r'''\scalebox{0.7}{%''')
                    write_line(tex_file, r'''\begin{tabular}{c'''+r'c'*len(self.other_variables)+r'''}''')
                    write_line(tex_file, r'''\textcolor{blue}{\bfseries{} pre-fit} & \multicolumn{'''+str(len(self.other_variables))+r'''}{c}{\textcolor{black}{Overview of important probe jet variables}}\\[.2cm]''')
                    write_line(tex_file, r'''& \multicolumn{'''+str(len(self.other_variables))+r'''}{c}{Probe jet $p_\text{T}$: '''+self.pt_bin_total_range.get_interval_tex()+r''' GeV (\textit{Total range})}\\[.2cm]''')
                    for year in self.years:
                        write_line(tex_file, r'''\makecell[b]{\textcolor{orange}{'''+(year.split('p')[0] if 'UL16' in year else year)+r'''} \\ \textcolor{orange}{'''+(year.split('UL16')[-1] if 'UL16' in year else r'''~''')+r'''} \\ ~\\} %''')
                        for variable in self.other_variables:
                            substring = '-'.join([self.tagger.name, self.null_wp.name, self.pt_bin_total_range.name, year, variable])
                            region_string = '_'.join(['Main', 'Pass', channel])
                            pdf_file_name = '-'.join(['Plot', 'prefitRaw', substring, region_string, self.mscSplitting, 'withLeg'])+'.pdf'
                            pdf_file_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', self.tagger.name, substring, 'plots', pdf_file_name)
                            write_line(tex_file, r'''& \includegraphics[height=1.8cm]{{'''+pdf_file_path+r'''}}''')
                        write_line(tex_file, r'''\\''')
                    write_line(tex_file, r'''\end{tabular}%''')
                    write_line(tex_file, r'''}''')
                    write_line(tex_file, r'''\end{center}''')
                    write_line(tex_file, r'''\end{frame}''')
            #__________________________________________________
            # iterate over slides for fit templates
            sorted_pt_bins = []
            for pt_bin in self.pt_bins.values():
                if not pt_bin.fit == True: continue
                sorted_pt_bins.append(pt_bin)
            # sorted_pt_bins = [x for x in self.pt_bins.values()]
            sorted_pt_bins = sorted(sorted_pt_bins, key=lambda x: x.var_min, reverse=False)
            for wp in self.tagger.get_wp():
                write_line(tex_file, r'''\begin{frame}''')
                write_line(tex_file, r'''\vfill''')
                write_line(tex_file, r'''\centering''')
                write_line(tex_file, r'''\begin{beamercolorbox}[sep=8pt,center,shadow=false,rounded=true]{title}''')
                write_line(tex_file, r'''\usebeamerfont{title}\bfseries{}\vphantom{Ag}Tagger: '''+escape_string_tex_style(self.tagger.name)+r'''\\ \vphantom{Ag}WP: '''+escape_string_tex_style(wp.name)+r'''\par%''')
                write_line(tex_file, r'''\end{beamercolorbox}''')
                write_line(tex_file, r'''\vfill''')
                write_line(tex_file, r'''\end{frame}''')
                for year in self.years:
                    for prepostfit in prepostfits:
                        # if wp.name != 'BkgEff0p001' and prepostfit != 'prefitRaw': continue # HACK
                        write_line(tex_file, r'''\begin{frame}[fragile]{\small{}\textbf{Tagger: '''+escape_string_tex_style(self.tagger.name)+r', WP: '+escape_string_tex_style(wp.name)+r'''},\\ \textcolor{orange}{'''+year+r'''} \bgroup\small(Summer20 MiniAODv2/NanoAODv9)\egroup}''')
                        write_line(tex_file, r'''\begin{center}''')
                        write_line(tex_file, r'''\scalebox{0.7}{%''')
                        write_line(tex_file, r'''\begin{tabular}{c'''+r'c'*(len(sorted_pt_bins))+r'''|c}''')
                        # write_line(tex_file, r'''\textcolor{blue}{\bfseries{}'''+(r'pre-fit' if prepostfit.lower().startswith('prefit') else r'post-fit')+r'''} & \multicolumn{'''+str(len(sorted_pt_bins))+r'''}{c|}{\textcolor{black}{\textit{lower} \quad{} $\longleftarrow$ \quad{} Probe jet $p_\text{T}$ [GeV] \quad{} $\longrightarrow$ \quad{} \textit{higher}}} & \textit{Total range}\\[.2cm]''')
                        if prepostfit.lower().startswith('postfit'):
                            write_line(tex_file, r'''\textcolor{red}{\bfseries{}\makebox[2cm][c]{postfit}}''')
                        else:
                            write_line(tex_file, r'''\textcolor{blue}{\bfseries{}\makebox[2cm][c]{prefit}}''')
                        write_line(tex_file, r'''& \multicolumn{'''+str(len(sorted_pt_bins))+r'''}{c|}{\textcolor{black}{\textit{lower} \quad{} $\longleftarrow$ \quad{} Probe jet $p_\text{T}$ [GeV] \quad{} $\longrightarrow$ \quad{} \textit{higher}}} & \textit{Total range}\\[.2cm]''')
                        write_line(tex_file, r'''& '''+r'''& '''.join([pt_bin.get_interval_tex() for pt_bin in sorted_pt_bins]+[self.pt_bin_total_range.get_interval_tex()])+r'\\[.2cm]')
                        for region in self.regions:
                            for channel in self.channels:
                                write_line(tex_file, r'''\makecell[b]{\textcolor{'''+(r'green' if region=='Pass' else r'red')+r'''}{'''+region.upper()+r'''} \\ '''+(r'\textmu{}' if channel=='muo' else r'e')+r'''+jets \\ ~\\} %''')
                                for pt_bin in sorted_pt_bins+[self.pt_bin_total_range]:
                                    # if not pt_bin.fit == True: continue
                                    # pdf_file_path = 'not_defined'
                                    substring = '-'.join([self.tagger.name, wp.name, pt_bin.name, year, self.tagger.fit_variable])
                                    region_string = '_'.join(['Main', region, channel])
                                    pdf_file_name = '-'.join(['Plot', prepostfit, substring, region_string, self.mscSplitting, 'withLeg'])+'.pdf'
                                    pdf_file_found = False
                                    if prepostfit == 'postfitCombine' or prepostfit == 'prefitCombine':
                                        # combine_task_name_suffix = ''

                                        # combine_task_name_suffix = '-'.join([('PtTotal' if pt_bin.total_range else 'PtSplit'), 'Run2']) # original
                                        combine_task_name_suffix = '-'.join([('PtTotal' if pt_bin.total_range else pt_bin.name), year]) # HACK

                                        combine_task_name = '-'.join(['combineTask', self.tagger.name, wp.name]) + ('-'+combine_task_name_suffix if len(combine_task_name_suffix) else '')
                                        pdf_file_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'workdirs', combine_task_name, 'plots_templates', pdf_file_name)
                                        pdf_file_found = os.path.isfile(pdf_file_path)
                                        if not pdf_file_found:
                                            print('Plot not found:', pdf_file_path, '\nFallback to prefitRaw plot.')
                                    if prepostfit == 'prefitRaw' or (not pdf_file_found):
                                        pdf_file_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', self.tagger.name, substring, 'plots', pdf_file_name)
                                        pdf_file_path = pdf_file_path.replace('prefitCombine', 'prefitRaw').replace('postfitCombine', 'prefitRaw') # the replace is a bit HACKy
                                        # print(pdf_file_path)
                                        pdf_file_found = os.path.isfile(pdf_file_path)
                                    if not pdf_file_found:
                                        print('Plot not found:', pdf_file_path, '\nPrefer to exit.')
                                        sys.exit()
                                    write_line(tex_file, r'''& \includegraphics[height=1.8cm]{{'''+pdf_file_path+r'''}}''')
                                write_line(tex_file, r'''\\''')
                        write_line(tex_file, r'''\end{tabular}%''')
                        write_line(tex_file, r'''}''')
                        write_line(tex_file, r'''\end{center}''')
                        write_line(tex_file, r'''\end{frame}''')
            #__________________________________________________
            # end document
            write_line(tex_file, r'''\end{document}''')

        print('Created', self.tex_file_path)

    def compile_tex(self):

        # need to compile twice
        for i in range(0, 2):
            # call("pdflatex {}".format(self.tex_file_path), shell=True)
            command = "pdflatex {}".format(self.tex_file_path)
            p = subprocess.Popen((command), shell=True, cwd=self.tex_file_dir)
            p.wait()
        print('Created', self.tex_file_path.replace('.tex', '.pdf'))

if __name__=='__main__':

    # # sc = SlideCreator('hotvr_t__tau')
    # sc = SlideCreator('ak8_t__tau')
    # sc.create_tex_file()
    # sc.compile_tex()

    for tagger in taggers.values():
        sc = SlideCreator(tagger.name)
        sc.create_tex_file()
        sc.compile_tex()

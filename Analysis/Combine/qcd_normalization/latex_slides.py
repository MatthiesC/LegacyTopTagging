import os
import sys


inDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/{year}/{channel}/qcd_normalization_study/plots')
inPath = os.path.join(inDir, '{filename}')

outName = 'qcd_normalization_study.pdf'
outDir = 'latex'
os.system('mkdir -p '+outDir)
outPath = os.path.join(outDir, outName)

with open(outPath, 'w') as outfile:
    outfile.write("""
\documentclass[9pt]{beamer]\n
\n
\begin{document}\n
\n
\begin{frame}\n
\begin{figure}\n
\begin{tabular}{cccc}\n
UL16preVFP & UL16postVFP & UL17 & UL18 \\\n
    """)
# \includegraphics{}
    outfile.write("""
\end{tabular}
\end{figure}
\end{frame}

\end{document}
    """)

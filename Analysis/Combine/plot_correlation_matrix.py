import ROOT as root
import os
import sys
from array import array
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json

workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine/')
workdir = os.path.join(workdir, 'ak8_t__tau/workdirs/combineTask-ak8_t__tau-BkgEff0p001-pt_480to600-UL16preVFP')
plot_name = 'correlation_matrix_nuisances.pdf'
plot_path = os.path.join(workdir, plot_name)
rootfile_name = 'robustHesse_obs.root'
rootfile_path = os.path.join(workdir, rootfile_name)

rootfile = root.TFile.Open(rootfile_path, 'READ')
matrix = rootfile.Get('h_correlation') #TH2F

correlation_matrix = np.zeros((matrix.GetNbinsX(), matrix.GetNbinsY()))
for i in range(matrix.GetNbinsX()):
    for j in range(matrix.GetNbinsY()):
        correlation_matrix[i][j] = matrix.GetBinContent(i+1, j+1)

# Get the bin labels from the histogram
x_labels = [matrix.GetXaxis().GetBinLabel(i) for i in range(1, matrix.GetNbinsX() + 1)]
y_labels = [matrix.GetYaxis().GetBinLabel(i) for i in range(1, matrix.GetNbinsY() + 1)]


# Define the label to remove (replace with the actual label you want to remove)
label_to_remove = "Scale_QCD"

# Find the index of the label to remove in x_labels and y_labels
if label_to_remove in x_labels:
    index_to_remove = x_labels.index(label_to_remove)
    x_labels.pop(index_to_remove)
    y_labels.pop(index_to_remove)
    correlation_matrix = np.delete(correlation_matrix, index_to_remove, axis=0)
    correlation_matrix = np.delete(correlation_matrix, index_to_remove, axis=1)
else:
    print("Label not found in x_labels")


# Load the JSON mapping file for labels
with open("rename_GoodExample-ak8_t__tau-BkgEff0p001-UL16preVFP-pt_480to600_latex.json", "r") as json_file:
    label_mapping = json.load(json_file)
# Replace x and y labels using the JSON mapping
x_labels = [label_mapping.get(label, label) for label in x_labels]
y_labels = [label_mapping.get(label, label) for label in y_labels]
# Enable LaTeX rendering for Matplotlib
plt.rcParams["text.usetex"] = True


print(correlation_matrix)
print(x_labels)
print(y_labels)






# Create a Seaborn heatmap with labels
plt.figure(figsize=(10, 6))
heatmap = sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap="coolwarm", xticklabels=x_labels, yticklabels=y_labels, vmin=-1, vmax=1, annot_kws={"size": 4}, cbar_kws={'ticks': [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]})

# Adjust x and y tick label font size
plt.xticks(fontsize=6)  # Adjust font size here
plt.yticks(fontsize=6)  # Adjust font size here

# # Set z-axis ticks in steps of 0.2 from -1 to +1
# z_ticks = np.arange(-1, 1.2, 0.2)
# heatmap.set_zticks(np.linspace(0, len(z_ticks) - 1, len(z_ticks)))
# heatmap.set_zticklabels(z_ticks)

plt.title("Postfit nuisance parameter correlations")
# plt.xlabel("", fontsize=6)  # Adjust x-axis label size here
# plt.ylabel("", fontsize=6)  # Adjust y-axis label size here

# Rotate x-axis labels for better readability (optional)
plt.xticks(rotation=45, ha='right')

# # Add custom text at the top left of the canvas
# plt.figtext(0.05, 0.95, r"\textbf{CMS} \textit{Private Work}", fontsize=14)
# # plt.figtext(0.05, 0.91, r"Private Work", fontsize=12, fontstyle="italic")

# Add custom text at the top left of the heatmap
# plt.annotate(r"CMS", xy=(0.01, 0.95), xycoords='axes fraction', fontsize=14, fontweight="bold")
# plt.annotate(r"Private Work", xy=(0.01, 0.91), xycoords='axes fraction', fontsize=12, fontstyle="italic")

# plt.annotate(r"\textbf{CMS} \textit{Private Work}", xy=(0.01, 0.95), xycoords='axes fraction', fontsize=14)

# # Add custom text at the top left of the canvas
# plt.figtext(0.01, 0.99, r"$\mathbf{CMS}$ Private Work", fontsize=14, fontweight="bold", verticalalignment="top")
# plt.figtext(0.01, 0.96, r"Private Work", fontsize=12, fontstyle="italic", verticalalignment="top")



# Save the heatmap plot as an image file (e.g., PNG)
plt.savefig(plot_path, bbox_inches="tight")

# Close the plot
plt.close()



# canvas = root.TCanvas('canvas', '', 1800, 1200)
# canvas.cd()
# canvas.SetTopMargin(0.01)
# canvas.SetLeftMargin(0.1)
# canvas.SetRightMargin(0.1)
# canvas.SetBottomMargin(0.3)
# matrix.GetXaxis().SetTickSize(0)
# matrix.GetYaxis().SetTickSize(0)
# canvas.SetTickx(0)
# canvas.SetTicky(0)
# canvas.SetLineWidth(5)
# #histogram settings
# # Red = array('d', [0.00, 1.00, 0.90])
# # Green = array('d', [0.90, 1.00, 0.20])
# # Blue = array('d', [ 0.00, 1.00, 0.20])
# Red = array('d', [0.00, 1.00, 0.90])
# Green = array('d', [0.00, 1.00, 0.20])
# Blue = array('d', [ 0.90, 1.00, 0.20])
# Length = array('d', [0.00, 0.50, 1.00])
# root.TColor.CreateGradientColorTable(3,Length,Red,Green,Blue,50)
#
# # root.gStyle.SetPalette(87)
# root.gStyle.SetOptStat(0)
#
# matrix.SetTitle('')
# matrix.GetZaxis().SetTitle('Correlation coefficient')
# matrix.GetZaxis().SetLabelSize(0.02)
# matrix.GetZaxis().SetLabelSize(0.02)
# matrix.SetMaximum(1.)
# matrix.SetMinimum(-1.)
# matrix.SetMarkerSize(0.5) # this adjusts the text size of the text drawn with the "text" draw option
# root.gStyle.SetPaintTextFormat('.2f')
# matrix.Draw('colz text')
# canvas.SaveAs(plot_path)

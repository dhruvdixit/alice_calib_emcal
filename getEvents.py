#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.environ['ROOTSYS'], 'lib'))
import ROOT
from ROOT import TH1D

import bad_channel_lhc15o_mine3
import bad_channel_lhc15o_mineallBut2


filename_list = []
for f in sys.argv[1:]:
    filename_list.append(f)

ncell = 17664

cell_ampHist = []
for filename in filename_list:
    f = ROOT.TFile.Open(filename)
    h = f.Get('AliAnalysisTaskCalibEmcal').Get('histogram').FindObject('_histogram_cell_id_amplitude')
    cell_ampHist.append(ROOT.TH2D(h))

cell_hists = []
for hist in cell_ampHist:
    root_histogram_spectrum = {}
    for i in range(ncell):
        name = 'cell'+str(i)
        root_histogram_spectrum[i] = ROOT.TH1D(hist.ProjectionY(name, i + 1, i + 1))
    cell_hists.append(root_histogram_spectrum)

events = []
#h = cell_hists[0]
#h = 
for hist in cell_hists:
    esum = 0
    for key in hist.keys():#key here is a dictionary of cell_id, energy_hist
        #if not (key in bad_channel_lhc15o_mineallBut2.hot):
        #if not (key in bad_channel_lhc15o_mine3.hot):
            h1 = hist[key]
            esum += int(h1.GetEntries())
    events.append(esum)

print esum

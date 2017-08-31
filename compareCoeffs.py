#! /usr/bin/env python
from os import path
import sys
import ROOT
from ROOT import TCanvas, TFile, TH1D, TH2D, TLine, TH1F, TColor, TStyle
from ROOT import gROOT, gSystem, gStyle


#gStyle.SetOptStat(0)

compare_file = 'allNew'

if compare_file == 'allNew':
    import bad_channel_lhc15o_all_new
    import bad_channel_lhc15o_mine_all_new
    bad_all = sorted(bad_channel_lhc15o_all_new.bad_all + bad_channel_lhc15o_mine_all_new.hot)
    filename = 'lhc15o_all_new_chiBeforeAfter.txt' 
    xaxisTitle = 'all_new'


if compare_file == 'allOld':
    import bad_channel_lhc15o_all
    import bad_channel_lhc15o_mineAll
    bad_all = sorted(bad_channel_lhc15o_all.bad_all + bad_channel_lhc15o_mineAll.hot)
    filename = 'lhc15o_all_chiBeforeAfter.txt' 
    xaxisTitle = 'all_old'

if compare_file == 'allBut2':
    import bad_channel_lhc15o_allBut2
    import bad_channel_lhc15o_mineallBut2
    bad_all = sorted(bad_channel_lhc15o_allBut2.bad_all + bad_channel_lhc15o_mineallBut2.hot)
    filename = 'lhc15o_allBut2_chiBeforeAfter.txt' 
    xaxisTitle = 'allBut2'

if compare_file == 'list3':
    import bad_channel_lhc15o_3
    import bad_channel_lhc15o_mine31
    bad_all = sorted(bad_channel_lhc15o_3.bad_all + bad_channel_lhc15o_mine3.hot)
    filename = 'lhc15o_list3_chiBeforeAfter.txt'
    xaxisTitle = 'c3'
    
nCells = 17664
nCells_emcal = 12288

#make some histogram
fhCellCoeff = TH1F("fhCellCoeff","fhCellCoeff",nCells+10,0,nCells+10)# cellID+1 = histogram bin

h2 = TH1D('h2', 'Coeff plot for list 2', 400, -0.5, 3.5)
h3 = TH1D('h3', 'Coeff plot for '+compare_file, 400, -0.5, 3.5)

cB = TH1D('cB', 'Chi Square before the correction for ' + xaxisTitle, 16100, -0.5, 16095.5)
cB.GetXaxis().SetTitle('#chi^{2}')
cA = TH1D('cA', 'Chi Square after the correction' + xaxisTitle, 890, -0.5, 889.5)
cA.GetXaxis().SetTitle('#chi^{2}')

rat = TH1D('rat', 'Differece between lists 2 and '+ xaxisTitle, 300, -1.5, 1.5)
rat.GetXaxis().SetTitle('c2-'+ xaxisTitle)

corrG = TH2D('corrG', 'Correlation between the two lists: good channels only', 400, -0.5, 3.5, 400, -0.5, 3.5)
corrG.GetXaxis().SetTitle('c2')
corrG.GetYaxis().SetTitle(xaxisTitle)
corrG_EMCAL = TH2D('corrG_EMCAL', 'Correlation between the two lists: good channels EMCAL only', 200, 0.5, 2.5, 200, 0.5, 2.5)
corrG_EMCAL.GetXaxis().SetTitle('c2')
corrG_EMCAL.GetYaxis().SetTitle(xaxisTitle)
corrG_DCAL = TH2D('corrG_DCAL', 'Correlation between the two lists: good channels DCAL only', 200, 0.5, 2.5, 200, 0.5, 2.5)
corrG_DCAL.GetXaxis().SetTitle('c2')
corrG_DCAL.GetYaxis().SetTitle(xaxisTitle)
corrB = TH2D('corrB', 'Correlation between the two lists: bad channels only', 400, -0.5, 3.5, 400, -0.5, 3.5)



list2_dict = {}
list3_dict = {}
with open('lhc15o_list2_chiBeforeAfter.txt','r') as list2_txt:
    for line2 in list2_txt:
        cell_id2, eb2, yb2, coeff2, chiA2, chiB2 = [float(i) for i in line2.split()]
        cell_id2 = int(cell_id2)
        eb2 = int(eb2)
        yb2 = int(yb2)
        list2 = [eb2, yb2, coeff2, chiA2, chiB2]
        list2_dict[cell_id2] = list2
        h2.Fill(coeff2)
        #print cell_id2, list2
        #if cell_id2 == 20:
        #    break

#print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

with open(filename, 'r') as list3_txt:
    for line3 in list3_txt:      
        cell_id3, eb3, yb3, coeff3, chiA3, chiB3 = [float(j) for j in line3.split()]
        cell_id3 = int(cell_id3)
        eb3 = int(eb3)
        yb3 = int(yb3)
        list3 = [eb3, yb3, coeff3, chiA3, chiB3]
        list3_dict[cell_id3] = list3
        h3.Fill(coeff3)
        #print cell_id3, list3
        #if cell_id3 == 20:
        #    break

    
bad_channels = []
missing_channels = []
outliers = []#1391, 2163, 2350, 2790, 4416, 6049, 6519, 7507, 8283, 9003, 9892, 10474, 11867, 13967, 14209, 14240, 14241, 14242, 14243, 14247, 14248, 14250, 14251, 14252, 14253, 14288, 14291, 14292, 14293, 14294, 14295, 14296, 14297, 14298, 14299, 14301, 14381, 14386, 14390, 14392, 16686, 17305, 17445, 17460, 17464, 17527]
chiA_min = 9999.9
chiB_min = 9999.9
chiA_max = 0.0
chiB_max = 0.0
for key in list2_dict.keys():
    if key in list3_dict:
        l2 = list2_dict[key]
        l3 = list3_dict[key]
        c2 = float(l2[2])
        c3 = float(l3[2])
        if key in bad_all:
            corrB.Fill(c2, c3)
            bad_channels.append(key)
            #print key, 'is a bad channel'
        else:
            if False:
                if(chiA_min > float(l3[3])):
                    chiA_min = float(l3[3])
                if(chiB_min > float(l3[4])):
                    chiB_min = float(l3[4])
                if(chiA_max < float(l3[3])):
                    chiA_max = float(l3[3])
                if(chiB_max < float(l3[4])):
                    chiB_max = float(l3[4])

            if c2 < 1.5 and c3 < 1.5:
                if key < nCells_emcal:
                    corrG_EMCAL.Fill(c2, c3)
                else:
                    corrG_DCAL.Fill(c2, c3)

                cA.Fill(float(l3[3]))
                cB.Fill(float(l3[4]))

            else:
                outliers.append(key)

                    
            if abs(c2-c3) <= 0.2:
                corrG.Fill(c2, c3)
                rat.Fill(c2-c3)
    else:
        #print key, 'is missing list3_dict'
        missing_channels.append(key)


print outliers, len(outliers)

cOutliers = TH1D('cOutliers', 'Chi Square of the outliers after the correction' + xaxisTitle, 890, -0.5, 889.5)
for n in outliers:
    cOutliers.Fill(list3_dict[n][3])


#print chiA_min, chiB_min, chiA_max, chiB_max

# #print bad_channels, '\n\n'
# print missing_channels

# for cell in range(nCells):
#     coeff = 0
#     if(cell in list3_dict):
#         coeff = (list3_dict[cell])[2]
#     fhCellCoeff.SetBinContent(cell+1,coeff)

# cvs = TCanvas('cb', 'cb', 600, 400)
# fhCellCoeff.Draw()
# eFile = TFile("cellCoeffHist.root", "RECREATE")
# if eFile.IsOpen :
#    fhCellCoeff.Write("fhCellCoeff")
# eFile.Close()


chiCanvas = TCanvas("chiCanvas", "chiCanvas", 600, 800)
chiCanvas.Divide(2,1)
chiCanvas.cd(1)
cB.Draw()
chiCanvas.cd(2)
cA.Draw()

EDCanvas = TCanvas("EDCanvas", "EDCanvas", 2100, 700)
EDCanvas.Divide(2,1)
EDCanvas.cd(1)
EDCanvas.SetLogz()
corrG_EMCAL.Draw('colz')
EDCanvas.cd(2)
EDCanvas.SetLogz()
corrG_DCAL.Draw('colz')


# lCorr = TLine(0,0,3.5,3.5)
cg = TCanvas('cg', 'cg', 600, 400)
cg.SetLogz()
corrG.Draw('colz')
# lnG = TLine(0,1,3.5,1)
# lnG.Draw()
# lCorr.Draw()

# cb = TCanvas('cb', 'cb', 600, 400)
# cb.SetLogz()
# corrB.Draw('colz')
# lnB = TLine(1,0,1,3)
# lnB.Draw()
# lCorr.Draw()

cr = TCanvas('cr', 'cr', 600, 400)
rat.Draw()
rat.Fit('gaus','','',-0.2,0.2)

cO = TCanvas('cO', 'cO', 600, 400)
cOutliers.Draw()


coeff_cvns = TCanvas('coeff_cvns','coeff_cvns',600,400)
#h2.SetLineColor(kRed)
#h3.SetLineColor(kBlue)
#h2.Draw()
h3.Draw()

#myFile = TFile("chiBnAForAll.root", "RECREATE")
#if myFile.IsOpen :
#    corrG.Write("corrG")
#    corrG_EMCAL.Write("corrG_EMCAL")
#    corrG_DCAL.Write("corrG_DCAL")
#    corrB.Write("corrB")
#    rat.Write("rat")
#    cA.Write("cA")
#    cB.Write("cB")
#myFile.Close()



# ---------------------------------------------------------------------------------- 
raw_input('Press Enter to exit')

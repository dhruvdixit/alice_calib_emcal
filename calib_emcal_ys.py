#!/usr/bin/env python

def factorized_model(hcount, hspectrum, bad_channel, gamma):
    ncell = hcount.GetNbinsX()
    content = [hcount.GetBinContent(i + 1) for i in range(ncell)]
    neta_max = 48
    nphi_max = 24
    nsm = 20
    mean_ieta = [[0 for i in range(neta_max)] for j in range(nsm)]
    count_ieta = [[0 for i in range(neta_max)] for j in range(nsm)]
    mean_iphi_2 = [[0 for i in range(2 * nphi_max)] for j in range(nsm)]
    count_iphi_2 = [[0 for i in range(2 * nphi_max)] for j in range(nsm)]
    mean_sm = [0 for j in range(nsm)]
    count_sm = [0 for j in range(nsm)]
    for i in range(ncell):
        if content[i] > 0 and i not in bad_channel:
            sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
            eta, phi = paint_emcal.to_eta_phi(sm, ieta, iphi)
            mean_ieta[sm][ieta] += content[i]
            count_ieta[sm][ieta] += 1
            iphi_2 = iphi + nphi_max * (ieta % 2)
            mean_iphi_2[sm][iphi_2] += content[i]
            count_iphi_2[sm][iphi_2] += 1
            mean_sm[sm] += content[i]
            count_sm[sm] += 1
    for sm in range(nsm):
        for ieta in range(neta_max):
            if count_ieta[sm][ieta] > 0:
                mean_ieta[sm][ieta] /= count_ieta[sm][ieta]
        for iphi_2 in range(2 * nphi_max):
            if count_iphi_2[sm][iphi_2] > 0:
                mean_iphi_2[sm][iphi_2] /= count_iphi_2[sm][iphi_2]
        if count_sm[sm] > 0:
            mean_sm[sm] /= count_sm[sm]
    print_warm = True
    function_inner = []
    for i in plot_cell_id + range(ncell):
        sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
        eta, phi = paint_emcal.to_eta_phi(sm, ieta, iphi)
        # n(ieta, iphi) ~ m(ieta) * m(iphi) / m
        content_model = 0
        if mean_sm[sm] > 0:
            content_model = mean_ieta[sm][ieta] * \
                            mean_iphi_2[sm][iphi + nphi_max * (ieta % 2)] / \
                            mean_sm[sm]
        # if i not in bad_channel and content[i] > 0 and content_model > 0:
        #     print (content_model / content[i])**(1.0 / 5.0), ','
        if content[i] > 0:
            a = 0
            if print_warm:
                if (not i in bad_channel_lhc15o_5.bad_all) and (content[i] - content_model) > 1000:
                    sys.stdout.write('%d, ' % i)
            elif content_model > 0:
                a = (content_model / content[i])**(-1.0 / gamma)
                if a != 0:
                    # scale_15_17, scale_17_20, scale_20_50 = [((a * r[1])**(gamma + 1) / (gamma + 1) - (a * r[0])**(gamma + 1) / (gamma + 1)) / (r[1]**(gamma + 1) / (gamma + 1) - r[0]**(gamma + 1) / (gamma + 1)) for r in ((1.5, 17), (1.7, 2.0), (2.0, 5.0))]
                    chi_square = []
                    if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                        histogram_model = [ROOT.TH1D(hspectrum[i]) for j in range(2)]
                        for h in histogram_model:
                            h.Reset()
                    for ia1 in range(2):
                        if ia1 == 0:
                            a1 = a
                        else:
                            a1 = 1
                        s = 0
                        dof = 0
                        for j in range(
                                hspectrum[i].GetXaxis().FindBin(1.5 / a1),
                                hspectrum[i].GetXaxis().FindBin(5 / a1) + 1):
                            content_bin = hspectrum[i].GetBinContent(j)
                            # content_model_bin = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5 / a1) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5 / a1))**(gamma + 1) / (gamma + 1))
                            # content_model_bin_orig = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5))**(gamma + 1) / (gamma + 1))
                            content_model_bin = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5 / a1) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5 / a1))**(gamma + 1) / (gamma + 1))
                            content_model_bin_orig = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / ((5 / a1)**(gamma + 1) / (gamma + 1) - (1.5 / a1)**(gamma + 1) / (gamma + 1))
                            if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                                if ia1 == 0:
                                    histogram_model[0].SetBinContent(j, content_model_bin)
                                else:
                                    histogram_model[1].SetBinContent(j, content_model_bin)
                            if content_bin > 0:
                                # if len(plot_cell_id) >= 1:
                                #     print i, ia1, a1, j, hspectrum[i].GetBinCenter(j), content_bin, content_model_bin, content_model_bin_orig
                                s += (content_bin -
                                               content_model_bin)**2 / \
                                    content_bin
                                dof += 1
                        if dof > 0:
                            s /= dof
                        chi_square.append(s)
                    if chi_square[0] / chi_square[1] < 1 and i in bad_channel_lhc15o_5.bad_all:
                        print i, i in bad_channel_lhc15o_5.bad_all and 1 or 0, i in bad_channel_lhc15o_mine.hot and 1 or 0, a, chi_square[0], chi_square[1], chi_square[0] / chi_square[1] #, content[i], content_model
                        #a = chi_square[0] / chi_square[1]
                        if a < 1.2 and chi_square[0] < 10:
                            a = 1
                        else:
                            a = 2
                    else:
                        a = 0
                        None
                    #a = chi_square[1]
                    if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                        canvas.cd()
                        canvas.SetLogx()
                        canvas.SetLogy()
                        hspectrum[i].SetMarkerStyle(20)
                        histogram_model[0].SetMarkerStyle(24)
                        histogram_model[1].SetMarkerStyle(25)
                        hspectrum[i].GetXaxis().SetTitle('E_T\\:(\\mathrm{GeV})')
                        hspectrum[i].GetYaxis().SetTitle('Counts')
                        hspectrum[i].GetYaxis().SetTitleFont(42)
                        hspectrum[i].GetYaxis().SetLabelFont(42)
                        hspectrum[i].SetMinimum(1)
                        hspectrum[i].Draw('e1x0')
                        histogram_model[0].Draw('e1x0same')
                        histogram_model[1].Draw('e1x0same')
                        canvas.Update()
                        ROOT.gApplication.Run()
        #content[i] = a
        content[i] = content[i] - content_model
        #content[i] -= content_model
        #if not (i in bad_channel_lhc15o_5.bad_all or i in bad_channel_lhc15o_mine.hot):
            #content[i] = 0
        content[i] = max(-2200, min(2200, content[i]))
    if print_warm:
        sys.stdout.write('\n')
    for i in range(ncell):
        hcount.SetBinContent(i + 1, content[i])
    return hcount

if __name__ == '__main__':
    import os, sys
    sys.path.append(os.path.join(os.environ['ROOTSYS'], 'lib'))
    import math, array, re, ROOT, paint_emcal
    application_name = ' '.join(sys.argv)
    paint_emcal.set_root_style()
    log_z = False
    filename_list = []
    plot_cell_id = []
    import_run_group = None
    for f in sys.argv[1:]:
        if f in ('-l', '--log-z'):
            log_z = True
        elif re.compile('^[0-9]+$').match(f):
            plot_cell_id.append(int(f))
        elif re.compile('^r[0-9]+$').match(f):
            import_run_group = int(f[1:])
        else:
            filename_list.append(f)
    # f, bad_channel_map = sys.argv[1:3]
    filename_list = [f]
    
    if len(plot_cell_id) > 0:
        canvas = canvas = ROOT.TCanvas(
            'canvas%d' % 0, application_name,
            800 + 4, 600 + 28)
        pad = []
    else:
        canvas, pad = paint_emcal.alice_emcal_canvas_pad(application_name)
    ncell = 17664

    # for i in range(ncell):
    #     sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
    #     if sm == 9 and ieta == 33 and iphi >= 3 and iphi <= 7:
    #         print i
    # sys.exit(0)

    root_histogram = []
    for i in range(1):
        root_histogram.append(ROOT.TH1D(
            'root_histogram%d' % len(root_histogram), '',
            ncell, -0.5, ncell + 0.5))
    hsum = None
    for filename in filename_list:
        f = ROOT.TFile.Open(filename)
        h = f.Get('AliAnalysisTaskCalibEmcal').Get(
            'histogram').FindObject('_histogram_cell_id_amplitude')
        if hsum == None:
            hsum = ROOT.TH2D(h)
        else:
            hsum.Add(h)
        hpx = ROOT.TH1D(h.ProjectionX(
            'hpx',
            h.GetYaxis().FindBin(1.5),
            h.GetYaxis().FindBin(5) - 1))
        print h.GetYaxis().FindBin(1.7), h.GetYaxis().FindBin(2)
        hpx_15_17 = ROOT.TH1D(h.ProjectionX(
            'hpx',
            h.GetYaxis().FindBin(1.5),
            h.GetYaxis().FindBin(1.7) - 1))
        hpx_17_20 = ROOT.TH1D(h.ProjectionX(
            'hpx',
            h.GetYaxis().FindBin(1.7),
            h.GetYaxis().FindBin(2) - 1))
        hpx_20_50 = ROOT.TH1D(h.ProjectionX(
            'hpx',
            h.GetYaxis().FindBin(2),
            h.GetYaxis().FindBin(5) - 1))
        #hpx.Divide(hpx_17_20, hpx_20_50)
        root_histogram[-1] = ROOT.TH1D(hpx)
    # exec('import ' + re.compile('\\.py').sub('', bad_channel_filename))
    # bad_channel_filename[:-3]
    import bad_channel_lhc15o_5y
    import bad_channel_lhc15o_mine
    bad_channel_lhc15o_mine.hot = [871, 1160, 1162, 1163, 1382, 1682, 1683, 1685, 1686, 1760, 1851, 1961, 2071, 2273, 3544, 3691, 3849, 3858, 3878, 3884, 4019, 4025, 4047, 4051, 4089, 4245, 4566, 4599, 4605, 5410, 5648, 5650, 5652, 5836, 6111, 6296, 6701, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6809, 7089, 7104, 7323, 7874, 8811, 8812, 8850, 8852, 8856, 8857, 8903, 8905, 8906, 8909, 8955, 8956, 8959, 8992, 8994, 8996, 9798, 9807, 10823, 11088, 11089, 11101, 11142, 11144, 11146, 11187, 11188, 11393, 11424, 11525, 11871, 11937, 12032, 12607, 12778, 12783, 12801, 12814, 13250, 13253, 13281, 13394, 13400, 13439, 13463, 13847, 13861, 13863, 13865, 14320, 14326, 14363, 14396, 14406, 14410, 14411, 14413, 14501, 14554, 14559, 14592, 14594, 14609, 14611, 14640, 14661, 14663, 14704, 14712, 14713, 14730, 14735, 14737, 14765, 14991, 15143, 15152, 15276, 15281, 15316, 15317, 15348, 15365, 15373, 15375, 15870, 15908, 15936, 15938, 15979, 15983, 16297, 16307, 16340, 16346, 16485, 16502, 16506, 16523, 16768, 16852, 16858, 16910, 16916, 16988, 17003, 17070, 17097, 17214, 17287, 17300, 17311, 17413, 17442, 17457, 17461, 17470, 17577]
    #bad_channel_lhc15o_mine.hot = [1061, 1276, 1682, 2778, 3135, 3544, 4026, 4047, 4605, 5608, 6340, 6701, 7089, 7874, 8812, 8904, 8907, 8911, 8944, 9338, 9798, 9801, 9807, 10091, 10138, 10823, 10829, 10831, 10981, 10988, 11034, 11088, 11095, 11099, 11102, 11141, 11144, 11146, 11148, 11150, 11197, 11411, 11630, 11647, 12360, 12607, 12614, 12783, 12848, 12856, 12874, 12916, 12919, 12988, 13122, 13156, 13248, 13249, 13260, 13263, 13280, 13281, 13288, 13293, 13294, 13394, 13403, 13433, 13439, 13467, 13509, 13510, 13574, 13578, 13786, 13847, 13903, 13921, 13926, 13933, 14058, 14060, 14157, 14306, 14311, 14325, 14362, 14373, 14501, 14503, 14505, 14554, 14559, 14605, 14617, 14620, 14622, 14640, 14659, 14660, 14662, 14663, 14668, 14672, 14699, 14704, 14712, 14731, 14733, 14735, 14738, 14741, 14754, 14767, 14773, 14839, 14872, 14977, 14991, 15009, 15015, 15020, 15027, 15033, 15094, 15137, 15143, 15164, 15282, 15293, 15318, 15326, 15333, 15348, 15350, 15398, 15476, 15670, 15949, 15969, 15972, 15983, 16026, 16388, 16422, 16428, 16440, 16443, 16446, 16468, 16481, 16492, 16509, 16768, 16819, 16821, 16826, 16830, 16858, 16907, 17097, 17107, 17208, 17396, 17452, 17463, 17531] + [1061, 1160, 1163, 1275, 1276, 1382, 1682, 1683, 1685, 1686, 1760, 1961, 2778, 3135, 3485, 3544, 3849, 3858, 3878, 3884, 3955, 4026, 4047, 4531, 4566, 4605, 5608, 5652, 5836, 6111, 6340, 6701, 6800, 6801, 6802, 6804, 6806, 6807, 6809, 7089, 7104, 7323, 7874, 8811, 8812, 8850, 8855, 8857, 8904, 8906, 8907, 8909, 8911, 8944, 8947, 8956, 9338, 9798, 9801, 9807, 10091, 10138, 10823, 10829, 10831, 10981, 10988, 11034, 11048, 11088, 11093, 11095, 11097, 11099, 11101, 11102, 11141, 11142, 11144, 11146, 11148, 11150, 11197, 11411, 11630, 11647, 11871, 11937, 12063, 12276, 12360, 12607, 12614, 12783, 12814, 12818, 12826, 12848, 12856, 12874, 12916, 12919, 12922, 12988, 13122, 13156, 13248, 13249, 13253, 13260, 13263, 13280, 13281, 13286, 13288, 13290, 13293, 13294, 13394, 13400, 13403, 13433, 13439, 13467, 13508, 13509, 13510, 13574, 13578, 13786, 13847, 13861, 13903, 13921, 13926, 13933, 13983, 14058, 14060, 14117, 14125, 14157, 14306, 14311, 14314, 14325, 14362, 14373, 14396, 14409, 14413, 14501, 14503, 14505, 14554, 14559, 14605, 14613, 14617, 14620, 14622, 14640, 14659, 14660, 14662, 14663, 14668, 14672, 14699, 14703, 14704, 14712, 14731, 14733, 14735, 14737, 14738, 14741, 14742, 14754, 14767, 14773, 14774, 14839, 14872, 14876, 14977, 14991, 15009, 15015, 15020, 15027, 15033, 15094, 15137, 15143, 15164, 15280, 15282, 15293, 15318, 15326, 15333, 15348, 15350, 15398, 15476, 15670, 15949, 15969, 15972, 15979, 15983, 16026, 16089, 16340, 16346, 16388, 16422, 16428, 16440, 16443, 16446, 16468, 16472, 16475, 16481, 16485, 16492, 16497, 16502, 16509, 16768, 16819, 16821, 16826, 16830, 16858, 16907, 16910, 17097, 17107, 17208, 17287, 17396, 17442, 17452, 17457, 17463, 17531, 17537]
    bad_all = sorted(bad_channel_lhc15o_5.bad_all + bad_channel_lhc15o_mine.hot)
    function = []
    hsumpy = None
    root_histogram_spectrum = []
    for i in range(ncell):
        root_histogram_spectrum.append(ROOT.TH1D(hsum.ProjectionY('e', i + 1, i + 1)))
        if not i in bad_all:
            if hsumpy == None:
                hsumpy = ROOT.TH1D(root_histogram_spectrum[-1])
            else:
                hsumpy.Add(root_histogram_spectrum[-1])
    function.append(ROOT.TF1(
        'function%d' % len(function), 'exp([0])*x^[1]', 1.5, 5))
    function[-1].SetParameter(0, 1)
    function[-1].SetParameter(1, -4)
    for i in range(1, hsumpy.GetNbinsX()):
        hsumpy.SetBinContent(i, hsumpy.GetBinContent(i) / hsumpy.GetXaxis().GetBinWidth(i))
        hsumpy.SetBinError(i, hsumpy.GetBinError(i) / hsumpy.GetXaxis().GetBinWidth(i))
    hsumpy.Fit(function[-1], 'r0')
    # for i in range(hsumpy.GetXaxis().FindBin(1.5),
    #                hsumpy.GetXaxis().FindBin(5)):
    #     print i, hsumpy.GetXaxis().GetBinCenter(i), hsumpy.GetBinContent(i)

    #bad_all = bad_channel_lhc15o_1.bad_all
    # if pad[0] == None: sys.exit(1)
    root_histogram[-1] = factorized_model(root_histogram[-1], root_histogram_spectrum, bad_all, function[-1].GetParameter(1))
    canvas.cd()
    for i in range(len(pad)):
        pad[i].Draw()
        if log_z:
            pad[i].SetLogz()
    paint_emcal.update(canvas, pad, root_histogram)
    ROOT.gApplication.Run()

import ROOT
import uproot
import numpy as np
ROOT.gROOT.SetBatch(True)

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

#set gRandom seed
ROOT.gRandom.SetSeed(1995)


def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for i in range(df.shape[0]):
        h.Fill(df[var1].iloc[i], df[var2].iloc[i])


def fill_res_hist(h, df, var1, var2):
    for var_val1, var_val2 in zip(df[var1], df[var2]): 
        h.Fill((var_val1 - var_val2)/var_val1)

def fill_res_hist_th2(h, df, var1, var2):
    for var_val1, var_val2 in zip(df[var1], df[var2]): 
        h.Fill(var_val1,(var_val1 - var_val2)/var_val1)


def fill_mass_weighted_hist(h, df, var, weight = [1,1]):

    for var_val, w in zip(df[var], df['isSignal']):
        if w == 1:
            h.Fill(var_val, weight[0])
        else:
            h.Fill(var_val, weight[1])

def significance_error(signal, background, signal_error, background_error):

    sb = signal + background + 1e-10
    sb_sqrt = np.sqrt(sb)

    s_propag = (sb_sqrt + signal / (2 * sb_sqrt))/sb * signal_error
    b_propag = signal / (2 * sb_sqrt)/sb * background_error

    if signal+background == 0:
        return 0

    return np.sqrt(s_propag * s_propag + b_propag * b_propag)

def scale_hist_content(h, scale):
    ## generate poissonian counts
    for i in range(1, h.GetNbinsX()+1):
        pois = ROOT.gRandom.Poisson(scale)
        pois_sqrt = np.sqrt(pois)
        h.SetBinContent(i, h.GetBinContent(i)+pois)
        h.SetBinError(i, np.sqrt(pois_sqrt*pois_sqrt + h.GetBinError(i)*h.GetBinError(i)))


def set_style():
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptDate(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetLabelSize(0.04,"xyz")
    ROOT.gStyle.SetTitleSize(0.05,"xyz") 
    ROOT.gStyle.SetTitleFont(42,"xyz")
    ROOT.gStyle.SetLabelFont(42,"xyz")
    ROOT.gStyle.SetTitleOffset(1.05,"x")
    ROOT.gStyle.SetTitleOffset(1.1,"y")
    ROOT.gStyle.SetCanvasDefW(800)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetPadBottomMargin(0.12) 
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadGridX(0) 
    ROOT.gStyle.SetPadGridY(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetPaperSize(20,24)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetEndErrorSize(0.)
    ROOT.gStyle.SetMarkerSize(1)




def fit_and_plot(histo, var, fit_function, signal, background, sigma, mu, n):

    data = ROOT.RooDataHist('data', 'data', ROOT.RooArgList(var), histo)
    fit_results = fit_function.fitTo(data, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True))
    output_file = ROOT.TFile('fit_results.root', 'recreate')
    frame = var.frame()
    frame.SetName('frame_tree')
    frame.SetTitle('')
    set_style()
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetMaxDigits(2)
    frame.GetXaxis().SetTitleOffset(1.1)

    data.plotOn(frame, ROOT.RooFit.Name('data'))
    fit_function.plotOn(frame, ROOT.RooFit.LineColor(kBlueC), ROOT.RooFit.Name('fit_func'))
    fit_function.plotOn(frame, ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(kOrangeC))
    # fit_function.plotOn(frame, ROOT.RooFit.Components('cb'), ROOT.RooFit.LineStyle(ROOT.kDashed))

    frame.SetMinimum(histo.GetMinimum()*0.95)
    frame.SetMaximum(histo.GetMaximum()*1.05)

    sigma_val = sigma.getVal()
    mu_val = mu.getVal()

    signal_counts = n.getVal()*histo.Integral()
    signal_counts_error = (n.getError()/n.getVal())*n.getVal()*histo.Integral()


    background_counts = (1-n.getVal())*histo.Integral()
    background_counts_error = (1-n.getVal())*histo.Integral()*n.getError()/n.getVal()

    #signal within 3 sigma
    var.setRange('signal', mu_val-3*sigma_val, mu_val+3*sigma_val)
    signal_int = signal.createIntegral(ROOT.RooArgSet(var), ROOT.RooArgSet(var), 'signal')
    signal_int_val_3s = signal_int.getVal()*signal_counts
    signal_int_val_3s_error = signal_int_val_3s*signal_counts_error/signal_counts
    #background within 3 sigma
    var.setRange('bkg', mu_val-3*sigma_val, mu_val+3*sigma_val)
    bkg_int = background.createIntegral(ROOT.RooArgSet(var), ROOT.RooArgSet(var), 'bkg')
    bkg_int_val_3s = bkg_int.getVal()*background_counts
    bkg_int_val_3s_error = bkg_int_val_3s*background_counts_error/background_counts
    significance = signal_int_val_3s/np.sqrt(signal_int_val_3s + bkg_int_val_3s)
    significance_err = significance_error(signal_int_val_3s, bkg_int_val_3s, signal_int_val_3s_error, bkg_int_val_3s_error)
    s_b_ratio_err = np.sqrt((signal_int_val_3s_error/signal_int_val_3s)**2 + (bkg_int_val_3s_error/bkg_int_val_3s)**2)*signal_int_val_3s/bkg_int_val_3s


    # rootchi2 = ROOT.RooChi2Var("chi2", "chi2", fit_function, data, ROOT.RooFit.DataError(ROOT.RooAbsData.Expected))
    # chi2 = rootchi2.getVal()
    # fit_probability = ROOT.TMath.Prob(chi2, frame.GetNbinsX())

    chi2 = frame.chiSquare("fit_func", "data", 6)
    fit_probability = ROOT.TMath.Prob(chi2*(frame.GetNbinsX() -6), frame.GetNbinsX() -6)

    print('chi2: ', chi2)
    print('fit probability: ', fit_probability)
    print('muv: ', mu_val)

    # print('chi2: ', chi2, 'frame.GetNbinsX(): ', frame.GetNbinsX())
    # print('fit probability: ', fit_probability)

    pinfo = ROOT.TPaveText(0.632, 0.5, 0.932, 0.85, 'NDC')
    pinfo.SetBorderSize(0)
    pinfo.SetFillStyle(0)
    pinfo.SetTextAlign(11)
    pinfo.SetTextFont(42)
    string_list = []

    string_list.append('3 < #eta < 5, p_{T} > 0.5 GeV/#it{c}')
    string_list.append('Fit Probability: ' + f'{fit_probability:.2f}')
    string_list.append(f'Signal (S): {signal_counts:.0f} #pm {signal_counts_error:.0f}')
    string_list.append(f'S/B (3 #sigma): {signal_int_val_3s/bkg_int_val_3s:.2f} #pm {s_b_ratio_err:.2f}')
    string_list.append('S/#sqrt{S+B} (3 #sigma): ' + f'{significance:.0f} #pm {significance_err:.0f}')
    string_list.append('#mu = ' + f'{mu_val*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
    for s in string_list:
        pinfo.AddText(s)

    string_list = []
    pinfo2 = ROOT.TPaveText(0.17, 0.6, 0.45, 0.85, "NDC")
    pinfo2.SetBorderSize(0)
    pinfo2.SetFillStyle(0)
    pinfo2.SetTextAlign(11)
    pinfo2.SetTextFont(42)

    string_list.append("NA60+ Performance")
    string_list.append("Pb#minusPb, #sqrt{#it{s}_{NN}} = 6.06 GeV")
    string_list.append('{}^{5}_{#Lambda}He #rightarrow ^{4}He+p+#pi^{-}')
    for s in string_list:
        pinfo2.AddText(s)


    frame.addObject(pinfo)
    frame.addObject(pinfo2)

    fit_stats = {"signal": [signal_counts, signal_counts_error],
    "significance": [significance, significance_err], "s_b_ratio": [signal_int_val_3s/bkg_int_val_3s, s_b_ratio_err]}

    frame.Write()

    return fit_stats, fit_results
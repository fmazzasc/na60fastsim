import ROOT
import uproot
import numpy as np
from utils import *
ROOT.gROOT.SetBatch(True)

filter_data = False
fit = True



if filter_data:

    histDecLRes = ROOT.TH1F('histDecL', 'Decay length', 100, -1, 1)
    histPtRes = ROOT.TH1F('histPt', 'Pt', 100, -1, 1)
    th2PtRes = ROOT.TH2D('th2Pt', ';#it{p}_{T}^{gen};(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}', 100, 0, 2, 100, -1, 1)
    th2PzRes = ROOT.TH2D('th2Pz', ';#it{p}_z^{gen};(#it{p}_z^{gen} - #it{p}_z^{rec})/#it{p}_z^{gen}', 100, 30, 60, 100, -1, 1)
    th2DecLRes = ROOT.TH2F('th2DecL', ';L^{gen};(L^{gen} - L^{rec})/L^{gen}', 100, 0, 7, 100, -1, 1)
    histCtGen = ROOT.TH1F('ct gen', '; ct; Counts', 100, 0, 40)
    histCtRec = ROOT.TH1F('ct rec', '; ct; Counts', 100, 0, 40)
    histPGen = ROOT.TH1F('p gen', '; p; Counts', 200, 0, 60)
    histPRec = ROOT.TH1F('p rec', '; p; Counts', 200, 0, 60)
    histMRec = ROOT.TH1F('m_rec', '; m; Counts', 200, 4.8, 4.88)
    histMBkg = ROOT.TH1F('m_bkg', '; m; Counts', 200, 4.8, 4.88)


    hypermass = 4.839961
    df = uproot.open('/data/fmazzasc/na60plus/Hypernuclei-Signal-ntuple.root')['hyper;13'].arrays(library='pd')
    df['pxMC'] = df['ptMC'] * np.cos(df['phiMC'])
    df['pyMC'] = df['ptMC'] * np.sin(df['phiMC'])
    df['pzMC'] = df['ptMC'] * np.sinh(df['etaMC'])
    df['pMC'] = np.sqrt(df['pxMC']**2 + df['pyMC']**2 + df['pzMC']**2)
    df['decLMC'] = np.sqrt(df['xMC']**2 + df['yMC']**2 + df['zMC']**2)
    df.eval('ctMC=decLMC*4.839961/(pMC)', inplace=True)
    df_rec = df.query('reconstructed==3')
    df_rec.reset_index(drop=True, inplace=True)
    df_rec['px'] = df_rec['pt'] * np.cos(df_rec['phi'])
    df_rec['py'] = df_rec['pt'] * np.sin(df_rec['phi'])
    df_rec['pz'] = df_rec['pt'] * np.sinh(df_rec['eta'])
    df_rec['p'] = np.sqrt(df_rec['px']**2 + df_rec['py']**2 + df_rec['pz']**2)
    df_rec['decL'] = np.sqrt(df_rec['x']**2 + df_rec['y']**2 + df_rec['z']**2)
    # df_rec['dcaNucl'] = np.sqrt(df_rec['dcaX']**2 + df_rec['dcaY']**2 )
    df_rec.eval('ct=decL*4.839961/p', inplace=True)
    df_rec['isSignal'] = 1



    df_bkg = uproot.open('/data/fmazzasc/na60plus/Hypernuclei-Background-ntuple.root')['hyper'].arrays(library='pd')
    df_bkg['px'] = df_bkg['pt'] * np.cos(df_bkg['phi'])
    df_bkg['py'] = df_bkg['pt'] * np.sin(df_bkg['phi'])
    df_bkg['pz'] = df_bkg['pt'] * np.sinh(df_bkg['eta'])
    df_bkg['p'] = np.sqrt(df_bkg['px']**2 + df_bkg['py']**2 + df_bkg['pz']**2)
    df_bkg['decL'] = np.sqrt(df_bkg['x']**2 + df_bkg['y']**2 + df_bkg['z']**2)
    # df_bkg['dcaNucl'] = np.sqrt(df_bkg['dcaX']**2 + df_bkg['dcaY']**2 )
    df_bkg.eval('ct=decL*4.839961/p', inplace=True)
    df_bkg['isSignal'] = 0
    df_full = df_rec.append(df_bkg, ignore_index=True)
    df_full.query('pt>0.5 and 3<eta<5 and z>0 and cosPA>0.99999', inplace=True)



    fill_th1_hist(histCtGen, df, 'ctMC')
    fill_th1_hist(histCtRec, df_rec, 'ct')
    fill_th1_hist(histPGen, df, 'pMC')
    fill_th1_hist(histPRec, df_rec, 'p')
    fill_mass_weighted_hist(histMRec, df_full, 'm', [0.34, 1])
    fill_th1_hist(histMBkg, df_full.query('isSignal==0'), 'm')
    fill_res_hist(histDecLRes, df_rec, 'decLMC', 'decL')
    fill_res_hist(histPtRes, df_rec, 'ptMC', 'pt')
    fill_res_hist_th2(th2PtRes, df_rec, 'ptMC', 'pt')
    fill_res_hist_th2(th2DecLRes, df_rec, 'decLMC', 'decL')
    fill_res_hist_th2(th2PzRes, df_rec, 'pzMC', 'pz')


    histMBkg.Fit('pol0', 'L', '', 4.82, 4.9)
    pol0_res = histMBkg.GetFunction('pol0').GetParameter(0)
    scale = pol0_res*(20e3)
    scale_hist_content(histMRec, scale)

    ffile = ROOT.TFile('res_histos.root', 'recreate')
    histCtGen.Write()
    histCtRec.Write()
    histPGen.Write()
    histPRec.Write()
    histMRec.Write()
    histMBkg.Write()
    histDecLRes.Write()
    histPtRes.Write()
    th2PtRes.Write()
    th2DecLRes.Write()
    th2PzRes.Write()
    ffile.Close()



if fit:
    ROOT.gROOT.LoadMacro('RooCustomPdfs/RooDSCBShape.cxx++')
    from ROOT import RooDSCBShape
    kBlueC = ROOT.TColor.GetColor('#1f78b4')
    kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

    mass = ROOT.RooRealVar('m', '#it{M}_{^{4}He+p+#pi^{-}}', 4.8, 4.9, 'GeV/c^{2}')
    mu = ROOT.RooRealVar('mu', 'hypernucl mass', 4.83, 4.85, 'GeV/c^{2}')
    sigma = ROOT.RooRealVar('sigma', 'hypernucl width', 0.001, 0.004, 'GeV/c^{2}')
    a1 = ROOT.RooRealVar('a1', 'a1', 0, 3.)
    a2 = ROOT.RooRealVar('a2', 'a2', 0.3, 0.8)
    n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
    n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
    signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)
    c0 = ROOT.RooRealVar('c0', 'constant c0', 0., 1e5)
    background = ROOT.RooPolynomial('bkg', 'pol1 bkg', mass, ROOT.RooArgList())
    n = ROOT.RooRealVar('n', 'n const', 0.01, 0.1)
    # define the fit funciton and perform the actual fit
    fit_function = ROOT.RooAddPdf('total_pdf', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(n))


    input_file = ROOT.TFile('res_histos.root', 'read')
    histMRec = input_file.Get('m_rec')
    histMRec.SetDirectory(0)
    input_file.Close()


    output_file = ROOT.TFile('fit_results.root', 'recreate')
    fitstat, fitres = fit_and_plot(histMRec, mass, fit_function, signal, background, sigma, mu, n)





    










# print(df_full.columns)


# print(df[['decLMC','ctMC', 'pMC']][0:50])
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <vector>
#include "RooPolynomial.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFormulaVar.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include "TColor.h"
#include "TAxis.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TLegend.h>
#include <TObject.h>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*);
void plot(TH1F* data_p_bn_var_name,RooRealVar r_yield_bn,RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,TString comp,RooAddPdf totalPdf_var_name,RooDataSet* data_set,RooDataSet* r_data_set,RooDataHist p_bn_var_name,TCanvas* TCanvas_var_name,TLegend* leg_var_name);
void weight_plot(RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,RooDataSet* data_set,TCanvas* TCanvas_var_name,TLegend* leg_var_name);

void splots()
{
    RooWorkspace* wspace = new RooWorkspace("mywspace");
   // add the signal and background models to the workspace.
    AddModel(wspace);
}


void AddModel(RooWorkspace* ws)
{

    TFile* file1 = new TFile("TMVA_output_splots.root");
    TTree* datatree = (TTree*)(file1->Get("tree"));

    TFile* r_file1 = new TFile("r_TMVA_output_splots.root");
    TTree* r_datatree = (TTree*)(r_file1->Get("tree"));
  
    TH1F* data_p_b_bdt = (TH1F*)(file1->Get("p_b_bdt"));
    TH1F* data_p_n_bdt = (TH1F*)(file1->Get("p_n_bdt"));
    TH1F* data_p_bdt = (TH1F*)(file1->Get("p_bdt"));
    TH1F* r_data_p_bdt = (TH1F*)(r_file1->Get("r_p_bdt"));

    TH1F* r_data_p_30_bdt = (TH1F*)(r_file1->Get("r_p_30_bdt"));
    TH1F* data_p_b_30_bdt = (TH1F*)(file1->Get("p_b_30_bdt"));
    TH1F* data_p_n_30_bdt = (TH1F*)(file1->Get("p_n_30_bdt"));
    TH1F* r_data_p_40_bdt = (TH1F*)(r_file1->Get("r_p_40_bdt"));
    TH1F* data_p_b_40_bdt = (TH1F*)(file1->Get("p_b_40_bdt"));
    TH1F* data_p_n_40_bdt = (TH1F*)(file1->Get("p_n_40_bdt"));  
    TH1F* data_p_b_60_bdt = (TH1F*)(file1->Get("p_b_60_bdt"));
    TH1F* data_p_n_60_bdt = (TH1F*)(file1->Get("p_n_60_bdt"));
    TH1F* data_p_60_bdt = (TH1F*)(file1->Get("p_60_bdt"));
    TH1F* r_data_p_60_bdt = (TH1F*)(r_file1->Get("r_p_60_bdt"));
    TH1F* data_p_b_80_bdt = (TH1F*)(file1->Get("p_b_80_bdt"));
    TH1F* data_p_n_80_bdt = (TH1F*)(file1->Get("p_n_80_bdt"));
    TH1F* data_p_80_bdt = (TH1F*)(file1->Get("p_80_bdt"));
    TH1F* r_data_p_80_bdt = (TH1F*)(r_file1->Get("r_p_80_bdt"));
    TH1F* r_data_p_100_bdt = (TH1F*)(r_file1->Get("r_p_100_bdt"));
    TH1F* data_p_b_100_bdt = (TH1F*)(file1->Get("p_b_100_bdt"));
    TH1F* data_p_n_100_bdt = (TH1F*)(file1->Get("p_n_100_bdt"));
    TH1F* data_p_b_120_bdt = (TH1F*)(file1->Get("p_b_120_bdt"));
    TH1F* data_p_n_120_bdt = (TH1F*)(file1->Get("p_n_120_bdt"));
    TH1F* data_p_120_bdt = (TH1F*)(file1->Get("p_120_bdt"));
    TH1F* r_data_p_120_bdt = (TH1F*)(r_file1->Get("r_p_120_bdt"));

    TH1F* data_p_b_60_120_bdt = (TH1F*)(data_p_b_60_bdt->Clone("data_p_b_60_120_bdt"));
    data_p_b_60_120_bdt->Add(data_p_b_80_bdt);
    data_p_b_60_120_bdt->Add(data_p_b_100_bdt);
    data_p_b_60_120_bdt->Add(data_p_b_120_bdt); 

    TH1F* data_p_n_60_120_bdt = (TH1F*)(data_p_n_60_bdt->Clone("data_p_n_60_120_bdt"));
    data_p_n_60_120_bdt->Add(data_p_n_80_bdt);
    data_p_n_60_120_bdt->Add(data_p_n_100_bdt);
    data_p_n_60_120_bdt->Add(data_p_n_120_bdt);  

    TH1F* r_data_p_60_120_bdt = (TH1F*)(r_data_p_60_bdt->Clone("r_data_p_60_120_bdt"));
    r_data_p_60_120_bdt->Add(r_data_p_80_bdt);
    r_data_p_60_120_bdt->Add(r_data_p_120_bdt);
    r_data_p_60_120_bdt->Add(r_data_p_120_bdt);  

//    TH1F* data_p_bdt = (TH1F*)data_p_b_bdt->Clone("p_bdt");
//    data_p_bdt->Add(data_p_n_bdt);

    TH1F* data_p_b_Jet_leptonDeltaR = (TH1F*)(file1->Get("p_b_Jet_leptonDeltaR"));
    TH1F* data_p_b_Jet_mass = (TH1F*)(file1->Get("p_b_Jet_mass"));
    TH1F* data_p_b_Jet_pt = (TH1F*)(file1->Get("p_b_Jet_pt"));
    TH1F* data_p_b_Jet_eta = (TH1F*)(file1->Get("p_b_Jet_eta"));
    TH1F* data_p_b_Jet_phi = (TH1F*)(file1->Get("p_b_Jet_phi"));
    TH1F* data_p_b_Jet_energyRing_dR0_em = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR0_em"));
    TH1F* data_p_b_Jet_energyRing_dR1_em = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR1_em"));
    TH1F* data_p_b_Jet_energyRing_dR2_em = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR2_em"));
    TH1F* data_p_b_Jet_energyRing_dR3_em = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR3_em"));
    TH1F* data_p_b_Jet_energyRing_dR4_em = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR4_em"));
    TH1F* data_p_b_Jet_energyRing_dR0_mu = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR0_mu"));
    TH1F* data_p_b_Jet_energyRing_dR1_mu = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR1_mu"));
    TH1F* data_p_b_Jet_energyRing_dR2_mu = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR2_mu"));
    TH1F* data_p_b_Jet_energyRing_dR3_mu = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR3_mu"));
    TH1F* data_p_b_Jet_energyRing_dR4_mu = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR4_mu"));
    TH1F* data_p_b_Jet_energyRing_dR0_ch = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR0_ch"));
    TH1F* data_p_b_Jet_energyRing_dR1_ch = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR1_ch"));
    TH1F* data_p_b_Jet_energyRing_dR2_ch = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR2_ch"));
    TH1F* data_p_b_Jet_energyRing_dR3_ch = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR3_ch"));
    TH1F* data_p_b_Jet_energyRing_dR4_ch = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR4_ch"));
    TH1F* data_p_b_Jet_energyRing_dR0_neut = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR0_neut"));
    TH1F* data_p_b_Jet_energyRing_dR1_neut = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR1_neut"));
    TH1F* data_p_b_Jet_energyRing_dR2_neut = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR2_neut"));
    TH1F* data_p_b_Jet_energyRing_dR3_neut = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR3_neut"));
    TH1F* data_p_b_Jet_energyRing_dR4_neut = (TH1F*)(file1->Get("p_b_Jet_energyRing_dR4_neut"));
    TH1F* data_p_b_Jet_rawEnergy = (TH1F*)(file1->Get("p_b_Jet_rawEnergy"));
    TH1F* data_p_b_Jet_numDaughters_pt03 = (TH1F*)(file1->Get("p_b_Jet_numDaughters_pt03"));
    TH1F* data_p_b_Jet_numberOfDaughters = (TH1F*)(file1->Get("p_b_Jet_numberOfDaughters"));
    TH1F* data_p_b_Jet_rawPt = (TH1F*)(file1->Get("p_b_Jet_rawPt"));
    TH1F* data_p_b_Jet_chHEF = (TH1F*)(file1->Get("p_b_Jet_chHEF"));
    TH1F* data_p_b_Jet_muEF = (TH1F*)(file1->Get("p_b_Jet_muEF"));
    TH1F* data_p_b_Jet_neHEF = (TH1F*)(file1->Get("p_b_Jet_neHEF"));
    TH1F* data_p_b_Jet_chEmEF = (TH1F*)(file1->Get("p_b_Jet_chEmEF"));
    TH1F* data_p_b_Jet_neEmEF = (TH1F*)(file1->Get("p_b_Jet_neEmEF"));
    TH1F* data_p_b_Jet_leadTrackPt = (TH1F*)(file1->Get("p_b_Jet_leadTrackPt"));
    TH1F* data_p_b_Jet_leptonPt = (TH1F*)(file1->Get("p_b_Jet_leptonPt"));
    TH1F* data_p_b_Jet_leptonPtRel = (TH1F*)(file1->Get("p_b_Jet_leptonPtRel"));
    TH1F* data_p_b_Jet_leptonPdgId = (TH1F*)(file1->Get("p_b_Jet_leptonPdgId"));
    TH1F* data_p_b_Jet_leptonPtRelInv = (TH1F*)(file1->Get("p_b_Jet_leptonPtRelInv"));
    TH1F* data_p_b_Jet_vtxMass = (TH1F*)(file1->Get("p_b_Jet_vtxMass"));
    TH1F* data_p_b_Jet_vtxNtracks = (TH1F*)(file1->Get("p_b_Jet_vtxNtracks"));
    TH1F* data_p_b_Jet_vtxPt = (TH1F*)(file1->Get("p_b_Jet_vtxPt"));
    TH1F* data_p_b_Jet_vtx3DSig = (TH1F*)(file1->Get("p_b_Jet_vtx3DSig"));
    TH1F* data_p_b_Jet_vtx3DVal = (TH1F*)(file1->Get("p_b_Jet_vtx3DVal"));
    TH1F* data_p_b_Jet_ptd = (TH1F*)(file1->Get("p_b_Jet_ptd"));
    TH1F* data_p_b_rho = (TH1F*)(file1->Get("p_b_rho"));

    TH1F* data_p_n_Jet_leptonDeltaR = (TH1F*)(file1->Get("p_n_Jet_leptonDeltaR"));
    TH1F* data_p_n_Jet_mass = (TH1F*)(file1->Get("p_n_Jet_mass"));
    TH1F* data_p_n_Jet_pt = (TH1F*)(file1->Get("p_n_Jet_pt"));
    TH1F* data_p_n_Jet_eta = (TH1F*)(file1->Get("p_n_Jet_eta"));
    TH1F* data_p_n_Jet_phi = (TH1F*)(file1->Get("p_n_Jet_phi"));
    TH1F* data_p_n_Jet_energyRing_dR0_em = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR0_em"));
    TH1F* data_p_n_Jet_energyRing_dR1_em = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR1_em"));
    TH1F* data_p_n_Jet_energyRing_dR2_em = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR2_em"));
    TH1F* data_p_n_Jet_energyRing_dR3_em = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR3_em"));
    TH1F* data_p_n_Jet_energyRing_dR4_em = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR4_em"));
    TH1F* data_p_n_Jet_energyRing_dR0_mu = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR0_mu"));
    TH1F* data_p_n_Jet_energyRing_dR1_mu = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR1_mu"));
    TH1F* data_p_n_Jet_energyRing_dR2_mu = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR2_mu"));
    TH1F* data_p_n_Jet_energyRing_dR3_mu = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR3_mu"));
    TH1F* data_p_n_Jet_energyRing_dR4_mu = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR4_mu"));
    TH1F* data_p_n_Jet_energyRing_dR0_ch = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR0_ch"));
    TH1F* data_p_n_Jet_energyRing_dR1_ch = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR1_ch"));
    TH1F* data_p_n_Jet_energyRing_dR2_ch = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR2_ch"));
    TH1F* data_p_n_Jet_energyRing_dR3_ch = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR3_ch"));
    TH1F* data_p_n_Jet_energyRing_dR4_ch = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR4_ch"));
    TH1F* data_p_n_Jet_energyRing_dR0_neut = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR0_neut"));
    TH1F* data_p_n_Jet_energyRing_dR1_neut = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR1_neut"));
    TH1F* data_p_n_Jet_energyRing_dR2_neut = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR2_neut"));
    TH1F* data_p_n_Jet_energyRing_dR3_neut = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR3_neut"));
    TH1F* data_p_n_Jet_energyRing_dR4_neut = (TH1F*)(file1->Get("p_n_Jet_energyRing_dR4_neut"));
    TH1F* data_p_n_Jet_rawEnergy = (TH1F*)(file1->Get("p_n_Jet_rawEnergy"));
    TH1F* data_p_n_Jet_numDaughters_pt03 = (TH1F*)(file1->Get("p_n_Jet_numDaughters_pt03"));
    TH1F* data_p_n_Jet_numberOfDaughters = (TH1F*)(file1->Get("p_n_Jet_numberOfDaughters"));
    TH1F* data_p_n_Jet_rawPt = (TH1F*)(file1->Get("p_n_Jet_rawPt"));
    TH1F* data_p_n_Jet_chHEF = (TH1F*)(file1->Get("p_n_Jet_chHEF"));
    TH1F* data_p_n_Jet_muEF = (TH1F*)(file1->Get("p_n_Jet_muEF"));
    TH1F* data_p_n_Jet_neHEF = (TH1F*)(file1->Get("p_n_Jet_neHEF"));
    TH1F* data_p_n_Jet_chEmEF = (TH1F*)(file1->Get("p_n_Jet_chEmEF"));
    TH1F* data_p_n_Jet_neEmEF = (TH1F*)(file1->Get("p_n_Jet_neEmEF"));
    TH1F* data_p_n_Jet_leadTrackPt = (TH1F*)(file1->Get("p_n_Jet_leadTrackPt"));
    TH1F* data_p_n_Jet_leptonPt = (TH1F*)(file1->Get("p_n_Jet_leptonPt"));
    TH1F* data_p_n_Jet_leptonPtRel = (TH1F*)(file1->Get("p_n_Jet_leptonPtRel"));
    TH1F* data_p_n_Jet_leptonPdgId = (TH1F*)(file1->Get("p_n_Jet_leptonPdgId"));
    TH1F* data_p_n_Jet_leptonPtRelInv = (TH1F*)(file1->Get("p_n_Jet_leptonPtRelInv"));
    TH1F* data_p_n_Jet_vtxMass = (TH1F*)(file1->Get("p_n_Jet_vtxMass"));
    TH1F* data_p_n_Jet_vtxNtracks = (TH1F*)(file1->Get("p_n_Jet_vtxNtracks"));
    TH1F* data_p_n_Jet_vtxPt = (TH1F*)(file1->Get("p_n_Jet_vtxPt"));
    TH1F* data_p_n_Jet_vtx3DSig = (TH1F*)(file1->Get("p_n_Jet_vtx3DSig"));
    TH1F* data_p_n_Jet_vtx3DVal = (TH1F*)(file1->Get("p_n_Jet_vtx3DVal"));
    TH1F* data_p_n_Jet_ptd = (TH1F*)(file1->Get("p_n_Jet_ptd"));
    TH1F* data_p_n_rho = (TH1F*)(file1->Get("p_n_rho"));

    TH1F* r_data_p_Jet_leptonDeltaR = (TH1F*)(r_file1->Get("r_p_Jet_leptonDeltaR"));
    TH1F* r_data_p_Jet_mass = (TH1F*)(r_file1->Get("r_p_Jet_mass"));
    TH1F* r_data_p_Jet_pt = (TH1F*)(r_file1->Get("r_p_Jet_pt"));
    TH1F* r_data_p_Jet_eta = (TH1F*)(r_file1->Get("r_p_Jet_eta"));
    TH1F* r_data_p_Jet_phi = (TH1F*)(r_file1->Get("r_p_Jet_phi"));
    TH1F* r_data_p_Jet_energyRing_dR0_em = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR0_em"));
    TH1F* r_data_p_Jet_energyRing_dR1_em = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR1_em"));
    TH1F* r_data_p_Jet_energyRing_dR2_em = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR2_em"));
    TH1F* r_data_p_Jet_energyRing_dR3_em = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR3_em"));
    TH1F* r_data_p_Jet_energyRing_dR4_em = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR4_em"));
    TH1F* r_data_p_Jet_energyRing_dR0_mu = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR0_mu"));
    TH1F* r_data_p_Jet_energyRing_dR1_mu = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR1_mu"));
    TH1F* r_data_p_Jet_energyRing_dR2_mu = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR2_mu"));
    TH1F* r_data_p_Jet_energyRing_dR3_mu = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR3_mu"));
    TH1F* r_data_p_Jet_energyRing_dR4_mu = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR4_mu"));
    TH1F* r_data_p_Jet_energyRing_dR0_ch = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR0_ch"));
    TH1F* r_data_p_Jet_energyRing_dR1_ch = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR1_ch"));
    TH1F* r_data_p_Jet_energyRing_dR2_ch = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR2_ch"));
    TH1F* r_data_p_Jet_energyRing_dR3_ch = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR3_ch"));
    TH1F* r_data_p_Jet_energyRing_dR4_ch = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR4_ch"));
    TH1F* r_data_p_Jet_energyRing_dR0_neut = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR0_neut"));
    TH1F* r_data_p_Jet_energyRing_dR1_neut = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR1_neut"));
    TH1F* r_data_p_Jet_energyRing_dR2_neut = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR2_neut"));
    TH1F* r_data_p_Jet_energyRing_dR3_neut = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR3_neut"));
    TH1F* r_data_p_Jet_energyRing_dR4_neut = (TH1F*)(r_file1->Get("r_p_Jet_energyRing_dR4_neut"));
    TH1F* r_data_p_Jet_rawEnergy = (TH1F*)(r_file1->Get("r_p_Jet_rawEnergy"));
    TH1F* r_data_p_Jet_numDaughters_pt03 = (TH1F*)(r_file1->Get("r_p_Jet_numDaughters_pt03"));
    TH1F* r_data_p_Jet_numberOfDaughters = (TH1F*)(r_file1->Get("r_p_Jet_numberOfDaughters"));
    TH1F* r_data_p_Jet_rawPt = (TH1F*)(r_file1->Get("r_p_Jet_rawPt"));
    TH1F* r_data_p_Jet_chHEF = (TH1F*)(r_file1->Get("r_p_Jet_chHEF"));
    TH1F* r_data_p_Jet_muEF = (TH1F*)(r_file1->Get("r_p_Jet_muEF"));
    TH1F* r_data_p_Jet_neHEF = (TH1F*)(r_file1->Get("r_p_Jet_neHEF"));
    TH1F* r_data_p_Jet_chEmEF = (TH1F*)(r_file1->Get("r_p_Jet_chEmEF"));
    TH1F* r_data_p_Jet_neEmEF = (TH1F*)(r_file1->Get("r_p_Jet_neEmEF"));
    TH1F* r_data_p_Jet_leadTrackPt = (TH1F*)(r_file1->Get("r_p_Jet_leadTrackPt"));
    TH1F* r_data_p_Jet_leptonPt = (TH1F*)(r_file1->Get("r_p_Jet_leptonPt"));
    TH1F* r_data_p_Jet_leptonPtRel = (TH1F*)(r_file1->Get("r_p_Jet_leptonPtRel"));
    TH1F* r_data_p_Jet_leptonPdgId = (TH1F*)(r_file1->Get("r_p_Jet_leptonPdgId"));
    TH1F* r_data_p_Jet_leptonPtRelInv = (TH1F*)(r_file1->Get("r_p_Jet_leptonPtRelInv"));
    TH1F* r_data_p_Jet_vtxMass = (TH1F*)(r_file1->Get("r_p_Jet_vtxMass"));
    TH1F* r_data_p_Jet_vtxNtracks = (TH1F*)(r_file1->Get("r_p_Jet_vtxNtracks"));
    TH1F* r_data_p_Jet_vtxPt = (TH1F*)(r_file1->Get("r_p_Jet_vtxPt"));
    TH1F* r_data_p_Jet_vtx3DSig = (TH1F*)(r_file1->Get("r_p_Jet_vtx3DSig"));
    TH1F* r_data_p_Jet_vtx3DVal = (TH1F*)(r_file1->Get("r_p_Jet_vtx3DVal"));
    TH1F* r_data_p_Jet_ptd = (TH1F*)(r_file1->Get("r_p_Jet_ptd"));
    TH1F* r_data_p_rho = (TH1F*)(r_file1->Get("r_p_rho"));

    RooRealVar BDT_response("BDT_response", "BDT_response",-1,1);

    RooRealVar m_nl_j("m_nl_j", "m_nl_j",0,3000);
    RooRealVar jet_pt("jet_pt", "jet_pt",0,3000);
    RooRealVar Jet_mcFlavour("Jet_mcFlavour","Jet_mcFlavour",-30,30);
    RooRealVar jet_indices("jet_indices", "jet_indices",-100,100);

    RooRealVar Jet_leptonDeltaR("Jet_leptonDeltaR","Jet_leptonDeltaR",-105,5);
    RooRealVar Jet_mass("Jet_mass","Jet_mass",0,230);
    RooRealVar Jet_pt("Jet_pt","Jet_pt",0,2000);
    RooRealVar Jet_eta("Jet_eta","Jet_eta",-3,3);
    RooRealVar Jet_phi("Jet_phi","Jet_phi",-4,4);
    RooRealVar Jet_energyRing_dR0_neut("Jet_energyRing_dR0_neut","Jet_energyRing_dR0_neut",0,1);
    RooRealVar Jet_energyRing_dR1_neut("Jet_energyRing_dR1_neut","Jet_energyRing_dR1_neut",0,1);
    RooRealVar Jet_energyRing_dR2_neut("Jet_energyRing_dR2_neut","Jet_energyRing_dR2_neut",0,1);
    RooRealVar Jet_energyRing_dR3_neut("Jet_energyRing_dR3_neut","Jet_energyRing_dR3_neut",0,1);
    RooRealVar Jet_energyRing_dR4_neut("Jet_energyRing_dR4_neut","Jet_energyRing_dR4_neut",0,1);
    RooRealVar Jet_energyRing_dR0_ch("Jet_energyRing_dR0_ch","Jet_energyRing_dR0_ch",0,1);
    RooRealVar Jet_energyRing_dR1_ch("Jet_energyRing_dR1_ch","Jet_energyRing_dR1_ch",0,1);
    RooRealVar Jet_energyRing_dR2_ch("Jet_energyRing_dR2_ch","Jet_energyRing_dR2_ch",0,1);
    RooRealVar Jet_energyRing_dR3_ch("Jet_energyRing_dR3_ch","Jet_energyRing_dR3_ch",0,1);
    RooRealVar Jet_energyRing_dR4_ch("Jet_energyRing_dR4_ch","Jet_energyRing_dR4_ch",0,1);
    RooRealVar Jet_energyRing_dR0_em("Jet_energyRing_dR0_em","Jet_energyRing_dR0_em",0,1);
    RooRealVar Jet_energyRing_dR1_em("Jet_energyRing_dR1_em","Jet_energyRing_dR1_em",0,1);
    RooRealVar Jet_energyRing_dR2_em("Jet_energyRing_dR2_em","Jet_energyRing_dR2_em",0,1);
    RooRealVar Jet_energyRing_dR3_em("Jet_energyRing_dR3_em","Jet_energyRing_dR3_em",0,1);
    RooRealVar Jet_energyRing_dR4_em("Jet_energyRing_dR4_em","Jet_energyRing_dR4_em",0,1);
    RooRealVar Jet_energyRing_dR0_mu("Jet_energyRing_dR0_mu","Jet_energyRing_dR0_mu",0,1);
    RooRealVar Jet_energyRing_dR1_mu("Jet_energyRing_dR1_mu","Jet_energyRing_dR1_mu",0,1);
    RooRealVar Jet_energyRing_dR2_mu("Jet_energyRing_dR2_mu","Jet_energyRing_dR2_mu",0,1);
    RooRealVar Jet_energyRing_dR3_mu("Jet_energyRing_dR3_mu","Jet_energyRing_dR3_mu",0,1);
    RooRealVar Jet_energyRing_dR4_mu("Jet_energyRing_dR4_mu","Jet_energyRing_dR4_mu",0,1);
    RooRealVar Jet_rawEnergy("Jet_rawEnergy","Jet_rawEnergy",0,4000);
    RooRealVar Jet_numDaughters_pt03("Jet_numDaughters_pt03","Jet_numDaughters_pt03",0,100);
    RooRealVar Jet_numberOfDaughters("Jet_numberOfDaughters","Jet_numberOfDaughters",0,1200000000);
    RooRealVar rho("rho","rho",-200,100);
    RooRealVar Jet_rawPt("Jet_rawPt","Jet_rawPt",0,1600);
    RooRealVar Jet_chHEF("Jet_chHEF","Jet_chHEF",0,12);
    RooRealVar Jet_muEF("Jet_muEF","Jet_muEF",0,12);
    RooRealVar Jet_neHEF("Jet_neHEF","Jet_neHEF",0,12);
    RooRealVar Jet_chEmEF("Jet_chEmEF","Jet_chEmEF",0,12);
    RooRealVar Jet_neEmEF("Jet_neEmEF","Jet_neEmEF",0,12);
    RooRealVar Jet_leadTrackPt("Jet_leadTrackPt","Jet_leadTrackPt",0,800);
    RooRealVar Jet_leptonPt("Jet_leptonPt","Jet_leptonPt",-200,1000);
    RooRealVar Jet_leptonPtRel("Jet_leptonPtRel","Jet_leptonPtRel",-200,100);
    RooRealVar Jet_leptonPdgId("Jet_leptonPdgId","Jet_leptonPdgId",-200,100);
    RooRealVar Jet_leptonPtRelInv("Jet_leptonPtRelInv","Jet_leptonPtRelInv",-200,200);
    RooRealVar Jet_vtxMass("Jet_vtxMass","Jet_vtxMass",0,10);
    RooRealVar Jet_vtxNtracks("Jet_vtxNtracks","Jet_vtxNtracks",0,16);
    RooRealVar Jet_vtxPt("Jet_vtxPt","Jet_vtxPt",0,600);
    RooRealVar Jet_vtx3DSig("Jet_vtx3DSig","Jet_vtx3DSig",0,400);
    RooRealVar Jet_vtx3DVal("Jet_vtx3DVal","Jet_vtx3DVal",0,18);
    RooRealVar Jet_ptd("Jet_ptd","Jet_ptd",-30,20);

    RooArgSet var_set;
    var_set.add(Jet_mcFlavour);
    var_set.add(BDT_response);
    var_set.add(Jet_mass);
    var_set.add(Jet_pt);
    var_set.add(Jet_eta);
    var_set.add(Jet_phi);
    var_set.add(Jet_rawEnergy);
    var_set.add(Jet_energyRing_dR0_neut);   
    var_set.add(Jet_energyRing_dR1_neut);
    var_set.add(Jet_energyRing_dR2_neut);
    var_set.add(Jet_energyRing_dR3_neut);
    var_set.add(Jet_energyRing_dR4_neut);
    var_set.add(Jet_energyRing_dR0_ch);
    var_set.add(Jet_energyRing_dR1_ch);
    var_set.add(Jet_energyRing_dR2_ch);
    var_set.add(Jet_energyRing_dR3_ch);
    var_set.add(Jet_energyRing_dR4_ch);
    var_set.add(Jet_energyRing_dR0_em);
    var_set.add(Jet_energyRing_dR1_em);
    var_set.add(Jet_energyRing_dR2_em);
    var_set.add(Jet_energyRing_dR3_em);
    var_set.add(Jet_energyRing_dR4_em);
    var_set.add(Jet_energyRing_dR0_mu);
    var_set.add(Jet_energyRing_dR1_mu);
    var_set.add(Jet_energyRing_dR2_mu);
    var_set.add(Jet_energyRing_dR3_mu);
    var_set.add(Jet_energyRing_dR4_mu);
    var_set.add(Jet_numDaughters_pt03);
    var_set.add(Jet_numberOfDaughters);
    var_set.add(rho);
    var_set.add(Jet_rawPt);
    var_set.add(Jet_chHEF);
    var_set.add(Jet_muEF);
    var_set.add(Jet_neHEF);
    var_set.add(Jet_chEmEF);
    var_set.add(Jet_neEmEF);
    var_set.add(Jet_leadTrackPt);
    var_set.add(Jet_leptonPt);
    var_set.add(Jet_leptonPtRel);
    var_set.add(Jet_leptonPdgId);
    var_set.add(Jet_leptonPtRelInv);
    var_set.add(Jet_vtxMass);
    var_set.add(Jet_vtxNtracks);
    var_set.add(Jet_vtxPt);
    var_set.add(Jet_vtx3DSig);
    var_set.add(Jet_vtx3DVal);
    var_set.add(Jet_ptd);

    RooArgSet r_var_set;
    r_var_set.add(BDT_response);
    r_var_set.add(Jet_mass);
    r_var_set.add(Jet_pt);
    r_var_set.add(Jet_eta);
    r_var_set.add(Jet_phi);
    r_var_set.add(Jet_rawEnergy);
    r_var_set.add(Jet_energyRing_dR0_neut);
    r_var_set.add(Jet_energyRing_dR1_neut);
    r_var_set.add(Jet_energyRing_dR2_neut);
    r_var_set.add(Jet_energyRing_dR3_neut);
    r_var_set.add(Jet_energyRing_dR4_neut);
    r_var_set.add(Jet_energyRing_dR0_ch);
    r_var_set.add(Jet_energyRing_dR1_ch);
    r_var_set.add(Jet_energyRing_dR2_ch);
    r_var_set.add(Jet_energyRing_dR3_ch);
    r_var_set.add(Jet_energyRing_dR4_ch);
    r_var_set.add(Jet_energyRing_dR0_em);
    r_var_set.add(Jet_energyRing_dR1_em);
    r_var_set.add(Jet_energyRing_dR2_em);
    r_var_set.add(Jet_energyRing_dR3_em);
    r_var_set.add(Jet_energyRing_dR4_em);
    r_var_set.add(Jet_energyRing_dR0_mu);
    r_var_set.add(Jet_energyRing_dR1_mu);
    r_var_set.add(Jet_energyRing_dR2_mu);
    r_var_set.add(Jet_energyRing_dR3_mu);
    r_var_set.add(Jet_energyRing_dR4_mu);
    r_var_set.add(Jet_numDaughters_pt03);
    r_var_set.add(Jet_numberOfDaughters);
    r_var_set.add(rho);
    r_var_set.add(Jet_rawPt);
    r_var_set.add(Jet_chHEF);
    r_var_set.add(Jet_muEF);
    r_var_set.add(Jet_neHEF);
    r_var_set.add(Jet_chEmEF);
    r_var_set.add(Jet_neEmEF);
    r_var_set.add(Jet_leadTrackPt);
    r_var_set.add(Jet_leptonPt);
    r_var_set.add(Jet_leptonPtRel);
    r_var_set.add(Jet_leptonPdgId);
    r_var_set.add(Jet_leptonPtRelInv);
    r_var_set.add(Jet_vtxMass);
    r_var_set.add(Jet_vtxNtracks);
    r_var_set.add(Jet_vtxPt);
    r_var_set.add(Jet_vtx3DSig);
    r_var_set.add(Jet_vtx3DVal);
    r_var_set.add(Jet_ptd);

    RooDataSet x_data_set("data_set","data_set",var_set,Import(*datatree));//,"Jet_pt<120 && Jet_pt>60");
    RooDataSet x_r_data_set("r_data_set","r_data_set",r_var_set,Import(*r_datatree));//,"Jet_pt<120 && Jet_pt>60");

    RooDataSet* ptr_data_set = (RooDataSet*) x_data_set.reduce(var_set);//,"Jet_pt>=120");// && Jet_pt<120");
    RooDataSet data_set = *ptr_data_set;
   
    RooDataSet* ptr_r_data_set = (RooDataSet*) x_r_data_set.reduce(var_set);//,"Jet_pt>=120");// && Jet_pt<120");
    RooDataSet r_data_set = *ptr_r_data_set;

    cout<<endl<<" is x_dataset events = "<<data_set.sumEntries()<<endl;
    cout<<endl<<" is x_r_dataset events = "<<r_data_set.sumEntries()<<endl;

    RooDataHist p_b_bdt("p_b_bdt","p_b_bdt",RooArgList(BDT_response),Import(*data_p_b_bdt));
    RooDataHist p_n_bdt("p_n_bdt","p_n_bdt",RooArgList(BDT_response),Import(*data_p_n_bdt));
    RooDataHist p_bdt("p_bdt","p_bdt",RooArgList(BDT_response),Import(*data_p_bdt));
    RooDataHist r_p_bdt("r_p_bdt","r_p_bdt",RooArgList(BDT_response),Import(*r_data_p_bdt));

    RooDataHist r_p_30_bdt("r_p_30_bdt","r_p_30_bdt",RooArgList(BDT_response),Import(*r_data_p_30_bdt));
    RooDataHist p_b_30_bdt("p_b_30_bdt","p_b_30_bdt",RooArgList(BDT_response),Import(*data_p_b_30_bdt));
    RooDataHist p_n_30_bdt("p_n_30_bdt","p_n_30_bdt",RooArgList(BDT_response),Import(*data_p_n_30_bdt));
    RooDataHist r_p_40_bdt("r_p_40_bdt","r_p_40_bdt",RooArgList(BDT_response),Import(*r_data_p_40_bdt));
    RooDataHist p_b_40_bdt("p_b_40_bdt","p_b_40_bdt",RooArgList(BDT_response),Import(*data_p_b_40_bdt));
    RooDataHist p_n_40_bdt("p_n_40_bdt","p_n_40_bdt",RooArgList(BDT_response),Import(*data_p_n_40_bdt));
    RooDataHist r_p_60_bdt("r_p_60_bdt","r_p_60_bdt",RooArgList(BDT_response),Import(*r_data_p_60_bdt));
    RooDataHist p_b_60_bdt("p_b_60_bdt","p_b_60_bdt",RooArgList(BDT_response),Import(*data_p_b_60_bdt));
    RooDataHist p_n_60_bdt("p_n_60_bdt","p_n_60_bdt",RooArgList(BDT_response),Import(*data_p_n_60_bdt));
    RooDataHist r_p_80_bdt("r_p_80_bdt","r_p_80_bdt",RooArgList(BDT_response),Import(*r_data_p_80_bdt));
    RooDataHist p_b_80_bdt("p_b_80_bdt","p_b_80_bdt",RooArgList(BDT_response),Import(*data_p_b_80_bdt));
    RooDataHist p_n_80_bdt("p_n_80_bdt","p_n_80_bdt",RooArgList(BDT_response),Import(*data_p_n_80_bdt));
    RooDataHist r_p_100_bdt("r_p_100_bdt","r_p_100_bdt",RooArgList(BDT_response),Import(*r_data_p_100_bdt));
    RooDataHist p_b_100_bdt("p_b_100_bdt","p_b_100_bdt",RooArgList(BDT_response),Import(*data_p_b_100_bdt));
    RooDataHist p_n_100_bdt("p_n_100_bdt","p_n_100_bdt",RooArgList(BDT_response),Import(*data_p_n_100_bdt));
    RooDataHist r_p_120_bdt("r_p_120_bdt","r_p_120_bdt",RooArgList(BDT_response),Import(*r_data_p_120_bdt));
    RooDataHist p_b_120_bdt("p_b_120_bdt","p_b_120_bdt",RooArgList(BDT_response),Import(*data_p_b_120_bdt));
    RooDataHist p_n_120_bdt("p_n_120_bdt","p_n_120_bdt",RooArgList(BDT_response),Import(*data_p_n_120_bdt));

    RooDataHist p_b_60_120_bdt("p_b_60_120_bdt","p_b_60_120_bdt",RooArgList(BDT_response),Import(*data_p_b_60_120_bdt));
    RooDataHist p_n_60_120_bdt("p_n_60_120_bdt","p_n_60_120_bdt",RooArgList(BDT_response),Import(*data_p_n_60_120_bdt));
    RooDataHist r_p_60_120_bdt("r_p_60_120_bdt","r_p_60_120_bdt",RooArgList(BDT_response),Import(*r_data_p_60_120_bdt));

    RooDataHist p_b_Jet_leptonDeltaR("p_b_Jet_leptonDeltaR","p_b_Jet_leptonDeltaR",RooArgList(Jet_leptonDeltaR),Import(*data_p_b_Jet_leptonDeltaR));
    RooDataHist p_n_Jet_leptonDeltaR("p_n_Jet_leptonDeltaR","p_n_Jet_leptonDeltaR",RooArgList(Jet_leptonDeltaR),Import(*data_p_n_Jet_leptonDeltaR));
    RooDataHist r_p_Jet_leptonDeltaR("r_p_Jet_leptonDeltaR","r_p_Jet_leptonDeltaR",RooArgList(Jet_leptonDeltaR),Import(*r_data_p_Jet_leptonDeltaR));

    RooDataHist p_b_Jet_pt("p_b_Jet_pt","p_b_Jet_pt",RooArgList(Jet_pt),Import(*data_p_b_Jet_pt));
    RooDataHist p_n_Jet_pt("p_n_Jet_pt","p_n_Jet_pt",RooArgList(Jet_pt),Import(*data_p_n_Jet_pt));
    RooDataHist r_p_Jet_pt("r_p_Jet_pt","r_p_Jet_pt",RooArgList(Jet_pt),Import(*r_data_p_Jet_pt));

    RooDataHist p_b_Jet_phi("p_b_Jet_phi","p_b_Jet_phi",RooArgList(Jet_phi),Import(*data_p_b_Jet_phi));
    RooDataHist p_n_Jet_phi("p_n_Jet_phi","p_n_Jet_phi",RooArgList(Jet_phi),Import(*data_p_n_Jet_phi));
    RooDataHist r_p_Jet_phi("r_p_Jet_phi","r_p_Jet_phi",RooArgList(Jet_phi),Import(*r_data_p_Jet_phi));

    RooDataHist p_b_Jet_eta("p_b_Jet_eta","p_b_Jet_eta",RooArgList(Jet_eta),Import(*data_p_b_Jet_eta));
    RooDataHist p_n_Jet_eta("p_n_Jet_eta","p_n_Jet_eta",RooArgList(Jet_eta),Import(*data_p_n_Jet_eta));
    RooDataHist r_p_Jet_eta("r_p_Jet_eta","r_p_Jet_eta",RooArgList(Jet_eta),Import(*r_data_p_Jet_eta));

    RooDataHist p_b_Jet_mass("p_b_Jet_mass","p_b_Jet_mass",RooArgList(Jet_mass),Import(*data_p_b_Jet_mass));
    RooDataHist p_n_Jet_mass("p_n_Jet_mass","p_n_Jet_mass",RooArgList(Jet_mass),Import(*data_p_n_Jet_mass));
    RooDataHist r_p_Jet_mass("r_p_Jet_mass","r_p_Jet_mass",RooArgList(Jet_mass),Import(*r_data_p_Jet_mass));

    RooDataHist p_b_Jet_rawEnergy("p_b_Jet_rawEnergy","p_b_Jet_rawEnergy",RooArgList(Jet_rawEnergy),Import(*data_p_b_Jet_rawEnergy));
    RooDataHist p_n_Jet_rawEnergy("p_n_Jet_rawEnergy","p_n_Jet_rawEnergy",RooArgList(Jet_rawEnergy),Import(*data_p_n_Jet_rawEnergy));
    RooDataHist r_p_Jet_rawEnergy("r_p_Jet_rawEnergy","r_p_Jet_rawEnergy",RooArgList(Jet_rawEnergy),Import(*r_data_p_Jet_rawEnergy));

    RooDataHist p_b_Jet_energyRing_dR0_em("p_b_Jet_energyRing_dR0_em","p_b_Jet_energyRing_dR0_em",RooArgList(Jet_energyRing_dR0_em),Import(*data_p_b_Jet_energyRing_dR0_em));
    RooDataHist p_n_Jet_energyRing_dR0_em("p_n_Jet_energyRing_dR0_em","p_n_Jet_energyRing_dR0_em",RooArgList(Jet_energyRing_dR0_em),Import(*data_p_n_Jet_energyRing_dR0_em));
    RooDataHist r_p_Jet_energyRing_dR0_em("r_p_Jet_energyRing_dR0_em","r_p_Jet_energyRing_dR0_em",RooArgList(Jet_energyRing_dR0_em),Import(*r_data_p_Jet_energyRing_dR0_em));
    RooDataHist p_b_Jet_energyRing_dR1_em("p_b_Jet_energyRing_dR1_em","p_b_Jet_energyRing_dR1_em",RooArgList(Jet_energyRing_dR1_em),Import(*data_p_b_Jet_energyRing_dR1_em));
    RooDataHist p_n_Jet_energyRing_dR1_em("p_n_Jet_energyRing_dR1_em","p_n_Jet_energyRing_dR1_em",RooArgList(Jet_energyRing_dR1_em),Import(*data_p_n_Jet_energyRing_dR1_em));
    RooDataHist r_p_Jet_energyRing_dR1_em("r_p_Jet_energyRing_dR1_em","r_p_Jet_energyRing_dR1_em",RooArgList(Jet_energyRing_dR1_em),Import(*r_data_p_Jet_energyRing_dR1_em));
    RooDataHist p_b_Jet_energyRing_dR2_em("p_b_Jet_energyRing_dR2_em","p_b_Jet_energyRing_dR2_em",RooArgList(Jet_energyRing_dR2_em),Import(*data_p_b_Jet_energyRing_dR2_em));
    RooDataHist p_n_Jet_energyRing_dR2_em("p_n_Jet_energyRing_dR2_em","p_n_Jet_energyRing_dR2_em",RooArgList(Jet_energyRing_dR2_em),Import(*data_p_n_Jet_energyRing_dR2_em));
    RooDataHist r_p_Jet_energyRing_dR2_em("r_p_Jet_energyRing_dR2_em","r_p_Jet_energyRing_dR2_em",RooArgList(Jet_energyRing_dR2_em),Import(*r_data_p_Jet_energyRing_dR2_em));
    RooDataHist p_b_Jet_energyRing_dR3_em("p_b_Jet_energyRing_dR3_em","p_b_Jet_energyRing_dR3_em",RooArgList(Jet_energyRing_dR3_em),Import(*data_p_b_Jet_energyRing_dR3_em));
    RooDataHist p_n_Jet_energyRing_dR3_em("p_n_Jet_energyRing_dR3_em","p_n_Jet_energyRing_dR3_em",RooArgList(Jet_energyRing_dR3_em),Import(*data_p_n_Jet_energyRing_dR3_em));
    RooDataHist r_p_Jet_energyRing_dR3_em("r_p_Jet_energyRing_dR3_em","r_p_Jet_energyRing_dR3_em",RooArgList(Jet_energyRing_dR3_em),Import(*r_data_p_Jet_energyRing_dR3_em));
    RooDataHist p_b_Jet_energyRing_dR4_em("p_b_Jet_energyRing_dR4_em","p_b_Jet_energyRing_dR4_em",RooArgList(Jet_energyRing_dR4_em),Import(*data_p_b_Jet_energyRing_dR4_em));
    RooDataHist p_n_Jet_energyRing_dR4_em("p_n_Jet_energyRing_dR4_em","p_n_Jet_energyRing_dR4_em",RooArgList(Jet_energyRing_dR4_em),Import(*data_p_n_Jet_energyRing_dR4_em));
    RooDataHist r_p_Jet_energyRing_dR4_em("r_p_Jet_energyRing_dR4_em","r_p_Jet_energyRing_dR4_em",RooArgList(Jet_energyRing_dR4_em),Import(*r_data_p_Jet_energyRing_dR4_em));


    RooDataHist p_b_Jet_energyRing_dR0_mu("p_b_Jet_energyRing_dR0_mu","p_b_Jet_energyRing_dR0_mu",RooArgList(Jet_energyRing_dR0_mu),Import(*data_p_b_Jet_energyRing_dR0_mu));
    RooDataHist p_n_Jet_energyRing_dR0_mu("p_n_Jet_energyRing_dR0_mu","p_n_Jet_energyRing_dR0_mu",RooArgList(Jet_energyRing_dR0_mu),Import(*data_p_n_Jet_energyRing_dR0_mu));
    RooDataHist r_p_Jet_energyRing_dR0_mu("r_p_Jet_energyRing_dR0_mu","r_p_Jet_energyRing_dR0_mu",RooArgList(Jet_energyRing_dR0_mu),Import(*r_data_p_Jet_energyRing_dR0_mu));
    RooDataHist p_b_Jet_energyRing_dR1_mu("p_b_Jet_energyRing_dR1_mu","p_b_Jet_energyRing_dR1_mu",RooArgList(Jet_energyRing_dR1_mu),Import(*data_p_b_Jet_energyRing_dR1_mu));
    RooDataHist p_n_Jet_energyRing_dR1_mu("p_n_Jet_energyRing_dR1_mu","p_n_Jet_energyRing_dR1_mu",RooArgList(Jet_energyRing_dR1_mu),Import(*data_p_n_Jet_energyRing_dR1_mu));
    RooDataHist r_p_Jet_energyRing_dR1_mu("r_p_Jet_energyRing_dR1_mu","r_p_Jet_energyRing_dR1_mu",RooArgList(Jet_energyRing_dR1_mu),Import(*r_data_p_Jet_energyRing_dR1_mu));
    RooDataHist p_b_Jet_energyRing_dR2_mu("p_b_Jet_energyRing_dR2_mu","p_b_Jet_energyRing_dR2_mu",RooArgList(Jet_energyRing_dR2_mu),Import(*data_p_b_Jet_energyRing_dR2_mu));
    RooDataHist p_n_Jet_energyRing_dR2_mu("p_n_Jet_energyRing_dR2_mu","p_n_Jet_energyRing_dR2_mu",RooArgList(Jet_energyRing_dR2_mu),Import(*data_p_n_Jet_energyRing_dR2_mu));
    RooDataHist r_p_Jet_energyRing_dR2_mu("r_p_Jet_energyRing_dR2_mu","r_p_Jet_energyRing_dR2_mu",RooArgList(Jet_energyRing_dR2_mu),Import(*r_data_p_Jet_energyRing_dR2_mu));
    RooDataHist p_b_Jet_energyRing_dR3_mu("p_b_Jet_energyRing_dR3_mu","p_b_Jet_energyRing_dR3_mu",RooArgList(Jet_energyRing_dR3_mu),Import(*data_p_b_Jet_energyRing_dR3_mu));
    RooDataHist p_n_Jet_energyRing_dR3_mu("p_n_Jet_energyRing_dR3_mu","p_n_Jet_energyRing_dR3_mu",RooArgList(Jet_energyRing_dR3_mu),Import(*data_p_n_Jet_energyRing_dR3_mu));
    RooDataHist r_p_Jet_energyRing_dR3_mu("r_p_Jet_energyRing_dR3_mu","r_p_Jet_energyRing_dR3_mu",RooArgList(Jet_energyRing_dR3_mu),Import(*r_data_p_Jet_energyRing_dR3_mu));
    RooDataHist p_b_Jet_energyRing_dR4_mu("p_b_Jet_energyRing_dR4_mu","p_b_Jet_energyRing_dR4_mu",RooArgList(Jet_energyRing_dR4_mu),Import(*data_p_b_Jet_energyRing_dR4_mu));
    RooDataHist p_n_Jet_energyRing_dR4_mu("p_n_Jet_energyRing_dR4_mu","p_n_Jet_energyRing_dR4_mu",RooArgList(Jet_energyRing_dR4_mu),Import(*data_p_n_Jet_energyRing_dR4_mu));
    RooDataHist r_p_Jet_energyRing_dR4_mu("r_p_Jet_energyRing_dR4_mu","r_p_Jet_energyRing_dR4_mu",RooArgList(Jet_energyRing_dR4_mu),Import(*r_data_p_Jet_energyRing_dR4_mu));


    RooDataHist p_b_Jet_energyRing_dR0_ch("p_b_Jet_energyRing_dR0_ch","p_b_Jet_energyRing_dR0_ch",RooArgList(Jet_energyRing_dR0_ch),Import(*data_p_b_Jet_energyRing_dR0_ch));
    RooDataHist p_n_Jet_energyRing_dR0_ch("p_n_Jet_energyRing_dR0_ch","p_n_Jet_energyRing_dR0_ch",RooArgList(Jet_energyRing_dR0_ch),Import(*data_p_n_Jet_energyRing_dR0_ch));
    RooDataHist r_p_Jet_energyRing_dR0_ch("r_p_Jet_energyRing_dR0_ch","r_p_Jet_energyRing_dR0_ch",RooArgList(Jet_energyRing_dR0_ch),Import(*r_data_p_Jet_energyRing_dR0_ch));
    RooDataHist p_b_Jet_energyRing_dR1_ch("p_b_Jet_energyRing_dR1_ch","p_b_Jet_energyRing_dR1_ch",RooArgList(Jet_energyRing_dR1_ch),Import(*data_p_b_Jet_energyRing_dR1_ch));
    RooDataHist p_n_Jet_energyRing_dR1_ch("p_n_Jet_energyRing_dR1_ch","p_n_Jet_energyRing_dR1_ch",RooArgList(Jet_energyRing_dR1_ch),Import(*data_p_n_Jet_energyRing_dR1_ch));
    RooDataHist r_p_Jet_energyRing_dR1_ch("r_p_Jet_energyRing_dR1_ch","r_p_Jet_energyRing_dR1_ch",RooArgList(Jet_energyRing_dR1_ch),Import(*r_data_p_Jet_energyRing_dR1_ch));
    RooDataHist p_b_Jet_energyRing_dR2_ch("p_b_Jet_energyRing_dR2_ch","p_b_Jet_energyRing_dR2_ch",RooArgList(Jet_energyRing_dR2_ch),Import(*data_p_b_Jet_energyRing_dR2_ch));
    RooDataHist p_n_Jet_energyRing_dR2_ch("p_n_Jet_energyRing_dR2_ch","p_n_Jet_energyRing_dR2_ch",RooArgList(Jet_energyRing_dR2_ch),Import(*data_p_n_Jet_energyRing_dR2_ch));
    RooDataHist r_p_Jet_energyRing_dR2_ch("r_p_Jet_energyRing_dR2_ch","r_p_Jet_energyRing_dR2_ch",RooArgList(Jet_energyRing_dR2_ch),Import(*r_data_p_Jet_energyRing_dR2_ch));
    RooDataHist p_b_Jet_energyRing_dR3_ch("p_b_Jet_energyRing_dR3_ch","p_b_Jet_energyRing_dR3_ch",RooArgList(Jet_energyRing_dR3_ch),Import(*data_p_b_Jet_energyRing_dR3_ch));
    RooDataHist p_n_Jet_energyRing_dR3_ch("p_n_Jet_energyRing_dR3_ch","p_n_Jet_energyRing_dR3_ch",RooArgList(Jet_energyRing_dR3_ch),Import(*data_p_n_Jet_energyRing_dR3_ch));
    RooDataHist r_p_Jet_energyRing_dR3_ch("r_p_Jet_energyRing_dR3_ch","r_p_Jet_energyRing_dR3_ch",RooArgList(Jet_energyRing_dR3_ch),Import(*r_data_p_Jet_energyRing_dR3_ch));
    RooDataHist p_b_Jet_energyRing_dR4_ch("p_b_Jet_energyRing_dR4_ch","p_b_Jet_energyRing_dR4_ch",RooArgList(Jet_energyRing_dR4_ch),Import(*data_p_b_Jet_energyRing_dR4_ch));
    RooDataHist p_n_Jet_energyRing_dR4_ch("p_n_Jet_energyRing_dR4_ch","p_n_Jet_energyRing_dR4_ch",RooArgList(Jet_energyRing_dR4_ch),Import(*data_p_n_Jet_energyRing_dR4_ch));
    RooDataHist r_p_Jet_energyRing_dR4_ch("r_p_Jet_energyRing_dR4_ch","r_p_Jet_energyRing_dR4_ch",RooArgList(Jet_energyRing_dR4_ch),Import(*r_data_p_Jet_energyRing_dR4_ch));

    RooDataHist p_b_Jet_energyRing_dR0_neut("p_b_Jet_energyRing_dR0_neut","p_b_Jet_energyRing_dR0_neut",RooArgList(Jet_energyRing_dR0_neut),Import(*data_p_b_Jet_energyRing_dR0_neut));
    RooDataHist p_n_Jet_energyRing_dR0_neut("p_n_Jet_energyRing_dR0_neut","p_n_Jet_energyRing_dR0_neut",RooArgList(Jet_energyRing_dR0_neut),Import(*data_p_n_Jet_energyRing_dR0_neut));
    RooDataHist r_p_Jet_energyRing_dR0_neut("r_p_Jet_energyRing_dR0_neut","r_p_Jet_energyRing_dR0_neut",RooArgList(Jet_energyRing_dR0_neut),Import(*r_data_p_Jet_energyRing_dR0_neut));
    RooDataHist p_b_Jet_energyRing_dR1_neut("p_b_Jet_energyRing_dR1_neut","p_b_Jet_energyRing_dR1_neut",RooArgList(Jet_energyRing_dR1_neut),Import(*data_p_b_Jet_energyRing_dR1_neut));
    RooDataHist p_n_Jet_energyRing_dR1_neut("p_n_Jet_energyRing_dR1_neut","p_n_Jet_energyRing_dR1_neut",RooArgList(Jet_energyRing_dR1_neut),Import(*data_p_n_Jet_energyRing_dR1_neut));
    RooDataHist r_p_Jet_energyRing_dR1_neut("r_p_Jet_energyRing_dR1_neut","r_p_Jet_energyRing_dR1_neut",RooArgList(Jet_energyRing_dR1_neut),Import(*r_data_p_Jet_energyRing_dR1_neut));
    RooDataHist p_b_Jet_energyRing_dR2_neut("p_b_Jet_energyRing_dR2_neut","p_b_Jet_energyRing_dR2_neut",RooArgList(Jet_energyRing_dR2_neut),Import(*data_p_b_Jet_energyRing_dR2_neut));
    RooDataHist p_n_Jet_energyRing_dR2_neut("p_n_Jet_energyRing_dR2_neut","p_n_Jet_energyRing_dR2_neut",RooArgList(Jet_energyRing_dR2_neut),Import(*data_p_n_Jet_energyRing_dR2_neut));
    RooDataHist r_p_Jet_energyRing_dR2_neut("r_p_Jet_energyRing_dR2_neut","r_p_Jet_energyRing_dR2_neut",RooArgList(Jet_energyRing_dR2_neut),Import(*r_data_p_Jet_energyRing_dR2_neut));
    RooDataHist p_b_Jet_energyRing_dR3_neut("p_b_Jet_energyRing_dR3_neut","p_b_Jet_energyRing_dR3_neut",RooArgList(Jet_energyRing_dR3_neut),Import(*data_p_b_Jet_energyRing_dR3_neut));
    RooDataHist p_n_Jet_energyRing_dR3_neut("p_n_Jet_energyRing_dR3_neut","p_n_Jet_energyRing_dR3_neut",RooArgList(Jet_energyRing_dR3_neut),Import(*data_p_n_Jet_energyRing_dR3_neut));
    RooDataHist r_p_Jet_energyRing_dR3_neut("r_p_Jet_energyRing_dR3_neut","r_p_Jet_energyRing_dR3_neut",RooArgList(Jet_energyRing_dR3_neut),Import(*r_data_p_Jet_energyRing_dR3_neut));
    RooDataHist p_b_Jet_energyRing_dR4_neut("p_b_Jet_energyRing_dR4_neut","p_b_Jet_energyRing_dR4_neut",RooArgList(Jet_energyRing_dR4_neut),Import(*data_p_b_Jet_energyRing_dR4_neut));
    RooDataHist p_n_Jet_energyRing_dR4_neut("p_n_Jet_energyRing_dR4_neut","p_n_Jet_energyRing_dR4_neut",RooArgList(Jet_energyRing_dR4_neut),Import(*data_p_n_Jet_energyRing_dR4_neut));
    RooDataHist r_p_Jet_energyRing_dR4_neut("r_p_Jet_energyRing_dR4_neut","r_p_Jet_energyRing_dR4_neut",RooArgList(Jet_energyRing_dR4_neut),Import(*r_data_p_Jet_energyRing_dR4_neut));

    RooDataHist p_b_Jet_numDaughters_pt03("p_b_Jet_numDaughters_pt03","p_b_Jet_numDaughters_pt03",RooArgList(Jet_numDaughters_pt03),Import(*data_p_b_Jet_numDaughters_pt03));
    RooDataHist p_n_Jet_numDaughters_pt03("p_n_Jet_numDaughters_pt03","p_n_Jet_numDaughters_pt03",RooArgList(Jet_numDaughters_pt03),Import(*data_p_n_Jet_numDaughters_pt03));
    RooDataHist r_p_Jet_numDaughters_pt03("r_p_Jet_numDaughters_pt03","r_p_Jet_numDaughters_pt03",RooArgList(Jet_numDaughters_pt03),Import(*r_data_p_Jet_numDaughters_pt03));

    RooDataHist p_b_Jet_numberOfDaughters("p_b_Jet_numberOfDaughters","p_b_Jet_numberOfDaughters",RooArgList(Jet_numberOfDaughters),Import(*data_p_b_Jet_numberOfDaughters));
    RooDataHist p_n_Jet_numberOfDaughters("p_n_Jet_numberOfDaughters","p_n_Jet_numberOfDaughters",RooArgList(Jet_numberOfDaughters),Import(*data_p_n_Jet_numberOfDaughters));
    RooDataHist r_p_Jet_numberOfDaughters("r_p_Jet_numberOfDaughters","r_p_Jet_numberOfDaughters",RooArgList(Jet_numberOfDaughters),Import(*r_data_p_Jet_numberOfDaughters));

    RooDataHist p_b_Jet_rawPt("p_b_Jet_rawPt","p_b_Jet_rawPt",RooArgList(Jet_rawPt),Import(*data_p_b_Jet_rawPt));
    RooDataHist p_n_Jet_rawPt("p_n_Jet_rawPt","p_n_Jet_rawPt",RooArgList(Jet_rawPt),Import(*data_p_n_Jet_rawPt));
    RooDataHist r_p_Jet_rawPt("r_p_Jet_rawPt","r_p_Jet_rawPt",RooArgList(Jet_rawPt),Import(*r_data_p_Jet_rawPt));

    RooDataHist p_b_Jet_chHEF("p_b_Jet_chHEF","p_b_Jet_chHEF",RooArgList(Jet_chHEF),Import(*data_p_b_Jet_chHEF));
    RooDataHist p_n_Jet_chHEF("p_n_Jet_chHEF","p_n_Jet_chHEF",RooArgList(Jet_chHEF),Import(*data_p_n_Jet_chHEF));
    RooDataHist r_p_Jet_chHEF("r_p_Jet_chHEF","r_p_Jet_chHEF",RooArgList(Jet_chHEF),Import(*r_data_p_Jet_chHEF));

    RooDataHist p_b_Jet_muEF("p_b_Jet_muEF","p_b_Jet_muEF",RooArgList(Jet_muEF),Import(*data_p_b_Jet_muEF));
    RooDataHist p_n_Jet_muEF("p_n_Jet_muEF","p_n_Jet_muEF",RooArgList(Jet_muEF),Import(*data_p_n_Jet_muEF));
    RooDataHist r_p_Jet_muEF("r_p_Jet_muEF","r_p_Jet_muEF",RooArgList(Jet_muEF),Import(*r_data_p_Jet_muEF));

    RooDataHist p_b_Jet_neHEF("p_b_Jet_neHEF","p_b_Jet_neHEF",RooArgList(Jet_neHEF),Import(*data_p_b_Jet_neHEF));
    RooDataHist p_n_Jet_neHEF("p_n_Jet_neHEF","p_n_Jet_neHEF",RooArgList(Jet_neHEF),Import(*data_p_n_Jet_neHEF));
    RooDataHist r_p_Jet_neHEF("r_p_Jet_neHEF","r_p_Jet_neHEF",RooArgList(Jet_neHEF),Import(*r_data_p_Jet_neHEF));

    RooDataHist p_b_Jet_chEmEF("p_b_Jet_chEmEF","p_b_Jet_chEmEF",RooArgList(Jet_chEmEF),Import(*data_p_b_Jet_chEmEF));
    RooDataHist p_n_Jet_chEmEF("p_n_Jet_chEmEF","p_n_Jet_chEmEF",RooArgList(Jet_chEmEF),Import(*data_p_n_Jet_chEmEF));
    RooDataHist r_p_Jet_chEmEF("r_p_Jet_chEmEF","r_p_Jet_chEmEF",RooArgList(Jet_chEmEF),Import(*r_data_p_Jet_chEmEF));

    RooDataHist p_b_Jet_neEmEF("p_b_Jet_neEmEF","p_b_Jet_neEmEF",RooArgList(Jet_neEmEF),Import(*data_p_b_Jet_neEmEF));
    RooDataHist p_n_Jet_neEmEF("p_n_Jet_neEmEF","p_n_Jet_neEmEF",RooArgList(Jet_neEmEF),Import(*data_p_n_Jet_neEmEF));
    RooDataHist r_p_Jet_neEmEF("r_p_Jet_neEmEF","r_p_Jet_neEmEF",RooArgList(Jet_neEmEF),Import(*r_data_p_Jet_neEmEF));

    RooDataHist p_b_Jet_leadTrackPt("p_b_Jet_leadTrackPt","p_b_Jet_leadTrackPt",RooArgList(Jet_leadTrackPt),Import(*data_p_b_Jet_leadTrackPt));
    RooDataHist p_n_Jet_leadTrackPt("p_n_Jet_leadTrackPt","p_n_Jet_leadTrackPt",RooArgList(Jet_leadTrackPt),Import(*data_p_n_Jet_leadTrackPt));
    RooDataHist r_p_Jet_leadTrackPt("r_p_Jet_leadTrackPt","r_p_Jet_leadTrackPt",RooArgList(Jet_leadTrackPt),Import(*r_data_p_Jet_leadTrackPt));

    RooDataHist p_b_Jet_leptonPt("p_b_Jet_leptonPt","p_b_Jet_leptonPt",RooArgList(Jet_leptonPt),Import(*data_p_b_Jet_leptonPt));
    RooDataHist p_n_Jet_leptonPt("p_n_Jet_leptonPt","p_n_Jet_leptonPt",RooArgList(Jet_leptonPt),Import(*data_p_n_Jet_leptonPt));
    RooDataHist r_p_Jet_leptonPt("r_p_Jet_leptonPt","r_p_Jet_leptonPt",RooArgList(Jet_leptonPt),Import(*r_data_p_Jet_leptonPt));

    RooDataHist p_b_Jet_leptonPtRel("p_b_Jet_leptonPtRel","p_b_Jet_leptonPtRel",RooArgList(Jet_leptonPtRel),Import(*data_p_b_Jet_leptonPtRel));
    RooDataHist p_n_Jet_leptonPtRel("p_n_Jet_leptonPtRel","p_n_Jet_leptonPtRel",RooArgList(Jet_leptonPtRel),Import(*data_p_n_Jet_leptonPtRel));
    RooDataHist r_p_Jet_leptonPtRel("r_p_Jet_leptonPtRel","r_p_Jet_leptonPtRel",RooArgList(Jet_leptonPtRel),Import(*r_data_p_Jet_leptonPtRel));

    RooDataHist p_b_Jet_leptonPdgId("p_b_Jet_leptonPdgId","p_b_Jet_leptonPdgId",RooArgList(Jet_leptonPdgId),Import(*data_p_b_Jet_leptonPdgId));
    RooDataHist p_n_Jet_leptonPdgId("p_n_Jet_leptonPdgId","p_n_Jet_leptonPdgId",RooArgList(Jet_leptonPdgId),Import(*data_p_n_Jet_leptonPdgId));
    RooDataHist r_p_Jet_leptonPdgId("r_p_Jet_leptonPdgId","r_p_Jet_leptonPdgId",RooArgList(Jet_leptonPdgId),Import(*r_data_p_Jet_leptonPdgId));

    RooDataHist p_b_Jet_leptonPtRelInv("p_b_Jet_leptonPtRelInv","p_b_Jet_leptonPtRelInv",RooArgList(Jet_leptonPtRelInv),Import(*data_p_b_Jet_leptonPtRelInv));
    RooDataHist p_n_Jet_leptonPtRelInv("p_n_Jet_leptonPtRelInv","p_n_Jet_leptonPtRelInv",RooArgList(Jet_leptonPtRelInv),Import(*data_p_n_Jet_leptonPtRelInv));
    RooDataHist r_p_Jet_leptonPtRelInv("r_p_Jet_leptonPtRelInv","r_p_Jet_leptonPtRelInv",RooArgList(Jet_leptonPtRelInv),Import(*r_data_p_Jet_leptonPtRelInv));

    RooDataHist p_b_Jet_vtxMass("p_b_Jet_vtxMass","p_b_Jet_vtxMass",RooArgList(Jet_vtxMass),Import(*data_p_b_Jet_vtxMass));
    RooDataHist p_n_Jet_vtxMass("p_n_Jet_vtxMass","p_n_Jet_vtxMass",RooArgList(Jet_vtxMass),Import(*data_p_n_Jet_vtxMass));
    RooDataHist r_p_Jet_vtxMass("r_p_Jet_vtxMass","r_p_Jet_vtxMass",RooArgList(Jet_vtxMass),Import(*r_data_p_Jet_vtxMass));

    RooDataHist p_b_Jet_vtxNtracks("p_b_Jet_vtxNtracks","p_b_Jet_vtxNtracks",RooArgList(Jet_vtxNtracks),Import(*data_p_b_Jet_vtxNtracks));
    RooDataHist p_n_Jet_vtxNtracks("p_n_Jet_vtxNtracks","p_n_Jet_vtxNtracks",RooArgList(Jet_vtxNtracks),Import(*data_p_n_Jet_vtxNtracks));
    RooDataHist r_p_Jet_vtxNtracks("r_p_Jet_vtxNtracks","r_p_Jet_vtxNtracks",RooArgList(Jet_vtxNtracks),Import(*r_data_p_Jet_vtxNtracks));

    RooDataHist p_b_Jet_vtxPt("p_b_Jet_vtxPt","p_b_Jet_vtxPt",RooArgList(Jet_vtxPt),Import(*data_p_b_Jet_vtxPt));
    RooDataHist p_n_Jet_vtxPt("p_n_Jet_vtxPt","p_n_Jet_vtxPt",RooArgList(Jet_vtxPt),Import(*data_p_n_Jet_vtxPt));
    RooDataHist r_p_Jet_vtxPt("r_p_Jet_vtxPt","r_p_Jet_vtxPt",RooArgList(Jet_vtxPt),Import(*r_data_p_Jet_vtxPt));

    RooDataHist p_b_Jet_vtx3DSig("p_b_Jet_vtx3DSig","p_b_Jet_vtx3DSig",RooArgList(Jet_vtx3DSig),Import(*data_p_b_Jet_vtx3DSig));
    RooDataHist p_n_Jet_vtx3DSig("p_n_Jet_vtx3DSig","p_n_Jet_vtx3DSig",RooArgList(Jet_vtx3DSig),Import(*data_p_n_Jet_vtx3DSig));
    RooDataHist r_p_Jet_vtx3DSig("r_p_Jet_vtx3DSig","r_p_Jet_vtx3DSig",RooArgList(Jet_vtx3DSig),Import(*r_data_p_Jet_vtx3DSig));

    RooDataHist p_b_Jet_vtx3DVal("p_b_Jet_vtx3DVal","p_b_Jet_vtx3DVal",RooArgList(Jet_vtx3DVal),Import(*data_p_b_Jet_vtx3DVal));
    RooDataHist p_n_Jet_vtx3DVal("p_n_Jet_vtx3DVal","p_n_Jet_vtx3DVal",RooArgList(Jet_vtx3DVal),Import(*data_p_n_Jet_vtx3DVal));
    RooDataHist r_p_Jet_vtx3DVal("r_p_Jet_vtx3DVal","r_p_Jet_vtx3DVal",RooArgList(Jet_vtx3DVal),Import(*r_data_p_Jet_vtx3DVal));

    RooDataHist p_b_Jet_ptd("p_b_Jet_ptd","p_b_Jet_ptd",RooArgList(Jet_ptd),Import(*data_p_b_Jet_ptd));
    RooDataHist p_n_Jet_ptd("p_n_Jet_ptd","p_n_Jet_ptd",RooArgList(Jet_ptd),Import(*data_p_n_Jet_ptd));
    RooDataHist r_p_Jet_ptd("r_p_Jet_ptd","r_p_Jet_ptd",RooArgList(Jet_ptd),Import(*r_data_p_Jet_ptd));

    RooDataHist p_b_rho("p_b_rho","p_b_rho",RooArgList(rho),Import(*data_p_b_rho));
    RooDataHist p_n_rho("p_n_rho","p_n_rho",RooArgList(rho),Import(*data_p_n_rho));
    RooDataHist r_p_rho("r_p_rho","r_p_rho",RooArgList(rho),Import(*r_data_p_rho));

    cout<<endl<<"Now it is x_dataset events = "<<data_set.sumEntries()<<endl;
    cout<<endl<<"Now it is x_r_dataset events = "<<r_data_set.sumEntries()<<endl;

//done only to see number of signal and bcg mc events in dataset based on Jet_mcFlavour

    RooDataSet* x_mc_sig_events = (RooDataSet*) data_set.reduce((var_set),"fabs(Jet_mcFlavour)==5");
    RooDataSet* x_mc_bg_events = (RooDataSet*) data_set.reduce((var_set),"fabs(Jet_mcFlavour)!=5");
  
    double max = std::numeric_limits<double>::infinity();

    RooRealVar yield_b("yield_b","yield_b",100,0,max);
    RooRealVar yield_n("yield_n","yield_n",100,0,max);

//Use any one of the pdfs sets below as pdf_b and pdf_n should be named as final pdfs:

    RooHistPdf pdf_b("pdf_b","",BDT_response,p_b_bdt);
    RooHistPdf pdf_n("pdf_n","",BDT_response,p_n_bdt);
/*
    RooHistPdf pdf_b("pdf_b","",BDT_response,p_b_60_120_bdt);
    RooHistPdf pdf_n("pdf_n","",BDT_response,p_n_60_120_bdt);

    RooHistPdf pdf_b("pdf_b","",BDT_response,p_b_120_bdt);
    RooHistPdf pdf_n("pdf_n","",BDT_response,p_n_120_bdt);
*/
    RooAddPdf totalPdf("totalPdf","totalPdf",RooArgList(pdf_b,pdf_n),RooArgList(yield_b,yield_n));
//    totalPdf.fitTo(data_set,Extended());

    RooRealVar r_yield_b("r_yield_b","r_yield_b",100,0,max);
    RooRealVar r_yield_n("r_yield_n","r_yield_n",100,0,max);

    RooAddPdf r_totalPdf("r_totalPdf","r_totalPdf",RooArgList(pdf_b,pdf_n),RooArgList(r_yield_b,r_yield_n));
//    r_totalPdf.fitTo(r_data_set,Extended());

    RooAddPdf* ptr_totalPdf = &totalPdf;
    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",data_set,ptr_totalPdf,RooArgList(yield_b,yield_n));

    RooAddPdf* r_ptr_totalPdf = &r_totalPdf;
    RooStats::SPlot* r_sData = new RooStats::SPlot("r_sData","An r_SPlot",r_data_set,r_ptr_totalPdf,RooArgList(r_yield_b,r_yield_n));
   
    RooFormulaVar rewei_b = RooFormulaVar("rewei_b","@0/@1*@2", RooArgList(r_yield_b,yield_b,*(data_set.get()->find("yield_b_sw"))));
    data_set.addColumn(rewei_b);

    RooFormulaVar rewei_n = RooFormulaVar("rewei_n","@0/@1*@2", RooArgList(r_yield_n,yield_n,*(data_set.get()->find("yield_n_sw"))));
    data_set.addColumn(rewei_n);


    RooDataSet* jjj = &data_set;
    RooDataSet* data_set_sig = new RooDataSet(jjj->GetName(),jjj->GetTitle(),jjj,*jjj->get(),"","rewei_b");
    RooDataSet* data_set_bg = new RooDataSet(jjj->GetName(),jjj->GetTitle(),jjj,*jjj->get(),"","rewei_n");

    RooDataSet* jj = &r_data_set;
    RooDataSet* r_data_set_sig = new RooDataSet(jj->GetName(),jj->GetTitle(),jj,*jj->get(),"","r_yield_b_sw");
    RooDataSet* r_data_set_bg = new RooDataSet(jj->GetName(),jj->GetTitle(),jj,*jj->get(),"","r_yield_n_sw");

    RooRealVar yield_b_sw("yield_b_sw","yield_b_sw",-5,5);
    RooRealVar yield_n_sw("yield_n_sw","yield_n_sw",-5,5);
    RooRealVar r_yield_b_sw("r_yield_b_sw","r_yield_b_sw",-5,5);
    RooRealVar r_yield_n_sw("r_yield_n_sw","r_yield_n_sw",-5,5);

    RooDataSet* w_ptr_data_set = (RooDataSet*) data_set.reduce(RooArgSet(yield_b_sw,yield_n_sw));
    RooDataSet w_data_set = *w_ptr_data_set;

    RooDataSet* w_ptr_r_data_set = (RooDataSet*) r_data_set.reduce(RooArgSet(r_yield_b_sw,r_yield_n_sw));
    RooDataSet w_r_data_set = *w_ptr_r_data_set;

    RooPlot frame_Jet_leptonDeltaR,frame_Jet_mass,frame_Jet_pt,frame_Jet_eta,frame_Jet_phi,frame_Jet_rawEnergy,frame_Jet_energyRing_dR0_neut,frame_Jet_energyRing_dR1_neut,frame_Jet_energyRing_dR2_neut,frame_Jet_energyRing_dR3_neut,frame_Jet_energyRing_dR4_neut,frame_Jet_energyRing_dR0_ch,frame_Jet_energyRing_dR1_ch,frame_Jet_energyRing_dR2_ch,frame_Jet_energyRing_dR3_ch,frame_Jet_energyRing_dR4_ch,frame_Jet_energyRing_dR0_em,frame_Jet_energyRing_dR1_em,frame_Jet_energyRing_dR2_em,frame_Jet_energyRing_dR3_em,frame_Jet_energyRing_dR4_em,frame_Jet_energyRing_dR0_mu,frame_Jet_energyRing_dR1_mu,frame_Jet_energyRing_dR2_mu,frame_Jet_energyRing_dR3_mu,frame_Jet_energyRing_dR4_mu,frame_Jet_numDaughters_pt03,frame_Jet_numberOfDaughters,frame_rho,frame_Jet_rawPt,frame_Jet_chHEF,frame_Jet_muEF,frame_Jet_neHEF,frame_Jet_chEmEF,frame_Jet_neEmEF,frame_Jet_leadTrackPt,frame_Jet_leptonPt,frame_Jet_leptonPtRel,frame_Jet_leptonPdgId,frame_Jet_leptonPtRelInv,frame_Jet_vtxMass,frame_Jet_vtxNtracks,frame_Jet_vtxPt,frame_Jet_vtx3DSig,frame_Jet_vtx3DVal,frame_Jet_ptd,frame_yield_b_sw,frame_yield_n_sw,frame_r_yield_b_sw,frame_r_yield_n_sw,frame_BDT_response;

    TCanvas c_Jet_leptonDeltaR,c_Jet_mass,c_Jet_pt,c_Jet_eta,c_Jet_phi,c_Jet_rawEnergy,c_Jet_energyRing_dR0_neut,c_Jet_energyRing_dR1_neut,c_Jet_energyRing_dR2_neut,c_Jet_energyRing_dR3_neut,c_Jet_energyRing_dR4_neut,c_Jet_energyRing_dR0_ch,c_Jet_energyRing_dR1_ch,c_Jet_energyRing_dR2_ch,c_Jet_energyRing_dR3_ch,c_Jet_energyRing_dR4_ch,c_Jet_energyRing_dR0_em,c_Jet_energyRing_dR1_em,c_Jet_energyRing_dR2_em,c_Jet_energyRing_dR3_em,c_Jet_energyRing_dR4_em,c_Jet_energyRing_dR0_mu,c_Jet_energyRing_dR1_mu,c_Jet_energyRing_dR2_mu,c_Jet_energyRing_dR3_mu,c_Jet_energyRing_dR4_mu,c_Jet_numDaughters_pt03,c_Jet_numberOfDaughters,c_rho,c_Jet_rawPt,c_Jet_chHEF,c_Jet_muEF,c_Jet_neHEF,c_Jet_chEmEF,c_Jet_neEmEF,c_Jet_leadTrackPt,c_Jet_leptonPt,c_Jet_leptonPtRel,c_Jet_leptonPdgId,c_Jet_leptonPtRelInv,c_Jet_vtxMass,c_Jet_vtxNtracks,c_Jet_vtxPt,c_Jet_vtx3DSig,c_Jet_vtx3DVal,c_Jet_ptd,c_yield_b_sw,c_yield_n_sw,c_r_yield_b_sw,c_r_yield_n_sw,c_BDT_response;

    TLegend leg_Jet_leptonDeltaR,leg_Jet_mass,leg_Jet_pt,leg_Jet_eta,leg_Jet_phi,leg_Jet_rawEnergy,leg_Jet_energyRing_dR0_neut,leg_Jet_energyRing_dR1_neut,leg_Jet_energyRing_dR2_neut,leg_Jet_energyRing_dR3_neut,leg_Jet_energyRing_dR4_neut,leg_Jet_energyRing_dR0_ch,leg_Jet_energyRing_dR1_ch,leg_Jet_energyRing_dR2_ch,leg_Jet_energyRing_dR3_ch,leg_Jet_energyRing_dR4_ch,leg_Jet_energyRing_dR0_em,leg_Jet_energyRing_dR1_em,leg_Jet_energyRing_dR2_em,leg_Jet_energyRing_dR3_em,leg_Jet_energyRing_dR4_em,leg_Jet_energyRing_dR0_mu,leg_Jet_energyRing_dR1_mu,leg_Jet_energyRing_dR2_mu,leg_Jet_energyRing_dR3_mu,leg_Jet_energyRing_dR4_mu,leg_Jet_numDaughters_pt03,leg_Jet_numberOfDaughters,leg_rho,leg_Jet_rawPt,leg_Jet_chHEF,leg_Jet_muEF,leg_Jet_neHEF,leg_Jet_chEmEF,leg_Jet_neEmEF,leg_Jet_leadTrackPt,leg_Jet_leptonPt,leg_Jet_leptonPtRel,leg_Jet_leptonPdgId,leg_Jet_leptonPtRelInv,leg_Jet_vtxMass,leg_Jet_vtxNtracks,leg_Jet_vtxPt,leg_Jet_vtx3DSig,leg_Jet_vtx3DVal,leg_Jet_ptd,leg_yield_b_sw,leg_yield_n_sw,leg_r_yield_b_sw,leg_r_yield_n_sw,leg_BDT_response;

    plot(data_p_b_Jet_chHEF,r_yield_b,&frame_Jet_chHEF,Jet_chHEF,0,1,100,TString("Jet_chHEF"),TString("signal"),totalPdf,data_set_sig,r_data_set_sig,p_b_Jet_chHEF,&c_Jet_chHEF,&leg_Jet_chHEF);
    plot(data_p_n_Jet_chHEF,r_yield_n,&frame_Jet_chHEF,Jet_chHEF,0,1,100,TString("Jet_chHEF"),TString("bcg"),totalPdf,data_set_bg,r_data_set_bg,p_n_Jet_chHEF,&c_Jet_chHEF,&leg_Jet_chHEF);

    plot(data_p_b_Jet_energyRing_dR0_ch,r_yield_b,&frame_Jet_energyRing_dR0_ch,Jet_energyRing_dR0_ch,0,1,100,TString("Jet_energyRing_dR0_ch"),TString("signal"),totalPdf,data_set_sig,r_data_set_sig,p_b_Jet_energyRing_dR0_ch,&c_Jet_energyRing_dR0_ch,&leg_Jet_energyRing_dR0_ch);
    plot(data_p_n_Jet_energyRing_dR0_ch,r_yield_n,&frame_Jet_energyRing_dR0_ch,Jet_energyRing_dR0_ch,0,1,100,TString("Jet_energyRing_dR0_ch"),TString("bcg"),totalPdf,data_set_bg,r_data_set_bg,p_n_Jet_energyRing_dR0_ch,&c_Jet_energyRing_dR0_ch,&leg_Jet_energyRing_dR0_ch);

    plot(data_p_b_Jet_chEmEF,r_yield_b,&frame_Jet_chEmEF,Jet_chEmEF,0,1,100,TString("Jet_chEmEF"),TString("signal"),totalPdf,data_set_sig,r_data_set_sig,p_b_Jet_chEmEF,&c_Jet_chEmEF,&leg_Jet_chEmEF);
    plot(data_p_n_Jet_chEmEF,r_yield_n,&frame_Jet_chEmEF,Jet_chEmEF,0,1,100,TString("Jet_chEmEF"),TString("background"),totalPdf,data_set_bg,r_data_set_bg,p_n_Jet_chEmEF,&c_Jet_chEmEF,&leg_Jet_chEmEF);


//Following are weight plots:
/*
    weight_plot(&frame_yield_b_sw,yield_b_sw,-5,5,100,TString("signal mc_weight"),&w_data_set,&c_yield_b_sw,&leg_yield_b_sw);
    weight_plot(&frame_yield_n_sw,yield_n_sw,-5,5,100,TString("bcg mc_weight"),&w_data_set,&c_yield_n_sw,&leg_yield_n_sw);

    weight_plot(&frame_r_yield_b_sw,r_yield_b_sw,-5,5,100,TString("signal r_weight"),&w_r_data_set,&c_r_yield_b_sw,&leg_r_yield_b_sw);
    weight_plot(&frame_r_yield_n_sw,r_yield_n_sw,-5,5,100,TString("bcg r_weight"),&w_r_data_set,&c_r_yield_n_sw,&leg_r_yield_n_sw);
*/
/*
    TCanvas* c0 = new TCanvas("c0","c0");
    RooPlot* frame0 = Jet_eta.frame(-3.2,3.2,100);  
    data_set_bg->plotOn(frame0);
    frame0->Draw();
    c0->Draw();
*/
/*
    TCanvas* c1 = new TCanvas("c1","c1");
    RooPlot* frame1 = Jet_eta.frame(-3.2,3.2,100);
    x_mc_bg_events->plotOn(frame1);
    frame1->Draw();
    c1->Draw();
*/
}

    void plot(TH1F* data_p_bn_var_name,RooRealVar r_yield_bn,RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,TString comp,RooAddPdf totalPdf_var_name,RooDataSet* data_set,RooDataSet* r_data_set,RooDataHist p_bn_var_name,TCanvas* TCanvas_var_name,TLegend* leg_var_name)
{
    data_p_bn_var_name->Scale(r_yield_bn.getVal()/data_p_bn_var_name->Integral());
  
    frame_var_name = var_name.frame(x_low, x_high, bins);
    frame_var_name->SetTitle(title + " (" + comp + ")");
    frame_var_name->GetXaxis()->CenterTitle();
    totalPdf_var_name.paramOn(frame_var_name,Format("NEU",AutoPrecision(2)),Layout(0.1,0.4,0.9));
    frame_var_name->getAttText()->SetTextSize(0.03);
    frame_var_name->getAttText()->SetTextColor(kRed+3);
    data_set->plotOn(frame_var_name,MarkerSize(0.7),MarkerColor(kRed),LineColor(kRed),Name("sweighted MC"));
    r_data_set->plotOn(frame_var_name,MarkerSize(0.7),MarkerColor(kGreen+2),LineColor(kGreen+2),Name("sweighted data"));
    p_bn_var_name.plotOn(frame_var_name,MarkerSize(0.7),MarkerColor(kBlue+2),LineColor(kBlue+2),Name("unweighted MC"));

    TCanvas_var_name = new TCanvas(title+ "(" + comp + ")",title+ "(" + comp + ")");
    string variable_name = std::string(var_name.GetName());
    double bin_width = frame_var_name->GetXaxis()->GetBinWidth(0);
    std::string str = std::to_string(bin_width);
    str.erase (str.find_last_not_of('0') + 1, std::string::npos );
    TString yaxis = "log(events)/" + str;
/*    if (variable_name.find("energyRing")==4 || variable_name.find("rawPt")==4 || variable_name.find("muEF")==4 ||variable_name.find("neHEF")==4 || variable_name.find("chEmEF")==4 || variable_name.find("neEmEF")==4 || variable_name.find("leadTrackPt")==4 || variable_name.find("leptonPt")==4 || variable_name.find("vtxMass")==4 || variable_name.find("vtx3DSig")==4 ||variable_name.find("vtx3DVal")==4 ||variable_name.find("ptd")==4 || variable_name.find("vtxPt")==4)
    {
        TCanvas_var_name->SetLogy();
        frame_var_name->GetYaxis()->SetTitle(yaxis);
    }
*/    frame_var_name->GetYaxis()->SetTitleOffset(1.5);
    frame_var_name->GetXaxis()->SetTitleOffset(1.1);
    leg_var_name = new TLegend(0.6,0.6,0.8,0.8);
    TObject* sMC_leg = (RooPlot*)frame_var_name->findObject("sweighted MC");
    TObject* sdata_leg = (RooPlot*)frame_var_name->findObject("sweighted data");
    TObject* MC_leg = (RooPlot*)frame_var_name->findObject("unweighted MC");
    leg_var_name->SetFillColor(kWhite);
    leg_var_name->SetTextSize(0.03);
    leg_var_name->AddEntry(sMC_leg,"sweighted MC","LEP");
    leg_var_name->AddEntry(sdata_leg,"sweighted data","LEP");
    leg_var_name->AddEntry(MC_leg,"unweighted MC","LEP");
    leg_var_name->SetBorderSize(0);
    frame_var_name->Draw();
    leg_var_name->Draw();
    TCanvas_var_name->Draw();

}

    void weight_plot(RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,RooDataSet* data_set,TCanvas* TCanvas_var_name,TLegend* leg_var_name)
{
    frame_var_name = var_name.frame(x_low, x_high, bins);
    frame_var_name->SetTitle(title);
    frame_var_name->GetXaxis()->CenterTitle();
    data_set->plotOn(frame_var_name,MarkerSize(0.7),MarkerColor(kRed),LineColor(kRed),Name("sweighted MC"));

    TCanvas_var_name = new TCanvas(title,title);
    frame_var_name->GetYaxis()->SetTitleOffset(1.5);
    frame_var_name->GetXaxis()->SetTitleOffset(1.1);
    leg_var_name = new TLegend(0.6,0.6,0.8,0.8);
    TObject* sMC_leg = (RooPlot*)frame_var_name->findObject("sweighted MC");
    leg_var_name->SetFillColor(kWhite);
    leg_var_name->SetTextSize(0.03);
    leg_var_name->AddEntry(sMC_leg,title,"LEP");
    leg_var_name->SetBorderSize(0);
    frame_var_name->Draw();
    leg_var_name->Draw();
    TCanvas_var_name->Draw();
}









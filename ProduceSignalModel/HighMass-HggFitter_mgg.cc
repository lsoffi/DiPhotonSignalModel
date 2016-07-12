#include "CMS_lumi.C"
//#include "ModelConfig.h"
#include "RooVoigtian.h"
#include "RooRealVar.h"
#include "RooBinning.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvolution.h"
#include "RooHistFunc.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "RooGenericPdf.h"
#include <fstream>
#include "RooDCBShape.h"
//#include "RooCurve.h"
//#include "RooNumConvPdf.h"
//#include "RooHist.h"

using namespace RooFit;
using namespace RooStats;
using namespace RooFit;
using namespace RooStats ;
static const Int_t NCAT = 4;  // chiara
Int_t MINmass= 300;
Int_t MAXmass= 6000;
std::string filepostfix="";
double signalscaler=1.00;

void sigModelResponseFcnFit(RooWorkspace* w,Float_t mass,TString coupling, TString whichRel);
void sigModelGenFcnFit(RooWorkspace* w,Float_t mass,TString coupling, TString whichRel);
void sigModelShapeFcnFit(RooWorkspace* w001,RooWorkspace* w,Float_t mass,TString coupling,TString whichRel);
RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov );
void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, TString coupling,TString whichRel);
void runAllfits();
double computePdfFHWM(RooDCBShape pdf,RooRealVar roobs, double MH);
// Definition of the variables in the input ntuple
RooArgSet* defineVariables() {

  RooRealVar* mgg        = new RooRealVar("mgg",        "M(gg)",       MINmass, MAXmass, "GeV");
  RooRealVar* mggGen     = new RooRealVar("mggGen",     "M(gg) gen",   MINmass, MAXmass, "GeV");
  RooRealVar* eventClass = new RooRealVar("eventClass", "eventClass",    -10,      10,   "");
  RooRealVar* weight     = new RooRealVar("puweight",     "weightings",      0,     1000,  "");   

  RooArgSet* ntplVars = new RooArgSet(*mgg, *mggGen, *eventClass);                  
  
  return ntplVars;
}

void SetConstantParams(const RooArgSet* params) {
  
  std::cout << std::endl; std::cout << "Entering SetConstantParams" << std::endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}


TPaveText* get_labelcms( int legendquadrant = 0 , std::string year="2012", bool sim=false) {

  if( legendquadrant!=0 && legendquadrant!=1 && legendquadrant!=2 && legendquadrant!=3 ) {
    std::cout << "warning! legend quadrant '" << legendquadrant << "' not yet implemented for cms label. using 2." << std::endl;
    legendquadrant = 2;
  }

  float x1, y1, x2, y2;
  if( legendquadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendquadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;

  } else if( legendquadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
  }
 
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brndc" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendquadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
  
    std::string lefttext;
   
     
    if (sim)  lefttext = "cms simulation"; 
    else {
     lefttext = "cms preliminary, 19.5 fb^{-1}";
    }
    cmslabel->AddText(lefttext.c_str());
    return cmslabel;

}




TPaveText* get_labelsqrt( int legendquadrant ) {

  if( legendquadrant!=0 && legendquadrant!=1 && legendquadrant!=2 && legendquadrant!=3 ) {
    std::cout << "warning! legend quadrant '" << legendquadrant << "' not yet implemented for sqrt label. using 2." << std::endl;
    legendquadrant = 2;
  }


  float x1, y1, x2, y2;
  if( legendquadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendquadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendquadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendquadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brndc");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); 
  label_sqrt->AddText("#sqrt{s} = 8 tev");
  return label_sqrt;

}



// loading signal data and making roodatasets
void AddSigData(RooWorkspace* w, Float_t mass, TString coupling, TString whichRel) {


  Int_t ncat = NCAT;
  
  // Variables
  RooArgSet* ntplVars = defineVariables();

  // -------------------------  
  // Files
  int iMass = abs(mass);   
  TString inDir = "./";
  TString mainCut1 = TString::Format("mgg>=300 && mgg<=6000 && mggGen>=300 && mggGen<=6000");  
  TChain* sigTreeGen = new TChain();
  cout << "reading file gen " << "root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/SigModelGen/FormSigModGenOnly_genGenIso_kpl"+coupling+TString::Format("_M%d.root/DiPhotonTree", iMass) <<endl; 
    
  sigTreeGen->Add("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/SigModelGen/FormSigModGenOnly_genGenIso_kpl"+coupling+TString::Format("_M%d.root/DiPhotonTree", iMass));
  sigTreeGen->SetTitle("sigTreeGen");
  sigTreeGen->SetName("sigTreeGen");

  cout << "SigTreeGen: " << endl;
  //sigTreeGen->Print("V");
  cout << " nX for SigTreeGen:  " << sigTreeGen->GetEntries() << endl;
  cout << endl;
 
  // common preselection cut on the reduced mass
    TString mainCut = TString::Format("(mgg-mggGen)>-600.&&(mgg-mggGen)<600."); 
    // -------------------------
    bool ok80X =  coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000));
    bool ok76X_38T = coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass ==4000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000));
    bool ok76X_0T = coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass ==3000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 750 ||  mass == 3000 ));
    if((whichRel== "80X" && ok80X )||(whichRel== "76X_38T" && ok76X_38T )||(whichRel== "76X_0T" && ok76X_0T )){
      TChain* sigTreeK = new TChain();
      RooDataSet* sigWeightedK;
    cout << "reading file " <<endl; //<< TString::Format("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl"+coupling+"_M%d.root/DiPhotonTree", iMass) << endl;
    if((mass==1000&&whichRel=="80X") || ((mass==1000||mass==1250 || mass==1750 || mass==2250 || mass==2500 || mass==2750)&&whichRel=="76X_38T")||  ((mass==1000||mass==1250 || mass==1500 || mass==1750 || mass==2250 || mass==2500 || mass==2750|| mass==3000)&&whichRel=="76X_0T"))sigTreeK->Add("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl01"+TString::Format("_M%d.root/DiPhotonTree", iMass));
    else sigTreeK->Add("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl"+coupling+TString::Format("_M%d.root/DiPhotonTree", iMass));
    sigTreeK->SetTitle("sigTreeK");
    sigTreeK->SetName("sigTreeK");
    sigWeightedK = new RooDataSet("sigWeightedK","datasetK",sigTreeK,*ntplVars,mainCut1, "puweight");
    
    cout << "preparing dataset with observable mgg correct k" << endl;
    RooDataSet* signalK[NCAT];
    std::cout<<mass<<std::endl;
    RooDataSet* signalAllK;
    std::cout<<coupling<<std::endl;
    
    for (int c=0; c<ncat; ++c) {
      if (c==0) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
      if (c==1) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));
      if (c==2) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==2"));
      if (c==3) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==3"));
      
      TString myCut;
      w->import(*signalK[c],Rename(TString::Format("SigWeightK_cat%d",c)));
      cout << "cat " << c << ", signalK[c]: " << endl;
      signalK[c]->Print("V");
      cout << "---- for category " << c << ", nX for signal[c]:  " << signalK[c]->sumEntries() << endl; 
      cout << endl;
    }
    
    // Create full weighted signal data set without categorization
    signalAllK = (RooDataSet*) sigWeightedK->reduce(RooArgList(*w->var("mgg")),mainCut);
    w->import(*signalAllK, Rename("SigWeightK"));
    cout << "now signalAllK" << endl;
    signalAllK->Print("V");
    cout << "---- nX for signalAll:  " << signalAllK->sumEntries() << endl; 
    cout << endl;
  }





  if(coupling=="001"){
     
    TChain* sigTree = new TChain();
    cout << "reading file " << "root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl"+coupling+TString::Format("_M%d.root/DiPhotonTree", iMass) << endl;
    if((mass==1000&&whichRel=="80X") || ((mass==1000||mass==1250 || mass==1750 || mass==2250 || mass==2500 || mass==2750)&&whichRel=="76X_38T") ||  ((mass==1000||mass==1250 || mass==1500 || mass==1750 || mass==2250 || mass==2500 || mass==2750|| mass==3000)&&whichRel=="76X_0T"))sigTree->Add("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl01_"+TString::Format("M%d.root/DiPhotonTree", iMass));
    else sigTree->Add("root://eoscms//eos/cms/store/group/phys_higgs/soffi/Pippone/RecoModel"+whichRel+"/FormSigMod_cicGenIso_kpl001_"+TString::Format("M%d.root/DiPhotonTree", iMass));
    sigTree->SetTitle("sigTree");
    sigTree->SetName("sigTree");
    cout <<  "SigTree: " << endl;
    //  sigTree->Print("V");
    cout << " nX for SigTree:  " << sigTree->GetEntries() << endl;
    cout << endl;   
    // -------------------------
    // -------------------------
    // common preselection cut on mgg and mggGen
    
    RooDataSet sigWeighted("sigWeighted","dataset",sigTree,*ntplVars,mainCut1, "puweight"); 
    // -------------------------
    // reduced mass
    RooFormulaVar *massReduced_formula = new RooFormulaVar("massReduced_formula","","@0-@1",RooArgList(*w->var("mgg"),*w->var("mggGen")));
    RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
    massReduced->SetName("massReduced");
    massReduced->SetTitle("massReduced");
    w->import(*massReduced);  
    massReduced->setRange(-600., 600.);
    
    // split in categories, wrt mgg - this is the dataset to be used for the convolution
    cout << endl;
    cout << "preparing dataset with observable mgg k = 001" << endl;
    RooDataSet* signal[NCAT];
    
    for (int c=0; c<ncat; ++c) {
      if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
      if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));
      if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==2"));
      if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==3"));
      
      TString myCut;
      w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
      cout << "cat " << c << ", signal[c]: " << endl;
      signal[c]->Print("V");
      cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
      cout << endl;
    }
    
    // Create full weighted signal data set without categorization
    RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(RooArgList(*w->var("mgg")),mainCut);
    w->import(*signalAll, Rename("SigWeight"));
    cout << "now signalAll" << endl;
    signalAll->Print("V");
    cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
    cout << endl;
    
   
    
    bool wantResponse = true;
    // -------------------------
    // split in categories, wrt massReduced - to study the detector response
    if (wantResponse) {
      cout << endl;
      cout << endl;
      cout << "preparing dataset with observable massReduced" << endl;
      RooDataSet* signalR[NCAT];
      for (int c=0; c<ncat; ++c) {
	if (c==0) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==0"));
	if (c==1) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==1"));
	if (c==2) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==2"));
	if (c==3) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==3"));
	
	TString myCut;
	w->import(*signalR[c],Rename(TString::Format("SigWeightReduced_cat%d",c)));
      }
      cout << endl;
      // Create full weighted signal data set without categorization
      RooDataSet* signalRAll = (RooDataSet*) sigWeighted.reduce(RooArgList(*w->var("massReduced")),mainCut);
      w->import(*signalRAll, Rename("SigWeightReduced"));
      cout << "now signalRAll" << endl;
      signalRAll->Print("V");
      cout << "---- nX for signalRAll:  " << signalRAll->sumEntries() << endl; 
      cout << endl;
      
    
    }
  }
  bool wantGenLevel=true;
  RooDataSet sigWeightedGen("sigWeightedGen","datasetGen",sigTreeGen,*ntplVars,mainCut1, "puweight");   
  // -------------------------
  // split in categories, wrt genMass - to study the theory width
  if (wantGenLevel) { 
    cout << endl;
    cout << endl;
    cout << "preparing dataset with observable mggGen, no split in categories since they're all the same" << endl;
    RooDataSet* signalG[NCAT];
    for (int c=0; c<ncat; ++c) {
      if (c==0) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==0"));
      if (c==1) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==1"));
      if (c==2) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==2"));
      if (c==3) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==3"));
      
      TString myCut;
      w->import(*signalG[c],Rename(TString::Format("SigWeightGen_cat%d",c)));
    }
    RooDataSet* signalGAll = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"));
    w->import(*signalGAll, Rename("SigWeightGen"));
    cout << "now signalGAll" << endl;
    signalGAll->Print("V");
    cout << "---- nX for signalGAll:  " << signalGAll->sumEntries() << endl; 
    cout << endl;
    cout << endl;
  }
  //w->writeToFile("data_001_5000.root");
  cout << "workspace summary" << endl;
  w->Print();

}


// To run the analysis. Pdfs obtained here
void runfits(const Float_t mass=750, TString coupling="001", TString whichRel="80X") {
  gSystem->Load(".L RooDCBShape_cxx.so");
  //******************************************************************//
  //  Steps:
  //     - create signal and background data sets 
  //     - make and fit signal and background  models 
  //     - write signal and background workspaces in root files
  //     - write data card
  //*******************************************************************//
  
  TString fileBaseName("HighMassGG");    
  TString fileBkgName("HighMassGG.inputbkg");
  HLFactory hlf("HLFactory", "HighMassGG.rs", false);
  RooWorkspace* w = hlf.GetWs();
  int iMass = (int) abs(mass);   
  // range for the variables
  w->var("mgg")->setMin(MINmass);
  w->var("mgg")->setMax(MAXmass);
  w->var("mggGen")->setMin(MINmass);
  w->var("mggGen")->setMax(MAXmass);
  w->Print("V");
  
  cout << endl; 
  cout << "Now add signal data" << endl;
  std::cout<<"------> "<<coupling<<std::endl;
  AddSigData(w, mass, coupling, whichRel);
  
  TCanvas* canv = new TCanvas("canv","c",1);
  canv->cd();
  RooPlot* p = w->var("mgg")->frame(Range(300,3000), Bins(260));

  bool makeWs=false; 
  if(makeWs){
    if(coupling=="001")sigModelResponseFcnFit(w, mass, coupling, whichRel);
    sigModelGenFcnFit(w, mass, coupling, whichRel);
    w->Print("V");
    TFile* f = new TFile(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+".root", iMass), "RECREATE");
    f->cd();
    w->Write();
    f->Write();
    f->Close();
  }else{
    TFile* f001 = TFile::Open(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k001.root", iMass));
    TFile* f = TFile::Open(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+".root", iMass));

    RooWorkspace* ws001 = (RooWorkspace*) f001->Get("HLFactory_ws"); 
    RooWorkspace* ws = (RooWorkspace*) f->Get("HLFactory_ws"); 
    sigModelShapeFcnFit(ws001,ws, mass, coupling, whichRel);
    asimovDatasetFcnFit(ws, mass, coupling, whichRel);
    TFile* fout = new TFile(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+"_final.root", iMass), "RECREATE");
    fout->cd();
    ws->Write();
    fout->Write();
    fout->Close();
    }
 
}




void runAllFits(TString whichRel){
  double masses[11] = {500,750, 1000,1250,1750,2000,2250,2500,2750,3000,4000};// k001: 3250, 3500
  for(int  m = 0; m< 11; m++)runfits(masses[m], "001", whichRel);
  for(int  m = 0; m< 11; m++)runfits(masses[m], "01", whichRel);
  for(int  m = 0; m< 11; m++)runfits(masses[m], "02", whichRel);
 
}



void sigModelResponseFcnFit(RooWorkspace* w, Float_t mass, TString coupling, TString whichRel) {
  int ncat = NCAT;
  RooDataSet* signal;
  RooCBShape* responsecbpos[NCAT+1];
  RooCBShape* responsecbneg[NCAT+1];
  RooDCBShape* responseadd[NCAT+1];
  int iMass = (int)  abs(mass);
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelcms(0, "2012", true);
  TPaveText* label_sqrt = get_labelsqrt(0);
  
  for(int c = 0; c<ncat+1; c++){
    if(c==2||c==3)continue;
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    TString myCut;
    if(c==0||c==1||c==2||c==3)signal = (RooDataSet*)w->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)signal = (RooDataSet*)w->data("SigWeightReduced");
   //cb pos                                                                                                                     
   RooFormulaVar cbpos_mean(TString::Format("reducedmass_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("reducedmass_sig_mean_cat%d",c)));
   RooFormulaVar cbpos_sigma(TString::Format("reducedmass_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)));
   RooFormulaVar cbpos_alphacb(TString::Format("reducedmass_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbpos_cat%d",c)));
   RooFormulaVar cbpos_n(TString::Format("reducedmass_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_npos_cat%d",c)));     
   //cb neg
   RooFormulaVar cbneg_n(TString::Format("reducedmass_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_nneg_cat%d",c)));
   RooFormulaVar cbneg_alphacb(TString::Format("reducedmass_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbneg_cat%d",c)));   
   
   responseadd[c]= new  RooDCBShape(TString::Format("responseaddpdf_cat%d",c),TString::Format("responseaddpdf_cat%d",c) , *w->var("massReduced"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;
   w->import(*responseadd[c]);
   bool fixPar = true;
   if(mass==2000 && fixPar){
     w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-3.);
     w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(25);
     w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(1);
     w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1);
     w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(5);
     w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(5);
   }else if(mass==1500 && fixPar){
     w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-12.);
     w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(25);
     w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(-2);
     w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1);
     w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(3);
     w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(7);
   }
   RooFitResult* fitresults = (RooFitResult* ) responseadd[c]->fitTo(*signal,SumW2Error(kTRUE), Range(-600, 600), RooFit::Save(kTRUE));
   std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
   fitresults->Print("V");
   
   RooRealVar* massreduced = (RooRealVar*) w->var("massReduced");
   RooPlot* plotg = massreduced->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    signal->plotOn(plotg);
   

    responseadd[c]->plotOn(plotg, LineColor(kBlue));
    // responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    //responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotg->GetYaxis()->SetRangeUser(0.01,plotg->GetMaximum()*10 );
    plotg->GetXaxis()->SetTitle("#Delta m [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    // legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    // legmc->SetHeader( "#splitline{m_{x}=150 gev}{#splitline{}{class 0}}");
    legmc->AddEntry(plotg->getObject(0),"Simulation","lp");    
   
    legmc->AddEntry(plotg->getObject(1),"Double-Sided CB ","l");
    // legmc->AddEntry(plotg->getObject(2),"cb 1","l");   
    //legmc->AddEntry(plotg->getObject(3),"cb 2","l");   

    
   
    plotg->Draw();
    
    // lat->draw("same");
    // latex->draw("same");
    legmc->Draw("same");
    //    int ipos=11 ;
    //cms_lumi( c1,true,ipos );
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    c1->SetLogy(0);


    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/responsefcnfitcbcb_cat%d_M%d_k"+coupling+".png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/responsefcnfitcbcb_cat%d_M%d_k"+coupling+".pdf",c,iMass)); 
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/responsefcnfitcbcb_cat%d_M%d_k"+coupling+"_log.png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/responsefcnfitcbcb_cat%d_M%d_k"+coupling+"_log.pdf",c,iMass)); 
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    c1->SetLogy(0);  
	
    w->defineSet(TString::Format("responseaddpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c)),
									  *w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c)),
									  *w->var(TString::Format("reducedmass_sig_npos_cat%d",c)),
									  *w->var(TString::Format("reducedmass_sig_nneg_cat%d",c)),	   
									  *w->var(TString::Format("reducedmass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("responseaddpdfparam_cat%d",c)));
  }
  

 

}



void sigModelGenFcnFit(RooWorkspace* w, Float_t mass, TString coupling, TString whichRel) {
  int ncat = NCAT;
  RooDataSet* signal;
  //  RooBreitWigner* bw[NCAT+1];
  RooGenericPdf* bw[NCAT+1];
  RooCBShape* mgencbpos[NCAT+1];
  RooCBShape* mgencbneg[NCAT+1];
  RooDCBShape* mgenadd[NCAT+1];
  int iMass =  (int) abs(mass);
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelcms(0, "2012", true);
  TPaveText* label_sqrt = get_labelsqrt(0);
  
  for(int c = 4; c<ncat+1; c++){
    
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    TString myCut;
    signal = (RooDataSet*)w->data("SigWeightGen");
    w->var(TString::Format("mgen_sig_mean_cat%d",c))->setVal(mass);
    w->var(TString::Format("mgen_sig_mean_cat%d",c))->setRange(mass*0.8, mass*1.2);
    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean(TString::Format("mgen_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("mgen_sig_mean_cat%d",c)));
    RooFormulaVar cbpos_sigma(TString::Format("mgen_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("mgen_sig_sigma_cat%d",c)));
    RooFormulaVar cbpos_alphacb(TString::Format("mgen_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n(TString::Format("mgen_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_npos_cat%d",c)));    
    //cb neg
    RooFormulaVar cbneg_n(TString::Format("mgen_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb(TString::Format("mgen_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbneg_cat%d",c)));    
    mgenadd[c]=  new  RooDCBShape(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c) , *w->var("mggGen"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;
    w->import(*mgenadd[c]);
    std::cout<<mass<<" "<<signal->sumEntries()<<std::endl;
     std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
     RooFitResult* fitresults = (RooFitResult* ) mgenadd[c]->fitTo(*signal,SumW2Error(kTRUE), Range(mass*0.8, mass*1.2), RooFit::Save(kTRUE));   
     fitresults->Print("V");

  
   RooPlot* plotg = w->var("mggGen")->frame(Range(mass*0.8,mass*1.2),Title("mass generated"),Bins(120));
   signal->plotOn(plotg);
   mgenadd[c]->plotOn(plotg, LineColor(kBlue));
   // mgenadd[c]->plotOn(plotg,Components(TString::Format("mgencbneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
   //mgenadd[c]->plotOn(plotg,Components(TString::Format("mgencbpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
   

   
   plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}^{gen} [GeV]");
   plotg->GetXaxis()->SetTitleFont(42);
   plotg->GetYaxis()->SetTitle("arbitrary scale");
   plotg->GetYaxis()->SetTitleFont(42);
   plotg->GetYaxis()->SetTitleSize(0.04);
   TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
   legmc->SetTextFont(42);
   legmc->SetBorderSize(0);
   legmc->SetFillStyle(0);
   legmc->AddEntry(plotg->getObject(0),"Gen Simulation","lp");    
   legmc->AddEntry(plotg->getObject(1),"Double-Sided CB ","l");
   //  legmc->AddEntry(plotg->getObject(2),"cb 1","l");   
   // legmc->AddEntry(plotg->getObject(3),"cb 2","l");   
 
   // bw[c]->plotOn(plotg, LineColor(kBlue));
   plotg->GetYaxis()->SetRangeUser(0.01,plotg->GetMaximum()*1.2 );
   plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
   plotg->GetXaxis()->SetTitleFont(42);
   // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
   plotg->GetYaxis()->SetTitle("arbitrary scale");
   plotg->GetYaxis()->SetTitleFont(42);
   plotg->GetYaxis()->SetTitleSize(0.04);
   plotg->Draw();
   legmc->Draw("same");
  c1->SetLogy(0);


    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mgen_cat%d_M%d_k"+coupling+".png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mgen_cat%d_M%d_k"+coupling+".pdf",c,iMass)); 
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mgen_cat%d_M%d_k"+coupling+"_log.png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mgen_cat%d_M%d_k"+coupling+"_log.pdf",c,iMass)); 
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    c1->SetLogy(0);  
	
    w->defineSet(TString::Format("mgenpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("mgen_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("mgen_sig_alphacbpos_cat%d",c)),
									  *w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c)),
									  *w->var(TString::Format("mgen_sig_npos_cat%d",c)),
									  *w->var(TString::Format("mgen_sig_nneg_cat%d",c)),	   
									  *w->var(TString::Format("mgen_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("mgenpdfparam_cat%d",c)));
							 }
		 
		 
		 
		 
}



void sigModelShapeFcnFit(RooWorkspace* w001, RooWorkspace* w,Float_t mass, TString coupling, TString whichRel) {
  int ncat = NCAT;
  RooDataSet* signalK;
  //  RooBreitWigner* bw[NCAT+1];
  RooDCBShape* reducedmass[NCAT+1];
  RooCBShape*  mgencbpos[NCAT+1] ;
  RooCBShape*  mgencbneg[NCAT+1] ;
  RooCBShape*  reducedmasscbpos[NCAT+1] ;
  RooCBShape*  reducedmasscbneg[NCAT+1] ;
  RooDCBShape* mgen[NCAT+1];
  RooFFTConvPdf* conv[NCAT+1];
  int iMass =  (int) abs(mass);
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelcms(0, "2012", true);
  TPaveText* label_sqrt = get_labelsqrt(0);
  
  for(int c = 0; c<ncat+1; c++){
    if(c==2||c==3)continue;
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TString myCut;
    
   

  
    //pickup the functions from the ws
    double massMin,massMax;
    if(coupling=="001"){
      massMin=0.8*mass;
      massMax=1.2*mass;
    }else if(coupling=="01"){
      massMin=0.6*mass;
      massMax=1.4*mass;
    }else if(coupling=="02"){
      massMin=0.4*mass;
      massMax=1.6*mass;
    }
    RooRealVar* var = new RooRealVar("var", "var", massMin, massMax);
    RooRealVar* MH = new RooRealVar("MH", "MH",mass);
    w->import(*MH);
    MH->setConstant();
    //response
    double mh =100.;
    //gen shape
    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean_mgen(TString::Format("cbpos_sig_mean_cat4_mgen",c),"", "@0", *w->var(TString::Format("mgen_sig_mean_cat4",c)));
    RooFormulaVar cbpos_sigma_mgen(TString::Format("cbpos_sig_sigma_cat4_mgen",c), "", "sqrt(@0*@0)", *w->var(TString::Format("mgen_sig_sigma_cat4",c)));
    RooFormulaVar cbpos_alphacb_mgen(TString::Format("cbpos_sig_alphacb_cat4_mgen",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbpos_cat4",c)));
    RooFormulaVar cbpos_n_mgen(TString::Format("cbpos_sig_n_cat4_mgen",c),"", "@0", *w->var( TString::Format("mgen_sig_npos_cat4",c)));
    //cb neg
    RooFormulaVar cbneg_n_mgen(TString::Format("cbneg_sig_n_cat4_mgen",c),"", "@0", *w->var( TString::Format("mgen_sig_nneg_cat4",c)));
    RooFormulaVar cbneg_alphacb_mgen(TString::Format("cbneg_sig_alphacb_cat4_mgen",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbneg_cat4",c)));
    mgen[c]= new  RooDCBShape(TString::Format("mgen_cat%d",c),TString::Format("mgen_cat%d",c) , *var, cbpos_mean_mgen, cbpos_sigma_mgen,  cbneg_alphacb_mgen, cbpos_alphacb_mgen,  cbneg_n_mgen,cbpos_n_mgen) ;
    //response shape
    //cb pos                                  
    double scale=1;
    RooFormulaVar cbpos_mean_reducedmass(TString::Format("cbpos_sig_mean_cat%d_reducedmass",c),"","@0",RooArgList(*w001->var(TString::Format("reducedmass_sig_mean_cat%d",c)),*w->var("MH"))  );
    RooFormulaVar cbpos_sigma_reducedmass(TString::Format("cbpos_sig_sigma_cat%d_reducedmass",c), "", "sqrt(@0*@0)", RooArgList(*w001->var(TString::Format("reducedmass_sig_sigma_cat%d",c)),*w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    RooFormulaVar cbpos_alphacb_reducedmass(TString::Format("cbpos_sig_alphacb_cat%d_reducedmass",c),"", "@0", *w001->var( TString::Format("reducedmass_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n_reducedmass(TString::Format("cbpos_sig_n_cat%d_reducedmass",c),"", "@0", *w001->var( TString::Format("reducedmass_sig_npos_cat%d",c)));
    //cb neg
    RooFormulaVar cbneg_n_reducedmass(TString::Format("cbneg_sig_n_cat%d_reducedmass",c),"", "@0", *w001->var( TString::Format("reducedmass_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb_reducedmass(TString::Format("cbneg_sig_alphacb_cat%d_reducedmass",c),"", "@0", *w001->var( TString::Format("reducedmass_sig_alphacbneg_cat%d",c)));
    reducedmass[c]= new  RooDCBShape(TString::Format("addpdf_cat%d_reducedmass",c), TString::Format("addpdf_cat%d_reducedmass",c), *var, cbpos_mean_reducedmass, cbpos_sigma_reducedmass,  cbneg_alphacb_reducedmass, cbpos_alphacb_reducedmass,  cbneg_n_reducedmass,cbpos_n_reducedmass) ;
    
    //convolution
    var->setBins(4000, "cache");
    conv[c] = new RooFFTConvPdf(TString::Format("convpdf_cat%d",c),TString::Format("convpdf_cat%d",c), *var,*mgen[c],*reducedmass[c]);
    conv[c]->Print("V");
    std::cout<<"******************************************************************************************"<<std::endl;
    RooArgSet* params = conv[c]->getParameters(*w->var("mgg")) ;
    params->Print("v") ;
   
    //asimov dataset
    double min;
    double max;
    min=mass*0.8;
    max=mass*1.2;
    if(massMin<300)massMin=300;
    if(massMax>6000)massMax=6000;
    
    RooBinning binning(160, massMin, massMax, "binning");
    TH1F* dummy = new TH1F("dummy", "dummy", 160, massMin, massMax);
    RooDataHist* datahist_asimov = new RooDataHist("datahist_asimov","datahist_asimov",*var, dummy, 1);
    RooDataHist* datahist_asimov_clone = new RooDataHist("datahist_asimov","datahist_asimov",*w->var("mgg"), dummy, 1);
      
    //  RooDataHist* datahist_asimov; //= mgen[c]->generateBinned(RooArgSet(*var),RooFit::NumEvents(signal->numEntries()), RooFit::Asimov());
    //datahist_asimov = signal->binnedClone();
    throwAsimov(10000, conv[c], var,datahist_asimov); 
   
    if(c<4)    datahist_asimov->SetName(TString::Format("signal_asimov_cat%d",c));
    if(c==4)    datahist_asimov->SetName(TString::Format("signal_asimov",c));
    std::cout<<"--------------------------------------------------- "<<datahist_asimov->sumEntries()<<std::endl;
    w->import(*datahist_asimov);

    bool doPlot=false;
    bool ok80X =  coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000));
    bool ok76X_38T = coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass ==4000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 2000 ||  mass == 3000 || mass ==4000));
    bool ok76X_0T = coupling=="001" || ((coupling=="01" ) && (mass==500 || mass==750 || mass==1000 || mass==1500 || mass == 2000 ||  mass ==3000)) || (( coupling=="02") && (mass==500 ||  mass==1000 || mass==1500 || mass == 750 ||  mass == 3000 ));
    if((whichRel== "80X" && ok80X )||(whichRel== "76X_38T" && ok76X_38T )||(whichRel== "76X_0T" && ok76X_0T ))doPlot=true;
    
    if(doPlot){
      //plotting...
      if(c==0||c==1||c==2||c==3)signalK = (RooDataSet*)w->data(TString::Format("SigWeightK_cat%d",c));
      if(c==4)signalK = (RooDataSet*)w->data("SigWeightK_cat0");
      std::cout<<"flag1"<<std::endl;
      w->Print("v");
      std::cout<<"flag2"<<std::endl;
      RooPlot* plotgg = (w->var("mgg"))->frame(Range(massMin,massMax),Title("mass generated"),Bins(160));
      signalK->Print("v");
      signalK->plotOn(plotgg);
      std::cout<<"flag3"<<std::endl;
      RooPlot* plotg = var->frame(Range(massMin,massMax),Title("mass generated"),Bins(160));
      RooDataSet* fakedata= conv[c]->generate(*var, signalK ->sumEntries());
      fakedata->plotOn(plotg, RooFit::Invisible());
      // reducedmass[c]->plotOn(plotg, LineColor(kRed));
      mgen[c]->plotOn(plotg, LineColor(kGreen));
      conv[c]->plotOn(plotg, LineColor(kBlue));
      // datahist_asimov->plotOn(plotg, MarkerColor(kRed), MarkerSize(0.5));
      //check parameters 
      RooArgSet* model_params = conv[c]->getParameters(*w->var("mgg")) ;
      model_params->Print("v") ;
      std::cout<<"flag4"<<std::endl;
    
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotgg->getObject(0),"Reco Simulation","lp");    
    //    legmc->AddEntry(plotg->getObject(0),"Response ","l");
    legmc->AddEntry(plotg->getObject(1),"Gen shape","l");   
    legmc->AddEntry(plotg->getObject(2),"Convolution","l");   
    
    // bw[c]->plotOn(plotg, LineColor(kBlue));
    plotg->GetYaxis()->SetRangeUser(0.01,30000 );
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    plotg->Draw();
    plotgg->Draw("same");
    legmc->Draw("same");
    c1->SetLogy(0);
    

    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mreco_cat%d_M%d_k"+coupling+".png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mreco_cat%d_M%d_k"+coupling+".pdf",c,iMass)); 
    					        
    c1->SetLogy();			        
    					        
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mreco_cat%d_M%d_k"+coupling+"_log.png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/mreco_cat%d_M%d_k"+coupling+"_log.pdf",c,iMass)); 
    c1->SetLogy(0);  
    }
    
  }
		 
		 
		 
		 
}






void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, TString coupling, TString whichRel) {
  int ncat = NCAT;
  RooDataHist* signal;
  RooVoigtian* voigt[NCAT+1];
  RooDCBShape* cb[NCAT+1];
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelcms(0, "2012", true);
  TPaveText* label_sqrt = get_labelsqrt(0);
  int iMass = (int)abs(mass);
  for(int c = 0; c<ncat+1; c++){
    if(c==2||c==3)continue;
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TString myCut;
    
     //pickup the functions from the ws

    RooRealVar* gH = new RooRealVar("gH", "gH", 0, 1000);
    
    if(c<4)    signal = (RooDataHist*)w->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   signal = (RooDataHist*)w->data (TString::Format("signal_asimov",c));
   
 
    //crystallball
    //cb pos                                  

    RooRealVar* cb_deltaMean = new RooRealVar(TString::Format("cb_deltaMean_cat%d",c),TString::Format("cb_deltaMean_cat%d",c), 0., -300., 300.);
    RooFormulaVar* cb_mean = new RooFormulaVar(TString::Format("cb_mean_cat%d",c),"@0+@1", RooArgList(*cb_deltaMean, *w->var("MH")));
    RooRealVar* cb_sigma = new RooRealVar(TString::Format("cb_sigma_cat%d",c),TString::Format("cb_sigma_cat%d",c),20, 0., 200);//90
    RooRealVar* cb_alphaL = new RooRealVar(TString::Format("cb_alphaL_cat%d",c),TString::Format("cb_alphaL_cat%d",c),1.5, 0.7, 1.7);
    RooRealVar* cb_alphaR = new RooRealVar(TString::Format("cb_alphaR_cat%d",c),TString::Format("cb_alphaR_cat%d",c),2, 1.5, 4);
    //cb_alphaR->setConstant();
    //cb_alphaL->setConstant();
    RooRealVar* cb_nL = new RooRealVar(TString::Format("cb_nL_cat%d",c),TString::Format("cb_nL_cat%d",c), 3.3, 3, 4.5);
    RooRealVar* cb_nR = new RooRealVar(TString::Format("cb_nR_cat%d",c),TString::Format("cb_nR_cat%d",c), 4., 3, 6);
    if(coupling!="001") cb_nR->setRange(0.,15);
    if(coupling!="001") cb_nR->setVal(2);
    if(coupling!="001")cb_nL->setRange(0., 15);
    if(coupling!="001")cb_nL->setVal(2);
    if(coupling!="001") cb_alphaR->setRange(0.,5);
    if(coupling!="001") cb_alphaR->setVal(2);
    if(coupling!="001")cb_alphaL->setRange(0., 5);
    if(coupling!="001")cb_alphaL->setVal(2);
    
    if(mass>3200&&c==0&&coupling=="001") cb_nR->setRange(4.,5);
    if(mass>3200&&c==0&&coupling=="001") cb_nR->setVal(4.3);
    if(mass>3200&&c==0&&coupling=="001")cb_nL->setRange(4.5,6.);
    if(mass>3200&&c==0&&coupling=="001")cb_nL->setVal(5);
    //if(mass>3200&&c==1)cb_nL->setRange(3,4.5);
    //cb_nR->setConstant();
    //cb_nL->setConstant();
    cb[c] =  new RooDCBShape(TString::Format("cbshape_cat%d",c),TString::Format("cbshape_cat%d",c) , *w->var("var"), *cb_mean, *cb_sigma,  *cb_alphaL, *cb_alphaR,  *cb_nL,*cb_nR) ;
    w->import(*cb[c]);
   
    //fitting
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    RooFitResult* fitresults = (RooFitResult* ) cb[c]->fitTo(*signal,SumW2Error(kTRUE), Range(mass*0.7,mass*1.3), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
  
    //plotting...
    double massMin,massMax;
    if(coupling=="001"){
      massMin=0.8*mass;
      massMax=1.2*mass;
    }else if(coupling=="01"){
      massMin=0.6*mass;
      massMax=1.4*mass;
    }else if(coupling=="02"){
      massMin=0.4*mass;
      massMax=1.6*mass;
    }
    if(massMin<300)massMin=300;
    if(massMax>6000)massMax=6000;
    RooPlot* plotg = (w->var("var"))->frame(Range(massMin,massMax),Title("mass"),Bins(60));
    signal->plotOn(plotg);
    
    // voigt[c]->plotOn(plotg, LineColor(kBlue));
    // cb[c]->plotOn(plotg, LineColor(kBlue));
    //check parameters 
    // RooArgSet* model_params = voigt[c]->getParameters(*w->var("mgg")) ;
    //model_params->Print("v") ;
    RooArgList model_params = fitresults->floatParsFinal();
    model_params.Print("v") ;
    //save a new function in the ws as vittorio said

    RooRealVar* dcb_dm = new RooRealVar(TString::Format("dcb_dm_cat%d",c),TString::Format("dcb_dm_cat%d",c), cb_deltaMean->getVal());
    RooRealVar* dcb_m = new RooRealVar(TString::Format("dcb_m_cat%d",c),TString::Format("dcb_m_cat%d",c), cb_mean->getVal());
    RooRealVar* dcb_s = new RooRealVar(TString::Format("dcb_s_cat%d",c),TString::Format("dcb_s_cat%d",c), cb_sigma->getVal());
    RooRealVar* dcb_aL = new RooRealVar(TString::Format("dcb_aL_cat%d",c),TString::Format("dcb_aL_cat%d",c),cb_alphaL->getVal());
    RooRealVar* dcb_aR = new RooRealVar(TString::Format("dcb_aR_cat%d",c),TString::Format("dcb_aR_cat%d",c),cb_alphaR->getVal());
    RooRealVar* dcb_nL = new RooRealVar(TString::Format("dcb_nL_cat%d",c),TString::Format("dcb_nL_cat%d",c),cb_nL->getVal() );
    RooRealVar* dcb_nR = new RooRealVar(TString::Format("dcb_nR_cat%d",c),TString::Format("dcb_nR_cat%d",c),cb_nR->getVal() );
    dcb_dm->setError(cb_deltaMean->getError());
    dcb_m->setError(cb_deltaMean->getError());
    dcb_s->setError(cb_sigma->getError());
    dcb_aL->setError(cb_alphaL->getError());
    dcb_aR->setError(cb_alphaR->getError());
    dcb_nL->setError(cb_nL->getError());
    dcb_nR->setError(cb_nR->getError());
    RooDCBShape*  dcb =  new RooDCBShape(TString::Format("dcbshape_cat%d",c),TString::Format("dcbshape_cat%d",c) , *w->var("var"), *dcb_m, *dcb_s,  *dcb_aL, *dcb_aR,  *dcb_nL,*dcb_nR) ;
    w->import(*dcb);
    w->import(*dcb_dm);
    dcb->plotOn(plotg, LineColor(kBlue));
    RooArgSet* model_params2 = dcb->getParameters(*w->var("var")) ;
    model_params2->Print("v") ;

    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotg->getObject(0),"Asimov Data","lp");    
    //  legmc->AddEntry(plotg->getObject(1),"Voigtian Fit ","l");
    legmc->AddEntry(plotg->getObject(1),"DCB Fit ","l");
    
    // bw[c]->plotOn(plotg, LineColor(kBlue));
    plotg->GetYaxis()->SetRangeUser(0.1,30000 );
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    plotg->Draw();
    legmc->Draw("same");
    c1->SetLogy(0);
    

    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/fitasimov_cat%d_M%d_k"+coupling+".png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/fitasimov_cat%d_M%d_k"+coupling+".pdf",c,iMass)); 
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/fitasimov_cat%d_M%d_k"+coupling+"_log.png",c,iMass)); 
    c1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/fitasimov_cat%d_M%d_k"+coupling+"_log.pdf",c,iMass)); 
    c1->SetLogy(0);  
  }
		 
		 
		 
		 
}

RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov )
{
  RooBinning binning(360, 300, 6000, "binning");
    RooArgSet mset( *x );
    if( asimov != 0 ) {
        asimov->reset();
    } else {
      asimov = new RooDataHist(Form("asimov_dataset_%s",pdf->GetName()),Form("asimov_dataset_%s",pdf->GetName()),mset, "binning");
    }
    pdf->fillDataHist( asimov, &mset, 1, false, true );

    for( int i = 0 ; i < asimov->numEntries() ; i++ ) {
        asimov->get( i ) ;

        // Expected data, multiply p.d.f by nEvents
        Double_t w = asimov->weight() * nexp;
        asimov->set( w, sqrt( w ) );
	std::cout<<i<<" "<<w<<std::endl;
    }
    std::cout<<nexp<<" "<<asimov->sumEntries()<<std::endl;
    Double_t corr = nexp / asimov->sumEntries();
    for( int i = 0 ; i < asimov->numEntries() ; i++ ) {
        RooArgSet theSet = *( asimov->get( i ) );
        asimov->set( asimov->weight()*corr, sqrt( asimov->weight()*corr ) );
    }
    
    return asimov;
}



void plotAllSignalsAsimov(TString coupling, std::string fcn, TString whichRel){


  RooDCBShape final_shape[5];
  std::cout<<"flag1"<<std::endl;
  TFile* f500 = TFile::Open(whichRel+"/ws_ResponseAndGen_M500_k"+coupling+"_final.root");
  TFile* f750 = TFile::Open(whichRel+"/ws_ResponseAndGen_M750_k"+coupling+"_final.root");
  TFile* f1000 = TFile::Open(whichRel+"/ws_ResponseAndGen_M1000_k"+coupling+"_final.root");
  TFile* f1250 = TFile::Open(whichRel+"/ws_ResponseAndGen_M1250_k"+coupling+"_final.root");
  // TFile* f1500 = TFile::Open(whichRel+"/ws_ResponseAndGen_M1500_k"+coupling+"_final.root");
  TFile* f1750 = TFile::Open(whichRel+"/ws_ResponseAndGen_M1750_k"+coupling+"_final.root");
  TFile* f2000 = TFile::Open(whichRel+"/ws_ResponseAndGen_M2000_k"+coupling+"_final.root");
  TFile* f2250 = TFile::Open(whichRel+"/ws_ResponseAndGen_M2250_k"+coupling+"_final.root");
  TFile* f2500 = TFile::Open(whichRel+"/ws_ResponseAndGen_M2500_k"+coupling+"_final.root");
  TFile* f2750 = TFile::Open(whichRel+"/ws_ResponseAndGen_M2750_k"+coupling+"_final.root");
  TFile* f3000 = TFile::Open(whichRel+"/ws_ResponseAndGen_M3000_k"+coupling+"_final.root");
  //TFile* f3250 = TFile::Open(whichRel+"/ws_ResponseAndGen_M3250_k"+coupling+"_final.root");
  //TFile* f3500 = TFile::Open(whichRel+"/ws_ResponseAndGen_M3500_k"+coupling+"_final.root");
  //TFile* f3750 = TFile::Open(whichRel+"/ws_ResponseAndGen_M3750_k"+coupling+"_final.root");
  TFile* f4000 = TFile::Open(whichRel+"/ws_ResponseAndGen_M4000_k"+coupling+"_final.root");

  
  RooWorkspace* ws500 = (RooWorkspace*) f500->Get("HLFactory_ws");
  RooWorkspace* ws750 =  (RooWorkspace*)f750->Get("HLFactory_ws");
  RooWorkspace* ws1000 =  (RooWorkspace*)f1000->Get("HLFactory_ws");
  RooWorkspace* ws1250 =  (RooWorkspace*)f1250->Get("HLFactory_ws");
  //  RooWorkspace* ws1500 =  (RooWorkspace*)f1500->Get("HLFactory_ws");
  RooWorkspace* ws1750 =  (RooWorkspace*)f1750->Get("HLFactory_ws");
  RooWorkspace* ws2000 =  (RooWorkspace*)f2000->Get("HLFactory_ws");
  RooWorkspace* ws2250 =  (RooWorkspace*)f2250->Get("HLFactory_ws");
  RooWorkspace* ws2500 =  (RooWorkspace*)f2500->Get("HLFactory_ws");
  RooWorkspace* ws2750 =  (RooWorkspace*)f2750->Get("HLFactory_ws");
  RooWorkspace* ws3000 =  (RooWorkspace*)f3000->Get("HLFactory_ws");
  //RooWorkspace* ws3250 =  (RooWorkspace*)f3250->Get("HLFactory_ws");
  //RooWorkspace* ws3500 =  (RooWorkspace*)f3500->Get("HLFactory_ws");
  //RooWorkspace* ws3750 =  (RooWorkspace*)f3750->Get("HLFactory_ws");
  RooWorkspace* ws4000 =  (RooWorkspace*)f4000->Get("HLFactory_ws");
  std::cout<<"flag1www"<<std::endl;
  TCanvas* cc1 = new TCanvas("cc1", "cc1",1);
  cc1->cd();
  double m[11];
  double mErr[11];
  double s[11];
  double sErr[11];
  double aL[11];
  double aLErr[11];
  double aR[11];
  double aRErr[11];
  double nL[11];
  double nLErr[11];
  double nR[11];
  double nRErr[11];
  double masses[11] = {500, 750,1000,1150, 1750, 2000,2250,2500,2750,3000,4000};//1500 3250 3500 3750,
  double massesErr[11] = {0., 0., 0., 0,0., 0., 0., 0,0.,  0.,  0.};
  TString svar;
  TString spdf;
  if(fcn=="response"){
    spdf= "responseaddpdf";
    svar= "responseaddpdf";
  }
  
  RooRealVar* MH = new RooRealVar("MH", "MH", 300, 6000);
  MH->setConstant();
  RooWorkspace* ws = new RooWorkspace(whichRel+"/ws_inputs");
   RooDCBShape* sigshape[NCAT+1];
   TF1* fm[NCAT+1];
   TF1* fs[NCAT+1];
   TF1* faL[NCAT+1];
   TF1* faR[NCAT+1];
   TF1* fnL[NCAT+1];
   TF1* fnR[NCAT+1];
   TF1* ffhmw[NCAT+1];
   for(int c =0; c< NCAT+1; c++){
    if(c==2||c==3)continue;
    std::cout<<"flag1"<<std::endl;
    RooDCBShape* res500 = (RooDCBShape*) ws500->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res750 = (RooDCBShape*) ws750->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res1000 = (RooDCBShape*) ws1000->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res1250 = (RooDCBShape*) ws1250->pdf(TString::Format("dcbshape_cat%d",c));
    //   RooDCBShape* res1500 = (RooDCBShape*) ws1500->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res1750 = (RooDCBShape*) ws1750->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res2000 = (RooDCBShape*) ws2000->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res2250 = (RooDCBShape*) ws2250->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res2500 = (RooDCBShape*) ws2500->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res2750 = (RooDCBShape*) ws2750->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res3000 = (RooDCBShape*) ws3000->pdf(TString::Format("dcbshape_cat%d",c));
    //RooDCBShape* res3250 = (RooDCBShape*) ws3250->pdf(TString::Format("dcbshape_cat%d",c));
    //RooDCBShape* res3500 = (RooDCBShape*) ws3500->pdf(TString::Format("dcbshape_cat%d",c));
    // RooDCBShape* res3750 = (RooDCBShape*) ws3750->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res4000 = (RooDCBShape*) ws4000->pdf(TString::Format("dcbshape_cat%d",c));
   

   
    //compute FWHM
    double fhwm[11];
    fhwm[0] = computePdfFHWM(*res500,*ws500->var("var"), 500);
    fhwm[1] = computePdfFHWM(*res750,*ws750->var("var"), 750);
    fhwm[2] = computePdfFHWM(*res1000,*ws1000->var("var"), 1000);
    fhwm[3] = computePdfFHWM(*res1250,*ws1250->var("var"), 1250);
    fhwm[4] = computePdfFHWM(*res1750,*ws1750->var("var"), 1750);
    fhwm[5] = computePdfFHWM(*res2000,*ws2000->var("var"), 2000);
    fhwm[6] = computePdfFHWM(*res2250,*ws2250->var("var"), 2250);
    fhwm[7] = computePdfFHWM(*res2500,*ws2500->var("var"), 2500);
    fhwm[8] = computePdfFHWM(*res2750,*ws2750->var("var"), 2750);
    fhwm[9] = computePdfFHWM(*res3000,*ws3000->var("var"), 3000);
    //fhwm[10] = computePdfFHWM(*res3750,*ws3750->var("var"), 3750);
    fhwm[10] = computePdfFHWM(*res4000,*ws4000->var("var"), 4000);
   
    for(int i = 0; i<11; i++)std::cout<<fhwm[i]<<" "<<masses[i]<<std::endl;
    TCanvas* cf = new TCanvas("cf", "",1);
    cf->cd();
    TGraph* gr = new TGraph(11,masses,fhwm);
    ffhmw[c] = new TF1(TString::Format("ffhmw_cat%d",c), "pol1", 500,4000);
    gr->Fit(TString::Format("ffhmw_cat%d",c), "R");
    gr->GetXaxis()->SetTitle("m_X[GeV]");
    gr->GetYaxis()->SetTitle("FHWM [GeV]");
    gr->Draw("AP");   
    cf->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_001/FHWM.png_cat%d",c));
    


    TH1F* h500=new TH1F("h500","",460, 300,6000);
    h500->SetLineColor(kBlue);
    h500->SetMarkerColor(kBlue);
    TH1F* h750=new TH1F("h750","",460, 300,6000);
    h750->SetLineColor(kMagenta);
    h750->SetMarkerColor(kMagenta);
    TH1F* h1000=new TH1F("h1000","",460, 300,6000);
    h1000->SetLineColor(kRed);
    h1000->SetMarkerColor(kRed);
    TH1F* h1250=new TH1F("h1250","",460, 300,6000);
    h1250->SetLineColor(kViolet);
    h1250->SetMarkerColor(kViolet);
    TH1F* h1500=new TH1F("h1500","",460, 300,6000);
    h1500->SetLineColor(kPink);
    h1500->SetMarkerColor(kPink);
    TH1F* h1750=new TH1F("h1750","",460, 300,6000);
    h1750->SetLineColor(kGreen);
    h1750->SetMarkerColor(kGreen);
    TH1F* h2000=new TH1F("h2000","",460, 300,6000);
    h2000->SetLineColor(kOrange);
    h2000->SetMarkerColor(kOrange);
    TH1F* h2250=new TH1F("h2250","",460, 300,6000);
    h2250->SetLineColor(kBlack);
    h2250->SetMarkerColor(kBlack);
    TH1F* h2500=new TH1F("h2500","",460, 300,6000);
    h2500->SetLineColor(kAzure);
    h2500->SetMarkerColor(kAzure);
    TH1F* h2750=new TH1F("h2750","",460, 300,6000);
    h2750->SetLineColor(kGray);
    h2750->SetMarkerColor(kGray);
    TH1F* h3000=new TH1F("h3000","",460, 300,6000);
    h3000->SetLineColor(kSpring);
    h3000->SetMarkerColor(kSpring);
    TH1F* h3250=new TH1F("h3250","",460, 300,6000);
    h3250->SetLineColor(kBlue+8);
    h3250->SetMarkerColor(kBlue+8);
    TH1F* h3500=new TH1F("h3500","",460, 300,6000);
    h3500->SetLineColor(kRed+9);
    h3500->SetMarkerColor(kRed+9);
    TH1F* h3750=new TH1F("h3750","",460, 300,6000);
    h3750->SetLineColor(kRed+2);
    h3750->SetMarkerColor(kRed+2);
    TH1F* h4000=new TH1F("h4000","",460, 300,6000);
    h4000->SetLineColor(kOrange+4);
    h4000->SetMarkerColor(kOrange+4);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(h500,"M_{X} = 500 GeV","pl");    
    legmc->AddEntry(h750,"M_{X} = 750 GeV","pl");    
    legmc->AddEntry(h1000,"M_{X} = 1000 GeV","pl");    
    legmc->AddEntry(h1250,"M_{X} = 1250 GeV","pl");    
    //  legmc->AddEntry(h1500,"M_{X} = 1500 GeV","pl");    
    legmc->AddEntry(h1750,"M_{X} = 1750 GeV","pl");    
    legmc->AddEntry(h2000,"M_{X} = 2000 GeV","pl");    
    legmc->AddEntry(h2250,"M_{X} = 2250 GeV","pl");    
    legmc->AddEntry(h2500,"M_{X} = 2500 GeV","pl");    
    legmc->AddEntry(h2750,"M_{X} = 2750 GeV","pl");    
    legmc->AddEntry(h3000,"M_{X} = 3000 GeV","pl");    
    //legmc->AddEntry(h3250,"M_{X} = 3250 GeV","pl");    
    //legmc->AddEntry(h3500,"M_{X} = 3500 GeV","pl");    
    //legmc->AddEntry(h3750,"M_{X} = 3750 GeV","pl");    
    legmc->AddEntry(h4000,"M_{X} = 4000 GeV","pl");    
   
    RooDataHist* resdata500;
    RooDataHist* resdata750;
    RooDataHist* resdata1000;
    RooDataHist* resdata1250;
    RooDataHist* resdata1500;
    RooDataHist* resdata1750;
    RooDataHist* resdata2000;
    RooDataHist* resdata2250;
    RooDataHist* resdata2500;
    RooDataHist* resdata2750;
    RooDataHist* resdata3000;
    RooDataHist* resdata3250;
    RooDataHist* resdata3500;
    RooDataHist* resdata3750;
    RooDataHist* resdata4000;

 
    if(c<4)    resdata500 = (RooDataHist*)ws500->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata500 = (RooDataHist*)ws500->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata750 = (RooDataHist*)ws750->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata750 = (RooDataHist*)ws750->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata1000 = (RooDataHist*)ws1000->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata1000 = (RooDataHist*)ws1000->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata1250 = (RooDataHist*)ws1250->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata1250 = (RooDataHist*)ws1250->data (TString::Format("signal_asimov",c));
    // if(c<4)    resdata1500 = (RooDataHist*)ws1500->data(TString::Format("signal_asimov_cat%d",c));
    // if(c==4)   resdata1500 = (RooDataHist*)ws1500->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata1750 = (RooDataHist*)ws1750->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata1750 = (RooDataHist*)ws1750->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata2000 = (RooDataHist*)ws2000->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata2000 = (RooDataHist*)ws2000->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata2250 = (RooDataHist*)ws2250->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata2250 = (RooDataHist*)ws2250->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata2500 = (RooDataHist*)ws2500->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata2500 = (RooDataHist*)ws2500->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata2750 = (RooDataHist*)ws2750->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata2750 = (RooDataHist*)ws2750->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata3000 = (RooDataHist*)ws3000->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata3000 = (RooDataHist*)ws3000->data (TString::Format("signal_asimov",c));
    //if(c<4)    resdata3250 = (RooDataHist*)ws3250->data(TString::Format("signal_asimov_cat%d",c));
    // if(c==4)   resdata3250 = (RooDataHist*)ws3250->data (TString::Format("signal_asimov",c));
    //if(c<4)    resdata3500 = (RooDataHist*)ws3500->data(TString::Format("signal_asimov_cat%d",c));
    //if(c==4)   resdata3500 = (RooDataHist*)ws3500->data (TString::Format("signal_asimov",c));
    //if(c<4)    resdata3750 = (RooDataHist*)ws3750->data(TString::Format("signal_asimov_cat%d",c));
    //if(c==4)   resdata3750 = (RooDataHist*)ws3750->data (TString::Format("signal_asimov",c));
    if(c<4)    resdata4000 = (RooDataHist*)ws4000->data(TString::Format("signal_asimov_cat%d",c));
    if(c==4)   resdata4000 = (RooDataHist*)ws4000->data (TString::Format("signal_asimov",c));
  









    std::cout<<"flag2"<<std::endl;
    RooPlot* pres500 = ws750->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    std::cout<<"flag3"<<std::endl;
    resdata500->plotOn(pres500,MarkerColor(kBlue),LineColor(kBlue));
    std::cout<<"flag3"<<std::endl;
    res500->plotOn(pres500,LineColor(kBlue));
    std::cout<<"flag3"<<std::endl;
    RooPlot* pres750 = ws750->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata750->plotOn(pres750,MarkerColor(kMagenta),LineColor(kMagenta));
    res750->plotOn(pres750,LineColor(kMagenta));
    std::cout<<"flag2"<<std::endl;
    RooPlot* pres1000 = ws1000->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata1000->plotOn(pres1000,MarkerColor(kRed),LineColor(kRed));
    res1000->plotOn(pres1000,LineColor(kRed));
    RooPlot* pres1250 = ws1250->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata1250->plotOn(pres1250,MarkerColor(kViolet),LineColor(kViolet));
    res1250->plotOn(pres1250,LineColor(kViolet));
    // RooPlot* pres1500 = ws1500->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    // resdata1500->plotOn(pres1500,MarkerColor(kPink),LineColor(kPink));
    // res1500->plotOn(pres1500,LineColor(kPink));
    RooPlot* pres1750 = ws1750->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata1750->plotOn(pres1750,MarkerColor(kGreen),LineColor(kGreen));
    res1750->plotOn(pres1750,LineColor(kGreen));
    RooPlot* pres2000 = ws2000->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata2000->plotOn(pres2000,MarkerColor(kOrange),LineColor(kOrange));
    res2000->plotOn(pres2000,LineColor(kOrange));
    std::cout<<"flag2"<<std::endl;
    RooPlot* pres2250 = ws2250->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata2250->plotOn(pres2250,MarkerColor(kBlack),LineColor(kBlack));
    res2250->plotOn(pres2250,LineColor(kBlack));
    RooPlot* pres2500 = ws2500->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata2500->plotOn(pres2500,MarkerColor(kAzure),LineColor(kAzure));
    res2500->plotOn(pres2500,LineColor(kAzure));
    RooPlot* pres2750 = ws2750->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata2750->plotOn(pres2750,MarkerColor(kGray),LineColor(kGray));
    res2750->plotOn(pres2750,LineColor(kGray));
    RooPlot* pres3000 = ws3000->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata3000->plotOn(pres3000,MarkerColor(kSpring),LineColor(kSpring));
    res3000->plotOn(pres3000,LineColor(kSpring));
    std::cout<<"flag2"<<std::endl;
    /*RooPlot* pres3250 = ws3250->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata3250->plotOn(pres3250,MarkerColor(kBlue+8),LineColor(kBlue+8));
    res3250->plotOn(pres3250,LineColor(kBlue+8));
    RooPlot* pres3500 = ws3500->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata3500->plotOn(pres3500,MarkerColor(kRed+9),LineColor(kRed+9));
    res3500->plotOn(pres3500,LineColor(kRed+9));
    RooPlot* pres3750 = ws3750->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata3750->plotOn(pres3750,MarkerColor(kRed+2),LineColor(kRed+2));
    res3750->plotOn(pres3750,LineColor(kRed+2));
    */ RooPlot* pres4000 = ws4000->var("var")->frame(Range(300,5000),Title("mass asimov"),Bins(420));
    resdata4000->plotOn(pres4000,MarkerColor(kOrange+4),LineColor(kOrange+4));
    res4000->plotOn(pres4000,LineColor(kOrange+4));
  
    std::cout<<"flag2"<<std::endl;
    pres500->GetXaxis()->SetTitle("#Delta m [GeV]");
    pres500->GetYaxis()->SetTitle("a.u.");
    pres500->GetYaxis()->SetRangeUser(0.001,10000);
    pres500->Draw();
    pres750->Draw("same");
    pres1000->Draw("same");
    pres1250->Draw("same");
    // pres1500->Draw("same");
    pres1750->Draw("same");
    pres2000->Draw("same");
    pres2250->Draw("same");
    pres2500->Draw("same");
    pres2750->Draw("same");
    pres3000->Draw("same");
    //    pres3250->Draw("same");
    //pres3500->Draw("same");
    //pres3750->Draw("same");
    pres4000->Draw("same");
    std::cout<<"flag2"<<std::endl;
    legmc->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/asimovAllMasses_k"+coupling+"_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/asimovAllMasses_k"+coupling+"_cat%d.pdf",c));
    cc1->SetLogy();
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/asimovAllMasses_k"+coupling+"_cat%d_log.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/asimovAllMasses_k"+coupling+"_cat%d_log.pdf",c));
    cc1->SetLogy(0);

    
    //save parameters 
    RooArgSet* model_params_M500 = res500->getParameters(*ws500->var("var")) ;
    model_params_M500->Print("v") ;
    RooArgSet* model_params_M750 = res750->getParameters(*ws750->var("var")) ;
    model_params_M750->Print("v") ;
    RooArgSet* model_params_M1000 = res1000->getParameters(*ws1000->var("var")) ;
    model_params_M1000->Print("v") ;
    RooArgSet* model_params_M1250 = res1250->getParameters(*ws1250->var("var")) ;
    model_params_M1250->Print("v") ;
    // RooArgSet* model_params_M1500 = res1500->getParameters(*ws1500->var("var")) ;
    //model_params_M1500->Print("v") ;
    RooArgSet* model_params_M1750 = res1750->getParameters(*ws1750->var("var")) ;
    model_params_M1750->Print("v") ;
    RooArgSet* model_params_M2000 = res2000->getParameters(*ws2000->var("var")) ;
    model_params_M2000->Print("v") ;
    RooArgSet* model_params_M2250 = res2250->getParameters(*ws2250->var("var")) ;
    model_params_M2250->Print("v") ;
    RooArgSet* model_params_M2500 = res2500->getParameters(*ws2500->var("var")) ;
    model_params_M2500->Print("v") ;
    RooArgSet* model_params_M2750 = res2750->getParameters(*ws2750->var("var")) ;
    model_params_M2750->Print("v") ;
    RooArgSet* model_params_M3000 = res3000->getParameters(*ws3000->var("var")) ;
    model_params_M3000->Print("v") ;
    /* RooArgSet* model_params_M3250 = res3250->getParameters(*ws3250->var("var")) ;
    model_params_M3250->Print("v") ;
    RooArgSet* model_params_M3500 = res3500->getParameters(*ws3500->var("var")) ;
    model_params_M3500->Print("v") ;
    RooArgSet* model_params_M3750 = res3750->getParameters(*ws3750->var("var")) ;
    model_params_M3750->Print("v") ;*/
    RooArgSet* model_params_M4000 = res4000->getParameters(*ws4000->var("var")) ;
    model_params_M4000->Print("v") ;

    
    m[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws500->var("MH"))->getVal();
    mErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws750->var("MH"))->getVal();
    mErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws1000->var("MH"))->getVal();
    mErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws1250->var("MH"))->getVal();
    mErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    //m[4] = ((RooRealVar*)ws1500->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws1500->var("MH"))->getVal();
    //mErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[4] = ((RooRealVar*)ws1750->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws1750->var("MH"))->getVal();
    mErr[4] =((RooRealVar*)ws1750->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[5] = ((RooRealVar*)ws2000->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws2000->var("MH"))->getVal();
    mErr[5] =((RooRealVar*)ws2000->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[6] = ((RooRealVar*)ws2250->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws2250->var("MH"))->getVal();
    mErr[6] =((RooRealVar*)ws2250->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[7] = ((RooRealVar*)ws2500->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws2500->var("MH"))->getVal();
    mErr[7] =((RooRealVar*)ws2500->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[8] = ((RooRealVar*)ws2750->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws2750->var("MH"))->getVal();
    mErr[8] =((RooRealVar*)ws2750->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[9] = ((RooRealVar*)ws3000->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws3000->var("MH"))->getVal();
    mErr[9] =((RooRealVar*)ws3000->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    /*m[11] = ((RooRealVar*)ws3250->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws3250->var("MH"))->getVal();
    mErr[11] =((RooRealVar*)ws3250->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    m[12] = ((RooRealVar*)ws3500->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws3500->var("MH"))->getVal();
    mErr[12] =((RooRealVar*)ws3500->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
     m[10] = ((RooRealVar*)ws3750->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws3750->var("MH"))->getVal();
    mErr[10] =((RooRealVar*)ws3750->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    */m[10] = ((RooRealVar*)ws4000->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)ws4000->var("MH"))->getVal();
    mErr[10] =((RooRealVar*)ws4000->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
    for(int i=0;i<11;i++)std::cout<<m[i]<<" "<<mErr[i]<<std::endl;
  
    TGraphErrors* gm = new TGraphErrors(11, masses, m, massesErr, mErr);
    fm[c] = new TF1(TString::Format("fm_cat%d",c), "pol2", 500,4000);
    gm->Fit(TString::Format("fm_cat%d",c), "R");
    gm->GetYaxis()->SetTitle("#Delta m= m - m_{H} [GeV]");
    gm->GetXaxis()->SetTitle("m_{X}[GeV]");
    gm->Draw("APE");
    fm[c]->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/Response_"+coupling+"/meanVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/Response_"+coupling+"/meanVsMass_cat%d.pdf",c));
   
    s[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    //s[4] =  (RooRealVar*)ws1500->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    //sErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[4] = ( (RooRealVar*)ws1750->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[4]  =((RooRealVar*)ws1750->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[5] = ( (RooRealVar*)ws2000->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[5]  =((RooRealVar*)ws2000->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[6] = ( (RooRealVar*)ws2250->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[6]  =((RooRealVar*)ws2250->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[7] = ( (RooRealVar*)ws2500->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[7]  =((RooRealVar*)ws2500->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[8] = ( (RooRealVar*)ws2750->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[8]  =((RooRealVar*)ws2750->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    s[9] = ((RooRealVar*)ws3000->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[9]   =((RooRealVar*)ws3000->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    /*  s[10] =   ((RooRealVar*)ws3250->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[10]  =((RooRealVar*)ws3250->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
   /*s[11]  ((RooRealVar*)ws3500->var(TString::Format("dcb_s_cat%d",c)))->getVal();
     sErr[11]  =((RooRealVar*)ws3500->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    
    s[10]   = ((RooRealVar*)ws3750->var(TString::Format("dcb_s_cat%d",c)))->getVal();
    sErr[10]  =((RooRealVar*)ws3750->var(TString::Format("dcb_s_cat%d",c)))->getError() ;*/
   s[10]=  ((RooRealVar*)ws4000->var(TString::Format("dcb_s_cat%d",c)))->getVal();
   sErr[10]  =((RooRealVar*)ws4000->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
    for(int i=0;i<11;i++)std::cout<<s[i]<<" "<<sErr[i]<<std::endl;
   
    TGraphErrors* gs = new TGraphErrors(11, masses, s, massesErr, sErr);
    gs->GetYaxis()->SetTitle("#sigma [GeV]");
    gs->GetXaxis()->SetTitle("m_{X}[GeV]");
    fs[c] = new TF1(TString::Format("fs_cat%d",c), "pol2", 500,4000);
    gs->Fit(TString::Format("fs_cat%d",c), "R");
    gs->Draw("APE");
    fs[c]->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/sigmaVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/sigmaVsMass_cat%d.pdf",c));
    
    aL[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    //aL[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    //aLErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[4] = ((RooRealVar*)ws1750->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[4] =((RooRealVar*)ws1750->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[5] = ((RooRealVar*)ws2000->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[5] =((RooRealVar*)ws2000->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[6] = ((RooRealVar*)ws2250->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[6] =((RooRealVar*)ws2250->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[7] = ((RooRealVar*)ws2500->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[7] =((RooRealVar*)ws2500->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[8] = ((RooRealVar*)ws2750->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[8] =((RooRealVar*)ws2750->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[9] =  ((RooRealVar*)ws3000->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[9] =((RooRealVar*)ws3000->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
   /*aL[10] = ((RooRealVar*)ws3250->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[10] =((RooRealVar*)ws3250->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[11] ((RooRealVar*)ws3500->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[11] =((RooRealVar*)ws3500->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    aL[10]  = ((RooRealVar*)ws3750->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[10] =((RooRealVar*)ws3750->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
   */aL[10]= ((RooRealVar*)ws4000->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
    aLErr[10] =((RooRealVar*)ws4000->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
    
    
    for(int i=0;i<11;i++)std::cout<<aL[i]<<" "<<aLErr[i]<<std::endl;
  
    TGraphErrors* gaL = new TGraphErrors(11, masses, aL, massesErr, aLErr);
    gaL->GetYaxis()->SetTitle("#alpha_{L} [GeV]");
    gaL->GetXaxis()->SetTitle("m_{X}[GeV]");
    faL[c] = new TF1(TString::Format("faL_cat%d",c), "pol2", 500,4000);
    gaL->Fit(TString::Format("faL_cat%d",c), "R");
    gaL->Draw("APE");
    faL[c]->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/aLVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/aLVsMass_cat%d.pdf",c));

      
    aR[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    // aR[4] = ((RooRealVar*)ws1500->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    // aRErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[4] = ((RooRealVar*)ws1750->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[4] =((RooRealVar*)ws1750->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[5] = ((RooRealVar*)ws2000->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[5] =((RooRealVar*)ws2000->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[6] = ((RooRealVar*)ws2250->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[6] =((RooRealVar*)ws2250->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[7] = ((RooRealVar*)ws2500->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[7] =((RooRealVar*)ws2500->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[8] = ((RooRealVar*)ws2750->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[8] =((RooRealVar*)ws2750->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[9] = ((RooRealVar*)ws3000->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[9] =((RooRealVar*)ws3000->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    /*aR[11] = ((RooRealVar*)ws3250->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[11] =((RooRealVar*)ws3250->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[12] = ((RooRealVar*)ws3500->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[12] =((RooRealVar*)ws3500->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    aR[10] = ((RooRealVar*)ws3750->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[10] =((RooRealVar*)ws3750->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    */ aR[10] = ((RooRealVar*)ws4000->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
    aRErr[10] =((RooRealVar*)ws4000->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
    for(int i=0;i<11;i++)std::cout<<aR[i]<<" "<<aRErr[i]<<std::endl;
  
    TGraphErrors* gaR = new TGraphErrors(11, masses, aR, massesErr, aRErr);
    gaR->GetYaxis()->SetTitle("#alpha_{R} [GeV]");
    gaR->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    faR[c] = new TF1(TString::Format("faR_cat%d",c), "pol2", 500,4000);
    gaR->Fit(TString::Format("faR_cat%d",c), "R");
    gaR->Draw("APE");
    faR[c]->Draw("same");
    
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/aRVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/aRVsMass_cat%d.pdf",c));
 
    nR[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    // nR[4] = ((RooRealVar*)ws1500->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    // nRErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[4] = ((RooRealVar*)ws1750->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[4] =((RooRealVar*)ws1750->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[5] = ((RooRealVar*)ws2000->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[5] =((RooRealVar*)ws2000->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[6] = ((RooRealVar*)ws2250->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[6] =((RooRealVar*)ws2250->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[7] = ((RooRealVar*)ws2500->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[7] =((RooRealVar*)ws2500->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[8] = ((RooRealVar*)ws2750->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[8] =((RooRealVar*)ws2750->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[9] = ((RooRealVar*)ws3000->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[9] =((RooRealVar*)ws3000->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    /*nR[11] = ((RooRealVar*)ws3250->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[11] =((RooRealVar*)ws3250->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[12] = ((RooRealVar*)ws3500->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[12] =((RooRealVar*)ws3500->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    nR[10] = ((RooRealVar*)ws3750->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[10] =((RooRealVar*)ws3750->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    */nR[10] = ((RooRealVar*)ws4000->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
    nRErr[10] =((RooRealVar*)ws4000->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
    for(int i=0;i<11;i++)std::cout<<nR[i]<<" "<<nRErr[i]<<std::endl;
  
    TGraphErrors* gnR = new TGraphErrors(11, masses, nR, massesErr, nRErr);
    gnR->GetYaxis()->SetTitle("n_{R} [GeV]");
    gnR->GetXaxis()->SetTitle("m_{X}[GeV]");
    fnR[c] = new TF1(TString::Format("fnR_cat%d",c), "pol2", 500,4000);
    gnR->Fit(TString::Format("fnR_cat%d",c), "R");
    gnR->Draw("APE");
    fnR[c]->Draw("same");   
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/nRVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/nRVsMass_cat%d.pdf",c));

    nL[0] = ((RooRealVar*)ws500->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[0] =((RooRealVar*)ws500->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[1] = ((RooRealVar*)ws750->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[1] =((RooRealVar*)ws750->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[2] = ((RooRealVar*)ws1000->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[2] =((RooRealVar*)ws1000->var(TString::Format("dcb_nL_cat%d",c)))->getError(); 
    nL[3] = ((RooRealVar*)ws1250->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[3] =((RooRealVar*)ws1250->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    // nL[4] = ((RooRealVar*)ws1500->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    //nLErr[4] =((RooRealVar*)ws1500->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[4] = ((RooRealVar*)ws1750->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[4] =((RooRealVar*)ws1750->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[5] = ((RooRealVar*)ws2000->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[5] =((RooRealVar*)ws2000->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[6] = ((RooRealVar*)ws2250->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[6] =((RooRealVar*)ws2250->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[7] = ((RooRealVar*)ws2500->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[7] =((RooRealVar*)ws2500->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[8] = ((RooRealVar*)ws2750->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[8] =((RooRealVar*)ws2750->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    nL[9] = ((RooRealVar*)ws3000->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[9] =((RooRealVar*)ws3000->var(TString::Format("dcb_nL_cat%d",c)))->getError();
    /*nL[11] = ((RooRealVar*)ws3250->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[11] =((RooRealVar*)ws3250->var(TString::Format("dcb_nL_cat%d",c)))->getError();
    nL[12] = ((RooRealVar*)ws3500->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[12] =((RooRealVar*)ws3500->var(TString::Format("dcb_nL_cat%d",c)))->getError();
    nL[10] = ((RooRealVar*)ws3750->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[10] =((RooRealVar*)ws3750->var(TString::Format("dcb_nL_cat%d",c)))->getError();
    */ nL[10] = ((RooRealVar*)ws4000->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
    nLErr[10] =((RooRealVar*)ws4000->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
    for(int i=0;i<11;i++)std::cout<<nL[i]<<" "<<nLErr[i]<<std::endl;
    TGraphErrors* gnL = new TGraphErrors(11, masses, nL, massesErr, nLErr);
    gnL->GetYaxis()->SetTitle("n_{L} [GeV]");
    gnL->GetXaxis()->SetTitle("m_{X} [GeV]");
    fnL[c] = new TF1(TString::Format("fnL_cat%d",c), "pol2", 500,4000);
    gnL->Fit(TString::Format("fnL_cat%d",c), "R");
    gnL->Draw("APE");
    fnL[c]->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/nLVsMass_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/nLVsMass_cat%d.pdf",c));
    


    //build parametric model
    //TF1: p0+p1*x+p2*x*x
    
    RooRealVar* MH = new RooRealVar("MH", "MH", 300, 6000);
    MH->setConstant();
    RooRealVar* p0m = new RooRealVar(TString::Format("p0m_cat%d",c), TString::Format("p0m_cat%d",c),fm[c]->GetParameter(0));
    RooRealVar* p1m = new RooRealVar(TString::Format("p1m_cat%d",c), TString::Format("p1m_cat%d",c),fm[c]->GetParameter(1));
    RooRealVar* p2m = new RooRealVar(TString::Format("p2m_cat%d",c), TString::Format("p2m_cat%d",c),fm[c]->GetParameter(2));
   		    		                     	                    
    RooRealVar* p0s = new RooRealVar(TString::Format("p0s_cat%d",c), TString::Format("p0s_cat%d",c),fs[c]->GetParameter(0));
    RooRealVar* p1s = new RooRealVar(TString::Format("p1s_cat%d",c), TString::Format("p1s_cat%d",c),fs[c]->GetParameter(1));
    RooRealVar* p2s = new RooRealVar(TString::Format("p2s_cat%d",c), TString::Format("p2s_cat%d",c),fs[c]->GetParameter(2));
     		    
    RooRealVar* p0aL = new RooRealVar(TString::Format("p0aL_cat%d",c), TString::Format("p0aL_cat%d",c),faL[c]->GetParameter(0));
    RooRealVar* p1aL = new RooRealVar(TString::Format("p1aL_cat%d",c), TString::Format("p1aL_cat%d",c),faL[c]->GetParameter(1));
    RooRealVar* p2aL = new RooRealVar(TString::Format("p2aL_cat%d",c), TString::Format("p2aL_cat%d",c),faL[c]->GetParameter(2));
     		    					                      
    RooRealVar* p0aR = new RooRealVar(TString::Format("p0aR_cat%d",c), TString::Format("p0aR_cat%d",c),faR[c]->GetParameter(0));
    RooRealVar* p1aR = new RooRealVar(TString::Format("p1aR_cat%d",c), TString::Format("p1aR_cat%d",c),faR[c]->GetParameter(1));
    RooRealVar* p2aR = new RooRealVar(TString::Format("p2aR_cat%d",c), TString::Format("p2aR_cat%d",c),faR[c]->GetParameter(2));  
    		    
    RooRealVar* p0nL = new RooRealVar(TString::Format("p0nL_cat%d",c), TString::Format("p0nL_cat%d",c),fnL[c]->GetParameter(0));
    RooRealVar* p1nL = new RooRealVar(TString::Format("p1nL_cat%d",c), TString::Format("p1nL_cat%d",c),fnL[c]->GetParameter(1));
    RooRealVar* p2nL = new RooRealVar(TString::Format("p2nL_cat%d",c), TString::Format("p2nL_cat%d",c),fnL[c]->GetParameter(2));
   		    					                      
    RooRealVar* p0nR = new RooRealVar(TString::Format("p0nR_cat%d",c), TString::Format("p0nR_cat%d",c),fnR[c]->GetParameter(0));
    RooRealVar* p1nR = new RooRealVar(TString::Format("p1nR_cat%d",c), TString::Format("p1nR_cat%d",c),fnR[c]->GetParameter(1));
    RooRealVar* p2nR = new RooRealVar(TString::Format("p2nR_cat%d",c), TString::Format("p2nR_cat%d",c),fnR[c]->GetParameter(2));
    p0m ->setConstant();
    p1m ->setConstant();
    p2m ->setConstant();
    
    p0s ->setConstant();
    p1s ->setConstant();
    p2s ->setConstant();
    
    p0aL->setConstant();
    p1aL->setConstant();
    p2aL->setConstant();
    
    p0aR->setConstant();
    p1aR->setConstant();
    p2aR->setConstant();
    
    p0nL->setConstant();
    p1nL->setConstant();
    p2nL->setConstant();
    
    p0nR->setConstant();
    p1nR->setConstant();
    p2nR->setConstant();
   
  
    std::cout<<"m: "<<fm[0]->GetParameter(0)<<" "<<fm[0]->GetParameter(1)<<" "<<fm[0]->GetParameter(2)<<" "<<fm[0]->GetParameter(0)+fm[0]->GetParameter(1)*1750+fm[0]->GetParameter(2)*1750*1750+1750<<std::endl;
    std::cout<<"s: "<<fs[0]->GetParameter(0)<<" "<<fs[0]->GetParameter(1)<<" "<<fs[0]->GetParameter(2)<<" "<<fs[0]->GetParameter(0)+fs[0]->GetParameter(1)*1750+fs[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"aL: "<<faL[0]->GetParameter(0)<<" "<<faL[0]->GetParameter(1)<<" "<<faL[0]->GetParameter(2)<<" "<<faL[0]->GetParameter(0)+faL[0]->GetParameter(1)*1750+faL[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"aR: "<<faR[0]->GetParameter(0)<<" "<<faR[0]->GetParameter(1)<<" "<<faR[0]->GetParameter(2)<<faR[0]->GetParameter(0)+faR[0]->GetParameter(1)*1750+faR[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"nL: "<<fnL[0]->GetParameter(0)<<" "<<fnL[0]->GetParameter(1)<<" "<<fnL[0]->GetParameter(2)<<fnL[0]->GetParameter(0)+fnL[0]->GetParameter(1)*1750+fnL[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"nR: "<<fnR[0]->GetParameter(0)<<" "<<fnR[0]->GetParameter(1)<<" "<<fnR[0]->GetParameter(2)<<fnR[0]->GetParameter(0)+fnR[0]->GetParameter(1)*1750+fnR[0]->GetParameter(2)*1750*1750<<std::endl;
  
    RooRealVar* mgg = new RooRealVar("mgg", "mgg", 300, 6000);
    RooRealVar* thetaSmearEBEB = new RooRealVar("thetaSmearEBEB", "thetaSmearEBEB", 0., -1, 1);
    RooRealVar* thetaSmearEBEE = new RooRealVar("thetaSmearEBEE", "thetaSmearEBEE", 0., -1, 1);
    RooRealVar* thetaSmearAll = new RooRealVar("thetaSmearAll", "thetaSmearAll", 0., -1, 1);
    RooRealVar* thetaScaleEBEB = new RooRealVar("thetaSmearEBEB", "thetaSmearEBEB", 0., -1, 1);
    RooRealVar* thetaScaleEBEE = new RooRealVar("thetaSmearEBEE", "thetaSmearEBEE", 0., -1, 1);
    RooRealVar* thetaScaleAll = new RooRealVar("thetaSmearAll", "thetaSmearAll", 0., -1, 1);
    RooRealVar* deltaSmear = new RooRealVar("deltaSmear", "deltaSmear",  0.001);
    RooRealVar* deltaScale = new RooRealVar("deltaSmear", "deltaSmear",  0.001);
    deltaSmear->setConstant();
    RooRealVar* s0_EBEB = new RooRealVar("smear0EBEB", "smear0EBEB",  0.01);
    RooRealVar* s0_EBEE = new RooRealVar("smear0EBEE", "smear0EBEE",  0.015);
    RooRealVar* s0_All = new RooRealVar("smear0All", "smear0All",  0.01);
    s0_EBEB->setConstant();
    s0_EBEE->setConstant();
    std::cout<<"flaggggg"<<std::endl;
    RooFormulaVar* DeltaSmearEBEB = new RooFormulaVar("DeltaSmearEBEB", "2*@0*@1*@2*@2",  RooArgList(*s0_EBEB,*deltaSmear,*MH));
    RooFormulaVar* DeltaSmearEBEE = new RooFormulaVar("DeltaSmearEBEE", "2*@0*@1*@2*@2",  RooArgList(*s0_EBEE,*deltaSmear,*MH));
    RooFormulaVar* DeltaSmearAll = new RooFormulaVar("DeltaSmearAll", "2*@0*@1*@2*@2",  RooArgList(*s0_All,*deltaSmear,*MH));
    RooPlot* plot = (RooPlot*)mgg->frame(Range(300, 6000));
    RooDCBShape* fin_shape[11];
    RooFormulaVar* mean; 
    RooFormulaVar* sigma0;
    RooFormulaVar* sigma; 
    RooFormulaVar* aL; 
    RooFormulaVar* aR; 
    RooFormulaVar* nL; 
    RooFormulaVar* nR; 
    std::cout<<"flaggggg"<<std::endl;
    if(c==0){
      mean= new RooFormulaVar("mean_EBEB", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_EBEB", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_EBEB", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearEBEB,*DeltaSmearEBEB));
      aL= new RooFormulaVar("aL_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_EBEB", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));
    }else if(c==1){
      mean= new RooFormulaVar("mean_EBEE", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_EBEE", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_EBEE", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearEBEE,*DeltaSmearEBEE));
      aL= new RooFormulaVar("aL_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_EBEE", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));

    }else if(c==4){
      mean= new RooFormulaVar("mean_All", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_All", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_All", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearAll,*DeltaSmearAll));
      aL= new RooFormulaVar("aL_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_All", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));

    }
    std::cout<<"flaggggg"<<std::endl;
    for(int M =0; M<11; M++){
      
      MH->setVal(masses[M]);
      // MH->setConstant();
      std::cout<<MH->getVal()<<std::endl;
      
      std::cout<<mean->getVal(*MH)<<std::endl;
      std::cout<<sigma->getVal(*MH)<<std::endl;
      std::cout<<aL->getVal(*MH)<<std::endl;
      std::cout<<aR->getVal(*MH)<<std::endl;
      std::cout<<nL->getVal(*MH)<<std::endl;
      std::cout<<nR->getVal(*MH)<<std::endl;
   
      fin_shape[M] =  new RooDCBShape(TString::Format("SignalShape_kMpl"+coupling+"_cat%d_M%f",c, masses[M]),TString::Format("final_shape_k001_cat%d_M%f",c, masses[M] ) ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
     
     if(masses[M]==500)   fin_shape[M]->plotOn(plot, LineColor(kBlue));
    if(masses[M]==750)   fin_shape[M]->plotOn(plot, LineColor(kMagenta));
    if(masses[M]==1000)   fin_shape[M]->plotOn(plot, LineColor(kRed));
    if(masses[M]==1250)   fin_shape[M]->plotOn(plot, LineColor(kViolet));
    if(masses[M]==1500)   fin_shape[M]->plotOn(plot, LineColor(kPink));
    if(masses[M]==1750)   fin_shape[M]->plotOn(plot, LineColor(kGreen));
    if(masses[M]==2000)   fin_shape[M]->plotOn(plot, LineColor(kOrange));
    if(masses[M]==2250)   fin_shape[M]->plotOn(plot, LineColor(kBlack));
    if(masses[M]==2500)   fin_shape[M]->plotOn(plot, LineColor(kAzure));
    if(masses[M]==2750)   fin_shape[M]->plotOn(plot, LineColor(kGray));
    if(masses[M]==3000)   fin_shape[M]->plotOn(plot, LineColor(kSpring));
    if(masses[M]==3250)   fin_shape[M]->plotOn(plot, LineColor(kBlue+8));
    if(masses[M]==3500)   fin_shape[M]->plotOn(plot, LineColor(kRed+9));
    if(masses[M]==3750)   fin_shape[M]->plotOn(plot, LineColor(kRed+2));
    if(masses[M]==4000)   fin_shape[M]->plotOn(plot, LineColor(kOrange+4));
    }
    plot->GetYaxis()->SetTitle("a.u.");
    plot->GetYaxis()->SetRangeUser(0.001, 5);
    plot->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plot->Draw();
    legmc->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/final_shapes_k001_cat%d.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/final_shapes_k001_cat%d.pdf",c));
    cc1->SetLogy(); 
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/final_shapes_k001_cat%d_log.png",c));
    cc1->SaveAs(TString::Format("~/www/Pippone/"+whichRel+"/Response_"+coupling+"/final_shapes_k001_cat%d_log.pdf",c));
    
    
    if(c==0)    sigshape[c] =  new RooDCBShape("SignalShape_kMpl"+coupling+"_EBEB","SignalShape_kMpl"+coupling+"_EBEB" ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==1)    sigshape[c] =  new RooDCBShape("SignalShape_kMpl"+coupling+"_EBEE","SignalShape_kMpl"+coupling+"_EBEE" ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==4)    sigshape[c] =  new RooDCBShape("SignalShape_kMpl"+coupling+"_All","SignalShape_kMpl"+coupling+"_All" ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    ws->import(*sigshape[c]);
    
    }
   
   RooRealVar* p0_cat0 = new RooRealVar("p0_cat0", "p0_cat0", ffhmw[0]->GetParameter(0) );
   RooRealVar* p1_cat0 = new RooRealVar("p1_cat0", "p1_cat0", ffhmw[0]->GetParameter(1) );
   RooRealVar* p0_cat1 = new RooRealVar("p0_cat1", "p0_cat1", ffhmw[1]->GetParameter(0) );
   RooRealVar* p1_cat1 = new RooRealVar("p1_cat1", "p1_cat1", ffhmw[1]->GetParameter(1) );

   RooFormulaVar* FHWM_EBEB= new RooFormulaVar("FHWM_EBEB", "(@0+@1*@2)", RooArgList(*p0_cat0,*p1_cat0, *MH));
   RooFormulaVar* FHWM_EBEE= new RooFormulaVar("FHWM_EBEE", "(@0+@1*@2)", RooArgList(*p0_cat1,*p1_cat1, *MH));
   ws->import(*FHWM_EBEB);
   ws->import(*FHWM_EBEE);
  
   TFile* fout= new TFile(whichRel+"/SignalParametericShapes80X_ws_kMpl"+coupling+".root", "RECREATE");
   fout->cd();
   fm[0]->Write();       fm[1]->Write(); 
   fs[0]->Write(); 	 fs[1]->Write(); 
   faL[0]->Write();	 faL[1]->Write();
   faR[0]->Write();	 faR[1]->Write();
   fnL[0]->Write();	 fnL[1]->Write();
   fnR[0]->Write();	 fnR[1]->Write();
   ffhmw[0]->Write();	 ffhmw[1]->Write();

   ws->Write();
   fout->Write();
   fout->Close();
}
   
   
   
   
   void testSystematics(TString coupling, Float_t mass){
     TFile* fin= TFile::Open("SignalParametericShapes80X_ws_kMpl"+coupling+".root");
  RooWorkspace* ws =( RooWorkspace*) fin->Get("SignalParametricShapes_kMpl"+coupling);
  RooRealVar* mgg = ws->var("mgg");
  RooRealVar* MH = ws->var("MH");
  RooRealVar* thetaSM = ws->var("thetaSM");
  MH->setVal(mass);
  std::cout<<MH->getVal()<<std::endl;
  MH->setConstant();
  
  for(int c = 0;c<2; c++){
    RooPlot* p = mgg->frame(RooFit::Range(mass*0.8, mass*1.2));  
    thetaSM->setVal(0);
    if(c==0)  {
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEB",c))->plotOn(p, LineColor(kBlack));
      std::cout<<ws->var("sigma0_EBEB")->getVal()<<std::endl;
      std::cout<<ws->var("sigma_EBEB")->getVal()<<std::endl;
      thetaSM->setVal(1);
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEB",c))->plotOn(p, LineColor(kRed));
      std::cout<<ws->var("sigma0_EBEB")->getVal()<<std::endl;
      std::cout<<ws->var("sigma_EBEB")->getVal()<<std::endl;
      thetaSM->setVal(-1);
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEB",c))->plotOn(p, LineColor(kBlue));
      std::cout<<ws->var("sigma0_EBEB")->getVal()<<std::endl;
      std::cout<<ws->var("sigma_EBEB")->getVal()<<std::endl;
     
    }else if(c==1)  {
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEE",c))->plotOn(p, LineColor(kBlack));
      thetaSM->setVal(1);
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEE",c))->plotOn(p, LineColor(kRed));
      thetaSM->setVal(-1);
      ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEE",c))->plotOn(p, LineColor(kBlue));
    }
    p->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    p->GetYaxis()->SetTitle("a.u.");
    TCanvas* cc1 = new TCanvas("cc1", "cc1",1);
    cc1->cd();
    p->Draw();
    cc1->SaveAs(TString::Format("~/www/Pippone/Response_"+coupling+"/smearing_syst_M%.0f_cat%d.png",mass,c));
  }
}


void testShapes(TString coupling, TString whichRel){
  TFile* fin= TFile::Open("SignalParametericShapes80X_ws_kMpl"+coupling+".root");
  RooWorkspace* ws =( RooWorkspace*) fin->Get("SignalParametricShapes_kMpl"+coupling);

  RooRealVar* mgg = ws->var("mgg");
  RooRealVar* MH = ws->var("MH");
  double masses[5]={500, 750,  2000, 3000,4000};
  for(int c=0;c<2;c++){
    RooPlot* p[5]  ;
    for(int m = 0;m<5;m++){

      MH->setVal(masses[m]);  
      
      int iMass = abs(masses[m]);
      TFile* ftrue = TFile::Open(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+".root", iMass));
      RooWorkspace* wtrue = (RooWorkspace*)ftrue->Get("HLFactory_ws");
      RooDataSet* data = (RooDataSet*)wtrue->data(TString::Format("SigWeightK_cat%d",c));
      data->Print("v");
      p[m] = (RooPlot*) mgg->frame(RooFit::Range(300, 5000));
      data->plotOn(p[m], MarkerColor(1+m));
      if(c==0)ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEB",c))->plotOn(p[m], LineColor(m+1));
      if(c==1)ws->pdf(TString::Format("SignalShape_kMpl"+coupling+"_EBEE",c))->plotOn(p[m], LineColor(m+1));
  }
    
    p[0]->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    p[0]->GetYaxis()->SetTitle("a.u.");
    TCanvas* cc1 = new TCanvas("cc1", "cc1",1);
    cc1->cd();
    p[0]->Draw();
    for(int m = 1;m<5;m++)p[m]->Draw("same");
    cc1->SaveAs(TString::Format("~/www/Pippone/Response_"+coupling+"/MassesCheck_cat%d.png",c));
    cc1->SetLogy(); 
    cc1->SaveAs(TString::Format("~/www/Pippone/Response_"+coupling+"/MassesCheck_cat%d.png",c));
  }
}



void plotAllSignalsResponse(TString coupling, std::string fcn){
  std::cout<<"flag1"<<std::endl;
  TFile* f500 = TFile::Open("ws_ResponseAndGen_M500_k"+coupling+".root");
  TFile* f750 = TFile::Open("ws_ResponseAndGen_M750_k"+coupling+".root");
  TFile* f1000 = TFile::Open("ws_ResponseAndGen_M1000_k"+coupling+".root");
  TFile* f1250 = TFile::Open("ws_ResponseAndGen_M1250_k"+coupling+".root");
  TFile* f1500 = TFile::Open("ws_ResponseAndGen_M1500_k"+coupling+".root");
  TFile* f1750 = TFile::Open("ws_ResponseAndGen_M1750_k"+coupling+".root");
  TFile* f2000 = TFile::Open("ws_ResponseAndGen_M2000_k"+coupling+".root");
  TFile* f2250 = TFile::Open("ws_ResponseAndGen_M2250_k"+coupling+".root");
  TFile* f2500 = TFile::Open("ws_ResponseAndGen_M2500_k"+coupling+".root");
  TFile* f2750 = TFile::Open("ws_ResponseAndGen_M2750_k"+coupling+".root");
  TFile* f3000 = TFile::Open("ws_ResponseAndGen_M3000_k"+coupling+".root");
  //TFile* f3250 = TFile::Open("ws_ResponseAndGen_M3250_k"+coupling+".root");
  //TFile* f3500 = TFile::Open("ws_ResponseAndGen_M3500_k"+coupling+".root");
  TFile* f3750 = TFile::Open("ws_ResponseAndGen_M3750_k"+coupling+".root");
  TFile* f4000 = TFile::Open("ws_ResponseAndGen_M4000_k"+coupling+".root");

  std::cout<<"flag1"<<std::endl;
  RooWorkspace* ws500 = (RooWorkspace*) f500->Get("HLFactory_ws");
  RooWorkspace* ws750 =  (RooWorkspace*)f750->Get("HLFactory_ws");
  RooWorkspace* ws1000 =  (RooWorkspace*)f1000->Get("HLFactory_ws");
  RooWorkspace* ws1250 =  (RooWorkspace*)f1250->Get("HLFactory_ws");
  RooWorkspace* ws1500 =  (RooWorkspace*)f1500->Get("HLFactory_ws");
  RooWorkspace* ws1750 =  (RooWorkspace*)f1750->Get("HLFactory_ws");
  RooWorkspace* ws2000 =  (RooWorkspace*)f2000->Get("HLFactory_ws");
  RooWorkspace* ws2250 =  (RooWorkspace*)f2250->Get("HLFactory_ws");
  RooWorkspace* ws2500 =  (RooWorkspace*)f2500->Get("HLFactory_ws");
  RooWorkspace* ws2750 =  (RooWorkspace*)f2750->Get("HLFactory_ws");
  RooWorkspace* ws3000 =  (RooWorkspace*)f3000->Get("HLFactory_ws");
  //RooWorkspace* ws3250 =  (RooWorkspace*)f3250->Get("HLFactory_ws");
  //RooWorkspace* ws3500 =  (RooWorkspace*)f3500->Get("HLFactory_ws");
  RooWorkspace* ws3750 =  (RooWorkspace*)f3750->Get("HLFactory_ws");
  RooWorkspace* ws4000 =  (RooWorkspace*)f4000->Get("HLFactory_ws");

  TCanvas* c1 = new TCanvas("c1", "c1",1);
  c1->cd();

  TString svar;
  TString spdf;
  if(fcn=="response"){
    spdf= "responseaddpdf";
    svar= "responseaddpdf";
  }
   for(int c =0; c< NCAT+1; c++){
    if(c==2||c==3)continue;
    std::cout<<"flag1"<<std::endl;
    RooDCBShape* res500 = (RooDCBShape*) ws500->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res750 = (RooDCBShape*) ws750->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res1000 = (RooDCBShape*) ws1000->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res1250 = (RooDCBShape*) ws1250->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res1500 = (RooDCBShape*) ws1500->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res1750 = (RooDCBShape*) ws1750->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res2000 = (RooDCBShape*) ws2000->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res2250 = (RooDCBShape*) ws2250->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res2500 = (RooDCBShape*) ws2500->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res2750 = (RooDCBShape*) ws2750->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res3000 = (RooDCBShape*) ws3000->pdf(TString::Format("responseaddpdf_cat%d",c));
    //RooDCBShape* res3250 = (RooDCBShape*) ws3250->pdf(TString::Format("responseaddpdf_cat%d",c));
    //RooDCBShape* res3500 = (RooDCBShape*) ws3500->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res3750 = (RooDCBShape*) ws3750->pdf(TString::Format("responseaddpdf_cat%d",c));
    RooDCBShape* res4000 = (RooDCBShape*) ws4000->pdf(TString::Format("responseaddpdf_cat%d",c));
   
    TH1F* h500=new TH1F("h500","",260, -600,600);
    h500->SetLineColor(kBlue);
    h500->SetMarkerColor(kBlue);
    TH1F* h750=new TH1F("h750","",260, -600,600);
    h750->SetLineColor(kMagenta);
    h750->SetMarkerColor(kMagenta);
    TH1F* h1000=new TH1F("h1000","",260, -600,600);
    h1000->SetLineColor(kRed);
    h1000->SetMarkerColor(kRed);
    TH1F* h1250=new TH1F("h1250","",260, -600,600);
    h1250->SetLineColor(kViolet);
    h1250->SetMarkerColor(kViolet);
    TH1F* h1500=new TH1F("h1500","",260, -600,600);
    h1500->SetLineColor(kPink);
    h1500->SetMarkerColor(kPink);
    TH1F* h1750=new TH1F("h1750","",260, -600,600);
    h1750->SetLineColor(kGreen);
    h1750->SetMarkerColor(kGreen);
    TH1F* h2000=new TH1F("h2000","",260, -600,600);
    h2000->SetLineColor(kOrange);
    h2000->SetMarkerColor(kOrange);
    TH1F* h2250=new TH1F("h2250","",260, -600,600);
    h2250->SetLineColor(kBlack);
    h2250->SetMarkerColor(kBlack);
    TH1F* h2500=new TH1F("h2500","",260, -600,600);
    h2500->SetLineColor(kAzure);
    h2500->SetMarkerColor(kAzure);
    TH1F* h2750=new TH1F("h2750","",260, -600,600);
    h2750->SetLineColor(kGray);
    h2750->SetMarkerColor(kGray);
    TH1F* h3000=new TH1F("h3000","",260, -600,600);
    h3000->SetLineColor(kSpring);
    h3000->SetMarkerColor(kSpring);
    TH1F* h3250=new TH1F("h3250","",260, -600,600);
    h3250->SetLineColor(kBlue+8);
    h3250->SetMarkerColor(kBlue+8);
    TH1F* h3500=new TH1F("h3500","",260, -600,600);
    h3500->SetLineColor(kRed+9);
    h3500->SetMarkerColor(kRed+9);
    TH1F* h3750=new TH1F("h3750","",260, -600,600);
    h3750->SetLineColor(kRed+2);
    h3750->SetMarkerColor(kRed+2);
    TH1F* h4000=new TH1F("h4000","",260, -600,600);
    h4000->SetLineColor(kOrange+4);
    h4000->SetMarkerColor(kOrange+4);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(h500,"M_{X} = 500 GeV","pl");    
    legmc->AddEntry(h750,"M_{X} = 750 GeV","pl");    
    legmc->AddEntry(h1000,"M_{X} = 1000 GeV","pl");    
    legmc->AddEntry(h1250,"M_{X} = 1250 GeV","pl");    
    legmc->AddEntry(h1500,"M_{X} = 1500 GeV","pl");    
    legmc->AddEntry(h1750,"M_{X} = 1750 GeV","pl");    
    legmc->AddEntry(h2000,"M_{X} = 2000 GeV","pl");    
    legmc->AddEntry(h2250,"M_{X} = 2250 GeV","pl");    
    legmc->AddEntry(h2500,"M_{X} = 2500 GeV","pl");    
    legmc->AddEntry(h2750,"M_{X} = 2750 GeV","pl");    
    legmc->AddEntry(h3000,"M_{X} = 3000 GeV","pl");    
    //legmc->AddEntry(h3250,"M_{X} = 3250 GeV","pl");    
    //legmc->AddEntry(h3500,"M_{X} = 3500 GeV","pl");    
    legmc->AddEntry(h3750,"M_{X} = 3750 GeV","pl");    
    legmc->AddEntry(h4000,"M_{X} = 4000 GeV","pl");    
   
    RooDataSet* resdata500;
    RooDataSet* resdata750;
    RooDataSet* resdata1000;
    RooDataSet* resdata1250;
    RooDataSet* resdata1500;
    RooDataSet* resdata1750;
    RooDataSet* resdata2000;
    RooDataSet* resdata2250;
    RooDataSet* resdata2500;
    RooDataSet* resdata2750;
    RooDataSet* resdata3000;
    RooDataSet* resdata3250;
    RooDataSet* resdata3500;
    RooDataSet* resdata3750;
    RooDataSet* resdata4000;

    TString myCut;
    
    if(c==0||c==1||c==2||c==3)resdata500 = (RooDataSet*)ws500->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata500 = (RooDataSet*)ws500->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata750 = (RooDataSet*)ws750->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata750 = (RooDataSet*)ws750->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata1000 = (RooDataSet*)ws1000->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata1000 = (RooDataSet*)ws1000->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata1250 = (RooDataSet*)ws1250->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata1250 = (RooDataSet*)ws1250->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata1500 = (RooDataSet*)ws1500->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata1500 = (RooDataSet*)ws1500->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata1750 = (RooDataSet*)ws1750->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata1750 = (RooDataSet*)ws1750->data("SigWeightReduced");    
    if(c==0||c==1||c==2||c==3)resdata2000 = (RooDataSet*)ws2000->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata2000 = (RooDataSet*)ws2000->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata2250 = (RooDataSet*)ws2250->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata2250 = (RooDataSet*)ws2250->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata2500 = (RooDataSet*)ws2500->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata2500 = (RooDataSet*)ws2500->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata2750 = (RooDataSet*)ws2750->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata2750 = (RooDataSet*)ws2750->data("SigWeightReduced");    
    if(c==0||c==1||c==2||c==3)resdata3000 = (RooDataSet*)ws3000->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata3000 = (RooDataSet*)ws3000->data("SigWeightReduced");
    /*if(c==0||c==1||c==2||c==3)resdata3250 = (RooDataSet*)ws3250->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata3250 = (RooDataSet*)ws3250->data("SigWeightReduced");
    if(c==0||c==1||c==2||c==3)resdata3500 = (RooDataSet*)ws3500->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata3500 = (RooDataSet*)ws3500->data("SigWeightReduced");
    */if(c==0||c==1||c==2||c==3)resdata3750 = (RooDataSet*)ws3750->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata3750 = (RooDataSet*)ws3750->data("SigWeightReduced");    
    if(c==0||c==1||c==2||c==3)resdata4000 = (RooDataSet*)ws4000->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==4)resdata4000 = (RooDataSet*)ws4000->data("SigWeightReduced");
    std::cout<<"flag2"<<std::endl;
    RooPlot* pres500 = ws500->var("massReduced")->frame(Range(-600,600),Title("mass reduced"),Bins(120));
    resdata500->plotOn(pres500,MarkerColor(kBlue),LineColor(kBlue));
    res500->plotOn(pres500,LineColor(kBlue));
    RooPlot* pres750 = ws750->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata750->plotOn(pres750,MarkerColor(kMagenta),LineColor(kMagenta));
    res750->plotOn(pres750,LineColor(kMagenta));
    std::cout<<"flag2"<<std::endl;
    RooPlot* pres1000 = ws1000->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata1000->plotOn(pres1000,MarkerColor(kRed),LineColor(kRed));
    res1000->plotOn(pres1000,LineColor(kRed));
    RooPlot* pres1250 = ws1250->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata1250->plotOn(pres1250,MarkerColor(kViolet),LineColor(kViolet));
    res1250->plotOn(pres1250,LineColor(kViolet));
    RooPlot* pres1500 = ws1500->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata1500->plotOn(pres1500,MarkerColor(kPink),LineColor(kPink));
    res1500->plotOn(pres1500,LineColor(kPink));
    RooPlot* pres1750 = ws1750->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata1750->plotOn(pres1750,MarkerColor(kGreen),LineColor(kGreen));
    res1750->plotOn(pres1750,LineColor(kGreen));
    RooPlot* pres2000 = ws2000->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata2000->plotOn(pres2000,MarkerColor(kOrange),LineColor(kOrange));
    res2000->plotOn(pres2000,LineColor(kOrange));
    std::cout<<"flag2"<<std::endl;
    RooPlot* pres2250 = ws2250->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata2250->plotOn(pres2250,MarkerColor(kBlack),LineColor(kBlack));
    res2250->plotOn(pres2250,LineColor(kBlack));
    RooPlot* pres2500 = ws2500->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata2500->plotOn(pres2500,MarkerColor(kAzure),LineColor(kAzure));
    res2500->plotOn(pres2500,LineColor(kAzure));
    RooPlot* pres2750 = ws2750->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata2750->plotOn(pres2750,MarkerColor(kGray),LineColor(kGray));
    res2750->plotOn(pres2750,LineColor(kGray));
    RooPlot* pres3000 = ws3000->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata3000->plotOn(pres3000,MarkerColor(kSpring),LineColor(kSpring));
    res3000->plotOn(pres3000,LineColor(kSpring));
    std::cout<<"flag2"<<std::endl;
    /*RooPlot* pres3250 = ws3250->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata3250->plotOn(pres3250,MarkerColor(kBlue+8),LineColor(kBlue+8));
    res3250->plotOn(pres3250,LineColor(kBlue+8));
    RooPlot* pres3500 = ws3500->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata3500->plotOn(pres3500,MarkerColor(kRed+9),LineColor(kRed+9));
    res3500->plotOn(pres3500,LineColor(kRed+9));
    */RooPlot* pres3750 = ws3750->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata3750->plotOn(pres3750,MarkerColor(kRed+2),LineColor(kRed+2));
    res3750->plotOn(pres3750,LineColor(kRed+2));
    RooPlot* pres4000 = ws4000->var("massReduced")->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    resdata4000->plotOn(pres4000,MarkerColor(kOrange+4),LineColor(kOrange+4));
    res4000->plotOn(pres4000,LineColor(kOrange+4));
  
    std::cout<<"flag2"<<std::endl;
    pres500->GetXaxis()->SetTitle("#Delta m [GeV]");
    pres500->GetYaxis()->SetTitle("a.u.");
    pres500->Draw();
    pres750->Draw("same");
    pres1000->Draw("same");
    pres1250->Draw("same");
    pres1500->Draw("same");
    pres1750->Draw("same");
    pres2000->Draw("same");
    pres2250->Draw("same");
    pres2500->Draw("same");
    pres2750->Draw("same");
    pres3000->Draw("same");
    //    pres3250->Draw("same");
    //pres3500->Draw("same");
    pres3750->Draw("same");
    pres4000->Draw("same");
    std::cout<<"flag2"<<std::endl;
    legmc->Draw("same");
    c1->SaveAs(TString::Format("~/www/Pippone/Response/responsesAllMasses_k"+coupling+"_cat%d.png",c));
    c1->SaveAs(TString::Format("~/www/Pippone/Response/responsesAllMasses_k"+coupling+"_cat%d.pdf",c));
    c1->SetLogy();
    c1->SaveAs(TString::Format("~/www/Pippone/Response/responsesAllMasses_k"+coupling+"_cat%d_log.png",c));
    c1->SaveAs(TString::Format("~/www/Pippone/Response/responsesAllMasses_k"+coupling+"_cat%d_log.pdf",c));
    c1->SetLogy(0);

  
}
}





   
double computePdfFHWM(RooDCBShape pdf,RooRealVar roobs, double MH){
      TCanvas* ccc = new TCanvas("ccc", "",1);
      ccc->cd();
      int nBins = 1000;
      double mean = MH;
      double  sigma  = mean*0.2;
      TH1F*  hist = (TH1F*) pdf.createHistogram("sigHist",roobs, RooFit::Binning(nBins,mean-4.*sigma,mean+4.*sigma) );
      double  halfMaxVal = 0.5*hist->GetMaximum();
      int  maxBin = hist->GetMaximumBin();
        
      int  binLeft,binRight;
      double xWidth,xLeft,xRight;
      
      int  last = 1;
      for(int ibin = 1; ibin<=maxBin;ibin++){
	  double   binVal = hist->GetBinContent(ibin);
	  if (binVal >= halfMaxVal){
	    binLeft = last;
	    break;
	  }
          if (binVal>0) last=ibin;
      }
      last = hist->GetXaxis()->GetNbins()+1;
      for(int ibin = hist->GetXaxis()->GetNbins()+1; ibin>maxBin;ibin--){ 
	double   binVal = hist->GetBinContent(ibin);
	if (binVal >= halfMaxVal){
	  binRight = last;
	  break;
	}
	if (binVal>0) last=ibin;
      }
      for(int ibin = 1; ibin<=maxBin;ibin++){ 
	double  binVal = hist->GetBinContent(ibin);
	if (binVal >= halfMaxVal){
	  binLeft = ibin;
	    break;
	}
      }
      for(int ibin = maxBin+1; ibin<=hist->GetXaxis()->GetNbins()+1;ibin++){  
	double binVal = hist->GetBinContent(ibin);
	if (binVal < halfMaxVal){
	  binRight = ibin-1;
	    break;
	}
      } 
      xWidth = 0.;
      if (binLeft > 0 && binRight > 0 ){
	xLeft = hist->GetXaxis()->GetBinCenter(binLeft);
	xRight = hist->GetXaxis()->GetBinCenter(binRight);
	xWidth = xRight-xLeft;
	std::cout<<TString::Format("FWHM = %f",xWidth)<<std::endl;
      }  else{
      std::cout<<"Did not succeed to compute the FWHM"<<std::endl;
      }
      bool plot =false;
      if(plot){
      hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinCenter(maxBin)-5*xWidth,hist->GetXaxis()->GetBinCenter(maxBin)+5*xWidth);
      hist->Draw("HIST");
       
      TLine* xL = new TLine(xLeft,0,xLeft,halfMaxVal);
      TLine* xR =  new TLine(xRight,0,xRight,halfMaxVal);
      xL->SetLineColor(kRed);
      xR->SetLineColor(kRed);
      xL->Draw("same");  
      xR->Draw("same");
     ccc->SaveAs("~/www/Pippone/test2.png");
      }
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      return xWidth;
    }







void plotAllSignalsRels(TString coupling, TString mass){
  std::cout<<"flag1"<<std::endl;
  TFile* f80x = TFile::Open("80X/ws_ResponseAndGen_M"+mass+"_k"+coupling+"_final.root");
  TFile* f76x_38T = TFile::Open("76X_38T/ws_ResponseAndGen_M"+mass+"_k"+coupling+"_final.root");
  TFile* f76x_0T = TFile::Open("76X_0T/ws_ResponseAndGen_M"+mass+"_k"+coupling+"_final.root");
 
  RooWorkspace* ws80x = (RooWorkspace*) f80x->Get("HLFactory_ws");
  RooWorkspace* ws76x_38T =  (RooWorkspace*)f76x_38T->Get("HLFactory_ws");
  RooWorkspace* ws76x_0T =  (RooWorkspace*)f76x_0T->Get("HLFactory_ws");

  TCanvas* c1 = new TCanvas("c1", "c1",1);
  c1->cd();

  TString svar;
  TString spdf;
  
  for(int c =0; c< NCAT+1; c++){
    if(c==2||c==3)continue;
    std::cout<<"flag1"<<std::endl;
    RooDCBShape* res80x = (RooDCBShape*) ws80x->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res76x_38T = (RooDCBShape*) ws76x_38T->pdf(TString::Format("dcbshape_cat%d",c));
    RooDCBShape* res76x_0T = (RooDCBShape*) ws76x_0T->pdf(TString::Format("dcbshape_cat%d",c));
   
    TH1F* h80x=new TH1F("h80x","",460, 300,6000);
    h80x->SetLineColor(kBlue);
    h80x->SetLineWidth(2);
    h80x->SetMarkerColor(kBlue);
    TH1F* h76x_38T=new TH1F("h76x_38T","",460, 300,6000);
    h76x_38T->SetLineColor(kMagenta);
    h76x_38T->SetLineWidth(2);
    h76x_38T->SetMarkerColor(kMagenta);
    TH1F* h76x_0T=new TH1F("h76x_0T","",460, 300,6000);
    h76x_0T->SetLineColor(kSpring+9);
    h76x_0T->SetLineWidth(2);
    h76x_0T->SetMarkerColor(kSpring+9);
    
    TLegend* legmc = new TLegend(0.58, 0.54, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(h80x," 80X","l");    
    legmc->AddEntry(h76x_38T,"76X 3.8T","l");    
    legmc->AddEntry(h76x_0T,"76X 0T","l");    
     
   

    std::cout<<"flag2"<<std::endl;
    RooPlot* pres = ws80x->var("var")->frame(Title("mass "),Bins(420));
    res80x->plotOn(pres,LineColor(kBlue));
    res76x_38T->plotOn(pres,LineColor(kMagenta));
    res76x_0T->plotOn(pres,LineColor(kSpring+9));
    std::cout<<"flag2"<<std::endl;
   
    std::cout<<"flag2"<<std::endl;
    pres->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    pres->GetYaxis()->SetTitle("a.u.");
    pres->Draw();
    legmc->Draw("same");
    c1->SaveAs(TString::Format("~/www/Pippone/CompareRelease_M_"+mass+"_k"+coupling+"_cat%d.png",c));
    c1->SaveAs(TString::Format("~/www/Pippone/CompareRelease_M_"+mass+"_kk"+coupling+"_cat%d.pdf",c));
    c1->SetLogy();
    c1->SaveAs(TString::Format("~/www/Pippone/CompareRelease_M_"+mass+"_kk"+coupling+"_cat%d_log.png",c));
    c1->SaveAs(TString::Format("~/www/Pippone/CompareRelease_M_"+mass+"_kk"+coupling+"_cat%d_log.pdf",c));
    c1->SetLogy(0);

  
}
}


# DiPhotonSignalModel


******************************************************************************************
STEP1: 
convert treees in the right format to be used for the signal model generation:
cd ConverTrees
source formaNtupleForSignalModel.sh
source formaNtupleForSignalModelGenOnly.sh
-> These two steps create the rootfiles w/ the trees for each cat
source makeAllHadd.sh
source makeallHddGen.sh

-> The output trees are the ones used to run the singal model machinery

******************************************************************************************
STEP2: 
run the signal model machinery

cd ProduceSignalmodel

->Here one needs two directories called as follows where the trees produced at the previous step ahve to be put accorndigly to their names:
GenSamples76X
RecoSamples80X

->Run the macro w/ the following command:
.L Produceworkspace.C

-> One has to call the runAllFits() macro which, for each width does:

  -> For each available reco mass:

    if makeWs==true:

        -> fit the response function w/ a DCB

        -> fit the gen level shape w/ a DCB

        -> write a workspace called: "ws_ResponseAndGen_M%d_k"+coupling+".root", iMass

    else:

        -> Read the previous ws

        -> make the convolution of the response and the gen level shape

        -> fit the convolution w/ a DCB

        -> throw an asimov dataset from the DCB 

        -> save it in: ("ws_ResponseAndGen_M%d_k"+coupling+"_final.root", iMass)

-> Call the plotAllSignalAsimov(coupling) funxtion which, for each width:

  -> reads the final ws

  -> look at all the final DCB shapes for each mass

  %compute the FHWM%

  -> compute the FHWM for each of them and parametrize it vs MH and fit it

  ->plot the FHWM vsMH and save the TF1 of the fit
  
  %compare asimov datasets and dcb%

  ->look at all the asimov dataset for each M

  ->plot for each MH the asimov and the DCB overlayed
  
  %parametrize DCB parametes vs MH%

  -> plot all the DCb paramters vs M and fit them
  
  %build parametric signal model%

  ->build rooformula with the paramterization of those paramters vs M

  ->include the systematics energy resolution

  ->plot the shapes to crosscheck
  
  %save everything in a ws%

  ->save the signal models for each category

  ->save all the paramterizations
    
******************************************************************************************
STEP 3: Run the machinery from Pasquale to get the datacards

->Look at the instructions in the diphoton repository

-> compile the RooDCBShape.cc class:
    
-> Put in the macro directory the root file from the previous step hadding all the widths in the same rootfile like this:
hadd SignalParametericShapes80X.root SignalParametericShapes80X_ws_kMpl001.root SignalParametericShapes80X_ws_kMpl01.root SignalParametericShapes80X_ws_kMpl02.root

Run the datacard production as follows (for sure there is a smarter way to do so):

root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include");
root [1] gROOT->GetInterpreter()->AddIncludePath("$ROOFITSYS/include"); 
root [2] gSystem->SetIncludePath("-I$ROOFITSYS/include"); 
root [3] gSystem->Load("$ROOFITSYS/lib/libRooFit.so") ;
root [4] gSystem->Load("$ROOFITSYS/lib/libRooFitCore.so") ; 
root [5] gSystem->Load("$ROOFITSYS/lib/libRooStats.so") ; 
root [6] .L RooDCBShape.cc+

common_opts="--plot-fit-bands --rescale-signal-to 1e-3 --minos-bands  --prepare-data  --do-parametric-signal-nuisances --generate-ws-bkgnbias-new"
dir38=full_analysis_spring16v1_sync_v4_cert_275125
coup=001

./combine_maker.sh $dir38 --lumi 3.99 --fit-name cic2 $common_opts --parametric-signal-new SignalParametericShapes80X_ws_kMpl001.root --parametric-signal-acceptance  acceptance_76_spin0.json --load lumi.json --only-coups ${coup} --label spin0_wnuis 

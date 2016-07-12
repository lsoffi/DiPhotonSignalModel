# DiPhotonSignalModel



STEP1: convert treees in the right format to be used for the signal model generation:
cd ConverTrees
source formaNtupleForSignalModel.sh
source formaNtupleForSignalModelGenOnly.sh
-> These two steps create the rootfiles w/ the trees for each cat
source makeAllHadd.sh
source makeallHddGen.sh

-> The output trees are the ones used to run the singal model machinery

STEP2: run the signal model machinery

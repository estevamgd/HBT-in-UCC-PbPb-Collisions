#include "../include/sig_dl_HF.h"  

void sigmix() {
    const char* inputFile = "data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root";  
    const char* treeName  = "demo/TreeMBUCC";

    int numSelectionVars = 7;
    int selectionVars[7] = { 1, 2, 10, 20, 60, 100, 140};
    int selectionVarsI = 3400;

    int poolSize = 10;   

    for (int i = 0; i < numSelectionVars - 1; i++){
        //sig_sl_utn(inputFile, treeName, selectionVars[i], selectionVars[i+1]);
    }
    
    sig_dl_HF(inputFile, treeName, selectionVarsI);

}
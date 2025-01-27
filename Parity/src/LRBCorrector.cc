/********************************************************************
File Name: LRBCorrector.cc

Created by: Michael Vallee
Email: mv836315@ohio.edu

Description:  This is the implementaion file of the LRBCorrector
              class, which is a child of the VQwDataHandler class.
              The functionality of this class is derived from
              LinRegBlue.

Last Modified: August 1, 2018 1:41 PM
********************************************************************/

#include <iostream>
using namespace std;

#include "QwHelicityPattern.h"

#include "LRBCorrector.h"

// Qweak headers
#include "VQwDataElement.h"
#include "QwVQWK_Channel.h"
#include "QwParameterFile.h"
#define MYSQLPP_SSQLS_NO_STATICS
#ifdef __USE_DATABASE__
#include "QwParitySSQLS.h"
#include "QwParityDB.h"
#endif // __USE_DATABASE__

#include <TFile.h>
#include <TH2.h>
#include <TTree.h> 

#include <TChain.h>
#include <TObjArray.h>
#include <TEventList.h> 

#include <TMatrixD.h>


LRBCorrector::LRBCorrector(QwOptions &options, QwHelicityPattern& helicitypattern, const TString &run) {
  
  run_label = run;
  ParseSeparator = "_";
  fEnableCorrection = false;
  ProcessOptions(options);
  LoadChannelMap(fRegressionMapFile);
  fHelicityPattern = &helicitypattern;
  QwSubsystemArrayParity& asym = helicitypattern.fAsymmetry;
  QwSubsystemArrayParity& diff = helicitypattern.fDifference;
  ConnectChannels(asym,diff);
  
}

/**
 * Defines configuration options using QwOptions functionality.
 * @param options Options object
 */
void LRBCorrector::DefineOptions(QwOptions &options)
{
  options.AddOptions("LRBCorrector")
    ("enable-lrbcorrection", po::value<bool>()->zero_tokens()->default_value(false),
     "enable lrb correction");
  options.AddOptions("LRBCorrector")
    ("lrbregression-map", po::value<std::string>()->default_value("regression_new.map"),
     "variables and sensitivities for lrb correction");
}

/**
 * Process configuration options using QwOptions functionality.
 * @param options Options object
 */
void LRBCorrector::ProcessOptions(QwOptions &options)
{
  fEnableCorrection = options.GetValue<bool>("enable-lrbcorrection");
  fRegressionMapFile = options.GetValue<std::string>("lrbregression-map");
  outPath = options.GetValue<std::string>("slope-file-path");
}


Int_t LRBCorrector::LoadChannelMap(const std::string& mapfile) {
  
  if (fEnableCorrection == false) {
    QwWarning << "enable-lrbcorrection is set to false.  Skipping LoadChannelMap for LRBCorrector" << QwLog::endl;
    return 0;
  }

  string TmpFilePath = run_label.Data();
  fRegressionMapFile = "blueR" + TmpFilePath + "new.slope.root";
  string MapFilePath = outPath + "/";
  string tmp = MapFilePath + fRegressionMapFile;
  TString corFileName(tmp.c_str());
  QwMessage << "Trying to open " << corFileName << QwLog::endl;
  TFile*  corFile=new TFile(corFileName);
  if( !corFile->IsOpen()) {
    printf("Failed to open %s, slopes NOT found\n",corFile->GetName());
    return 0;
  }

  TMatrixD *alphasM=0;
  alphasM=(TMatrixD *) corFile->Get("slopes");
  assert(alphasM);

  TH1 *dvnames = (TH1 *) corFile->Get("DVname");
  assert(dvnames);
  TH1 *ivnames = (TH1 *) corFile->Get("IVname");
  assert(ivnames);

  pair<EQwRegType, string> type_name_dv;
  pair<EQwRegType, string> type_name_iv;

  //  Loop through ivnames to get IV type and name
  //    Loop over # of dep variables
  //      Push-back the sensitiivity, IV type and IVnames into their respective vectors for each DV

  for (Int_t i = 0; i < dvnames->GetXaxis()->GetNbins(); ++i){
    type_name_dv = ParseRegressionVariable(dvnames->GetXaxis()->GetBinLabel(i+1));
    fDependentType.push_back(type_name_dv.first);
    fDependentName.push_back(type_name_dv.second);
  }

  fSensitivity.resize(fDependentType.size());

  for (Int_t i = 0; i < ivnames->GetXaxis()->GetNbins(); ++i) {
    type_name_iv = ParseRegressionVariable(ivnames->GetXaxis()->GetBinLabel(i+1));
    fIndependentType.push_back(type_name_iv.first);
    fIndependentName.push_back(type_name_iv.second);
    for (Int_t j = 0; j < dvnames->GetXaxis()->GetNbins(); ++j) {
      fSensitivity[j].push_back(-1.0*(*alphasM)(i,j));
    }
  }
  
  //printf("opened %s, slopes found, dump:\n",corFile->GetName());
  //alphasM->Print();
  corFile->Close();
  
}


Int_t LRBCorrector::ConnectChannels(
    QwSubsystemArrayParity& asym,
    QwSubsystemArrayParity& diff)
{
  VQwDataHandler::ConnectChannels(asym, diff);

  if (fEnableCorrection == false) {return 0;}

  // Add independent variables
  for (size_t iv = 0; iv < fIndependentName.size(); iv++) {
    // Get the independent variables
    const VQwHardwareChannel* iv_ptr = 0;
    //QwMessage << "fInpedententType[" << iv << "] = " << fIndependentType.at(iv) 
    //          << "; fInpedententName[" << iv << "] = " << fIndependentName.at(iv)
    //          << QwLog::endl;
    switch (fIndependentType.at(iv)) {
      case kRegTypeAsym:
        iv_ptr = asym.ReturnInternalValue(fIndependentName.at(iv));
        break;
      case kRegTypeDiff:
        iv_ptr = diff.ReturnInternalValue(fIndependentName.at(iv));
        break;
      default:
        QwWarning << "Independent variable for regression has unknown type."
                  << QwLog::endl;
        break;
    }
    if (iv_ptr) {
      //QwMessage << " iv: " << fIndependentName.at(iv) /*<< " (sens = " << fSensitivity.at(dv).at(iv) << ")"*/ << QwLog::endl;
      fIndependentVar.push_back(iv_ptr);
    } else {
      QwWarning << "Independent variable " << fIndependentName.at(iv) << " could not be found."
                << QwLog::endl;
    }
  }

  QwMessage << "In LRBCorrector::ConnectChannels; Number of IVs: " << fIndependentVar.size()
            << " Number of DVs: " << fDependentVar.size() << QwLog::endl;

}

void LRBCorrector::CalcOneOutput(const VQwHardwareChannel* dv, VQwHardwareChannel* output,
                                  vector< const VQwHardwareChannel* > &ivs,
                                  vector< Double_t > &sens) {
  
  // if second is NULL, can't do regression
  if (output == NULL){
    QwError<<"Second is value is NULL, unable to calculate regression."<<QwLog::endl;
    return;
  }
  // For correct type (asym, diff, mps)
  // if (fDependentType.at(dv) != type) continue;

  // Clear data in second, if first is NULL
  if (dv == NULL){
    output->ClearEventData();
  }else{
    // Update second value
    output->AssignValueFrom(dv);
  }

  // Add corrections
  for (size_t iv = 0; iv < ivs.size(); iv++) {
    output->ScaledAdd(sens.at(iv), ivs.at(iv));
  }
  
}


void LRBCorrector::ProcessData() {
  
  for (size_t i = 0; i < fDependentVar.size(); ++i) {
    CalcOneOutput(fDependentVar[i], fOutputVar[i], fIndependentVar, fSensitivity[i]);
  }
  
}


void LRBCorrector::LinearRegression(EQwRegType type)
{
  // Return if regression is not enabled
  if (! fEnableCorrection){
    QwDebug << "Regression is not enabled!" << QwLog::endl;
    return;
  }
  // Get error flag from QwHelicityPattern
  if (fHelicityPattern != NULL){
    fErrorFlag = fHelicityPattern->GetEventcutErrorFlag();
  } else if (fSubsystemArray != NULL){
    fErrorFlag = fSubsystemArray->GetEventcutErrorFlag();
  } else {
    QwError << "LRBCorrector::LinearRegression: Can't set fErrorFlag" << QwLog::endl;
    fErrorFlag = 0;
  }
  
  ProcessData();
}

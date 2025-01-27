/********************************************************************
File Name: QwCorrelator.cc

Created by: Michael Vallee
Email: mv836315@ohio.edu

Description:  This is the implementation file of the QwCorrelator
              class, which is a child of the VQwDataHandler class.
              The functionality of this class is derived from
              LinRegBlue.

Last Modified: August 1, 2018 1:43 PM
********************************************************************/

#include <iostream>
using namespace std;

#include "QwHelicityPattern.h"

#include "QwCorrelator.h"

// Qweak headers
#include "VQwDataElement.h"
#include "QwVQWK_Channel.h"
#include "QwParameterFile.h"
#define MYSQLPP_SSQLS_NO_STATICS
#ifdef __USE_DATABASE__
#include "QwParitySSQLS.h"
#include "QwParityDB.h"
#endif // __USE_DATABASE__


//******************************************************************************************************************************************************

QwCorrelator::QwCorrelator(QwOptions &options, QwHelicityPattern& helicitypattern, const TString &run): corA("input") {

  run_label = run;
  ParseSeparator = "_";
  fEnableCorrelation = false;
  ProcessOptions(options);
  init("blueReg.conf");
  corA.SetDisableHistogramFlag(fDisableHistos);
  QwSubsystemArrayParity& asym = helicitypattern.fAsymmetry;
  QwSubsystemArrayParity& diff = helicitypattern.fDifference;
  ConnectChannels(asym, diff);

  if (fEnableCorrelation == true) {

    vector<TString> fIndependentName_t;
    vector<TString> fDependentName_t;
 
    for (size_t i = 0; i < fIndependentName.size(); ++i) {
      fIndependentName_t.push_back(TString(fIndependentName.at(i)));
    }
    for (size_t i = 0; i < fDependentName.size(); ++i) {
      fDependentName_t.push_back(TString(fDependentName.at(i)));
    }
 
    corA.init(fIndependentName_t, fDependentName_t);

  }
  
}

void QwCorrelator::init(const std::string configFName) {
	
  LoadChannelMap(configFName);

}


void QwCorrelator::FillCorrelator() {

  if (! fEnableCorrelation) return;

  UInt_t error;

  for (size_t i = 0; i < fDependentVar.size(); ++i) {
    error |= fDependentVar.at(i)->GetErrorCode();
    fDependentValues.at(i) = (fDependentVar[i]->GetValue());
    //QwMessage << "Loading DV " << fDependentVar.at(i) << " into fDependentValues." << QwLog::endl;
  }
  for (size_t i = 0; i < fIndependentVar.size(); ++i) {
    error |= fIndependentVar.at(i)->GetErrorCode();
    fIndependentValues.at(i) = (fIndependentVar[i]->GetValue());
    //QwMessage << "Loading IV " << fIndependentVar.at(i) << " into fIndependentValues." << QwLog::endl;
  }

  //QwMessage << "fDependentVar has a size of: " << fDependentVar.size() << QwLog::endl;
  //QwMessage << "fIndependentVar has a size of: " << fIndependentVar.size() << QwLog::endl;

  if (error == 0) {
    corA.addEvent(&fIndependentValues[0],&fDependentValues[0]);
  }
  
}


void QwCorrelator::CalcCorrelations() {

  if (! fEnableCorrelation) return;

	corA.finish();
	
  std::string TmpRunLabel = run_label.Data();
  std::string fSlopeFileName = "blueR" + TmpRunLabel + "new.slope.root";
  std::string fSlopeFilePath = fAlphaOutputPath + "/";
  std::string tmp = fSlopeFilePath + fSlopeFileName;

  TString outAlphas=Form(tmp.c_str());
  corA.exportAlphas(outAlphas, fIndependentFull, fDependentFull);
  corA.exportAlias(fAliasOutputPath, "/regalias_"+run_label, fIndependentFull, fDependentFull);

}


//******************************************************************************************************************************************************


void QwCorrelator::DefineOptions(QwOptions &options) {
	
  options.AddOptions("Correlator")
    ("enable-correlator", po::value<bool>()->zero_tokens()->default_value(false),
     "enables correlator");
  options.AddOptions("Correlator")
    ("correlator-map", po::value<std::string>()->default_value("regression.map"),
     "variables and sensitivities for regression");
  
  options.AddOptions("Correlator")
    ("slope-file-path", po::value<std::string>()->default_value("."),
     "Path for the slop files, also used for lrb (alpha)");
  options.AddOptions("Correlator")
    ("alias-output-path", po::value<std::string>()->default_value("."),
     "Path for the final correlation output file (alias)");
	
}


void QwCorrelator::ProcessOptions(QwOptions &options) {
	
  fEnableCorrelation = options.GetValue<bool>("enable-correlator");
  fCorrelatorMapFile = options.GetValue<std::string>("correlator-map");
	  
  fAlphaOutputPath = options.GetValue<std::string>("slope-file-path");
  fAliasOutputPath = options.GetValue<std::string>("alias-output-path");

  fDisableHistos = options.GetValue<bool>("disable-histos");
	  
}


/** Load the channel map
 *
 * @param mapfile Filename of map file
 * @return Zero when success
 */
Int_t QwCorrelator::LoadChannelMap(const std::string& mapfile)
{
  // Return if regression is not enabled
  if (! fEnableCorrelation) return 0;

  // Open the file
  QwParameterFile map(mapfile);

  // Read the sections of dependent variables
  std::pair<EQwRegType,std::string> type_name;

  // Add independent variables and sensitivities
  while (map.ReadNextLine()) {
    // Throw away comments, whitespace, empty lines
    map.TrimComment();
    map.TrimWhitespace();
    if (map.LineIsEmpty()) continue;
    // Get first token: label (dv or iv), second token is the name like "asym_blah"
    string primary_token = map.GetNextToken(" ");
    string current_token = map.GetNextToken(" ");
    // Parse current token into independent variable type and name
    type_name = ParseRegressionVariable(current_token);

    if (primary_token == "iv") {
      fIndependentType.push_back(type_name.first);
      fIndependentName.push_back(type_name.second);
      fIndependentFull.push_back(current_token);
      //QwMessage << "IV Type: " << type_name.first << QwLog::endl;
      //QwMessage << "IV Name: " << type_name.second << QwLog::endl;
    }
    else if (primary_token == "dv") {
      fDependentType.push_back(type_name.first);
      fDependentName.push_back(type_name.second);
      fDependentFull.push_back(current_token);
      //QwMessage << "DV Type: " << type_name.first << QwLog::endl;
      //QwMessage << "DV Name: " << type_name.second << QwLog::endl;
   	}
    else if (primary_token == "treetype") {
      QwMessage << "Tree Type read, ignoring." << QwLog::endl;
    }
    else {
   	  QwError << "Function LoadChannelMap in QwCorrelator.cc read in invalid primary_token." << QwLog::endl;
    }
  }
  
  //QwMessage << "fIndependentType has a size of: " << fIndependentType.size() << QwLog::endl;
  //QwMessage << "fIndependentName has a size of: " << fIndependentName.size() << QwLog::endl;
  //QwMessage << "fDependentType has a size of: " << fDependentType.size() << QwLog::endl;
  //QwMessage << "fDependentName has a size of: " << fDependentName.size() << QwLog::endl;

}


Int_t QwCorrelator::ConnectChannels(QwSubsystemArrayParity& asym, QwSubsystemArrayParity& diff) {
	
	// Return if regression is not enabled
  if (! fEnableCorrelation) return 0;

  /// Fill vector of pointers to the relevant data elements
  for (size_t dv = 0; dv < fDependentName.size(); dv++) {
    // Get the dependent variables

    VQwHardwareChannel* dv_ptr = 0;
    QwVQWK_Channel* new_vqwk = NULL;
    QwVQWK_Channel* vqwk = NULL;
    string name = "";
    string reg = "reg_";
    
    if (fDependentType.at(dv)==kRegTypeMps){
      //  Quietly ignore the MPS type when we're connecting the asym & diff
      continue;
    } else if(fDependentName.at(dv).at(0) == '@' ){
      name = fDependentName.at(dv).substr(1,fDependentName.at(dv).length());
    }else{
      switch (fDependentType.at(dv)) {
        case kRegTypeAsym:
          dv_ptr = asym.ReturnInternalValueForFriends(fDependentName.at(dv));
          break;
        case kRegTypeDiff:
          dv_ptr = diff.ReturnInternalValueForFriends(fDependentName.at(dv));
          break;
        default:
          QwWarning << "QwCombiner::ConnectChannels(QwSubsystemArrayParity& asym, QwSubsystemArrayParity& diff):  Dependent variable, "
	          	      << fDependentName.at(dv)
		                << ", for asym/diff regression does not have proper type, type=="
		                << fDependentType.at(dv) << "."<< QwLog::endl;
          break;
        }

      vqwk = dynamic_cast<QwVQWK_Channel*>(dv_ptr);
      name = vqwk->GetElementName().Data();
      name.insert(0, reg);
      new_vqwk = new QwVQWK_Channel(*vqwk, VQwDataElement::kDerived);
      new_vqwk->SetElementName(name);
    }

    // alias
    if(fDependentName.at(dv).at(0) == '@'){
      //QwMessage << "dv: " << name << QwLog::endl;
      new_vqwk = new QwVQWK_Channel(name, VQwDataElement::kDerived);
    }
    // defined type
    else if(dv_ptr!=NULL){
      //QwMessage << "dv: " << fDependentName.at(dv) << QwLog::endl;
    }else {
      QwWarning << "Dependent variable " << fDependentName.at(dv) << " could not be found, "
                << "or is not a VQWK channel." << QwLog::endl;
      continue; 
    }

    // pair creation
    if(vqwk != NULL){
      // fDependentVarType.push_back(fDependentType.at(dv));
      fDependentVar.push_back(vqwk);
    }

  }
  
    // Add independent variables
    // fIndependentVar.resize(fDependentVar.size());
  for (size_t iv = 0; iv < fIndependentName.size(); iv++) {
    // Get the independent variables
    const VQwHardwareChannel* iv_ptr = 0;
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
      //QwMessage << " iv: " << fIndependentName.at(iv) << QwLog::endl;
      fIndependentVar.push_back(iv_ptr);
    } else {
      QwWarning << "Independent variable " << fIndependentName.at(iv) << " for regression of could not be found."
                << QwLog::endl;
    }
      
  }
  fIndependentValues.resize(fIndependentVar.size());
  fDependentValues.resize(fDependentVar.size());
  return 0;
	
}




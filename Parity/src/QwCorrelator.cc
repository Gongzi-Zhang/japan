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

// System headers
#include <iostream>

// ROOT headers
#include <TFile.h>
#include <TH2D.h>

// Qweak headers
#include "QwCorrelator.h"
#include "QwHelicityPattern.h"
#include "VQwDataElement.h"
#include "QwVQWK_Channel.h"
#include "QwParameterFile.h"
#define MYSQLPP_SSQLS_NO_STATICS
#ifdef __USE_DATABASE__
#include "QwParitySSQLS.h"
#include "QwParityDB.h"
#endif // __USE_DATABASE__

// Register this handler with the factory
RegisterHandlerFactory(QwCorrelator);

//******************************************************************************************************************************************************

QwCorrelator::QwCorrelator(const TString& name)
: VQwDataHandler(name),
  fDisableHistos(true),
  fAlphaOutputPath("."),
  fAliasOutputPath("."),
  mCore("input"),
  nP(0), nY(0),
  h1iv(0),h2iv(0),
  h1dv(0),h2dv(0)
{
  ParseSeparator = "_";
  fTotalCount = 0;
  fGoodCount  = 0;
}

QwCorrelator::~QwCorrelator()
{
  QwMessage << "destructor QwCorrelator=" << mCore << QwLog::endl;

  // only if previously allocated
  if (h1iv) {
    delete  [] h1iv;
    delete  [] h2iv;
    delete  [] h1dv;
    delete  [] h2dv;
  }

  QwMessage << "destructor QwCorrelator done" << QwLog::endl;
}


void QwCorrelator::ParseConfigFile(QwParameterFile& file)
{
  VQwDataHandler::ParseConfigFile(file);
  file.PopValue("slope-path", fAlphaOutputPath);
  file.PopValue("alias-path", fAliasOutputPath);
  file.PopValue("disable-histos", fDisableHistos);
}

void QwCorrelator::AccumulateRunningSum()
{
  UInt_t error = 0;

  fTotalCount++;

  for (size_t i = 0; i < fDependentVar.size(); ++i) {
    error |= fDependentVar.at(i)->GetErrorCode();
    fDependentValues.at(i) = (fDependentVar[i]->GetValue());
    QwVerbose << "Loading DV " << fDependentVar.at(i)
              << " into fDependentValues." << QwLog::endl;
    if ( fDependentVar.at(i)->GetErrorCode() !=0)  (fErrCounts_DV.at(i))++;
  }
  for (size_t i = 0; i < fIndependentVar.size(); ++i) {
    error |= fIndependentVar.at(i)->GetErrorCode();
    fIndependentValues.at(i) = (fIndependentVar[i]->GetValue());
    QwVerbose << "Loading IV " << fIndependentVar.at(i)
              << " into fIndependentValues." << QwLog::endl;
    if ( fIndependentVar.at(i)->GetErrorCode() !=0)  (fErrCounts_IV.at(i))++;
  }

  QwVerbose << "fDependentVar has a size of: " << fDependentVar.size() << QwLog::endl;
  QwVerbose << "fIndependentVar has a size of: " << fIndependentVar.size() << QwLog::endl;

  if (error == 0) {
    fGoodCount++;
    addEvent(&fIndependentValues[0],&fDependentValues[0]);
  }
}


void QwCorrelator::CalcCorrelations()
{
  QwMessage << "QwCorrelator:  Total entries: " << fTotalCount
            << ", good entries: "<< fGoodCount << QwLog::endl;
  for (size_t i = 0; i < fDependentVar.size(); ++i) {
    if (fErrCounts_DV.at(i) > 0)
      QwMessage << "   Entries failed due to "
                << fDependentVar.at(i)->GetElementName()
                << ": " <<  fErrCounts_DV.at(i) << QwLog::endl;
  }
  for (size_t i = 0; i < fIndependentVar.size(); ++i) {
    if (fErrCounts_IV.at(i) > 0)
      QwMessage << "   Entries failed due to "
                << fIndependentVar.at(i)->GetElementName()
                << ": " <<  fErrCounts_IV.at(i) << QwLog::endl;
  }

  if (linReg.failed()) {
    QwMessage << " abnormal finish of linReg" << QwLog::endl;
    return;
  }

  linReg.printSummaryP();
  linReg.printSummaryY();
  linReg.solve();
  linReg.printSummaryAlphas();

  std::string RunLabel = run_label.Data();
  std::string slope_file_name = "blueR" + RunLabel + "new.slope.root";
  std::string slope_file_path = fAlphaOutputPath + "/";
  std::string slope_file = slope_file_path + slope_file_name;
  exportAlphas(slope_file, fIndependentFull, fDependentFull);

  std::string alias_file_path = fAlphaOutputPath + "/";
  std::string alias_file_name = "regalias_" + RunLabel;
  exportAlias(alias_file_path, alias_file_name, fIndependentFull, fDependentFull);
}


/** Load the channel map
 *
 * @param mapfile Filename of map file
 * @return Zero when success
 */
Int_t QwCorrelator::LoadChannelMap(const std::string& mapfile)
{
  // Open the file
  QwParameterFile map(mapfile);

  // Read the sections of dependent variables
  std::pair<EQwHandleType,std::string> type_name;

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
    type_name = ParseHandledVariable(current_token);

    if (primary_token == "iv") {
      fIndependentType.push_back(type_name.first);
      fIndependentName.push_back(type_name.second);
      fIndependentFull.push_back(current_token);
      QwVerbose << "IV Type: " << type_name.first << QwLog::endl;
      QwVerbose << "IV Name: " << type_name.second << QwLog::endl;
      QwVerbose << "IV Full: " << current_token << QwLog::endl;
    }
    else if (primary_token == "dv") {
      fDependentType.push_back(type_name.first);
      fDependentName.push_back(type_name.second);
      fDependentFull.push_back(current_token);
      QwVerbose << "DV Type: " << type_name.first << QwLog::endl;
      QwVerbose << "DV Name: " << type_name.second << QwLog::endl;
      QwVerbose << "DV Full: " << current_token << QwLog::endl;
    }
    else if (primary_token == "treetype") {
      QwMessage << "Tree Type read, ignoring." << QwLog::endl;
    }
    else {
   	  QwError << "Function LoadChannelMap in QwCorrelator.cc read in invalid primary_token." << QwLog::endl;
    }
  }
  
  QwVerbose << "fIndependentType has a size of: " << fIndependentType.size() << QwLog::endl;
  QwVerbose << "fIndependentName has a size of: " << fIndependentName.size() << QwLog::endl;
  QwVerbose << "fDependentType has a size of: " << fDependentType.size() << QwLog::endl;
  QwVerbose << "fDependentName has a size of: " << fDependentName.size() << QwLog::endl;

  return 0;
}


Int_t QwCorrelator::ConnectChannels(QwSubsystemArrayParity& asym, QwSubsystemArrayParity& diff)
{
  /// Fill vector of pointers to the relevant data elements
  for (size_t dv = 0; dv < fDependentName.size(); dv++) {
    // Get the dependent variables

    VQwHardwareChannel* dv_ptr = 0;
    QwVQWK_Channel* vqwk = NULL;
    string name = "";
    string reg = "reg_";
    
    if (fDependentType.at(dv)==kHandleTypeMps){
      //  Quietly ignore the MPS type when we're connecting the asym & diff
      continue;
    } else if(fDependentName.at(dv).at(0) == '@' ){
      name = fDependentName.at(dv).substr(1,fDependentName.at(dv).length());
    } else {
      switch (fDependentType.at(dv)) {
        case kHandleTypeAsym:
          dv_ptr = asym.ReturnInternalValueForFriends(fDependentName.at(dv));
          break;
        case kHandleTypeDiff:
          dv_ptr = diff.ReturnInternalValueForFriends(fDependentName.at(dv));
          break;
        default:
          QwWarning << "QwCombiner::ConnectChannels(QwSubsystemArrayParity& asym, QwSubsystemArrayParity& diff):  Dependent variable, "
                                << fDependentName.at(dv)
		                << ", for asym/diff correlator does not have proper type, type=="
		                << fDependentType.at(dv) << "."<< QwLog::endl;
          break;
        }

      vqwk = dynamic_cast<QwVQWK_Channel*>(dv_ptr);
      name = vqwk->GetElementName().Data();
      name.insert(0, reg);
    }

    // pair creation
    if (vqwk != NULL) {
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
      case kHandleTypeAsym:
        iv_ptr = asym.ReturnInternalValue(fIndependentName.at(iv));
        break;
      case kHandleTypeDiff:
        iv_ptr = diff.ReturnInternalValue(fIndependentName.at(iv));
        break;
      default:
        QwWarning << "Independent variable for correlator has unknown type."
                  << QwLog::endl;
        break;
    }
    if (iv_ptr) {
      QwVerbose << " iv: " << fIndependentName.at(iv) << QwLog::endl;
      fIndependentVar.push_back(iv_ptr);
    } else {
      QwWarning << "Independent variable " << fIndependentName.at(iv) << " for correlator could not be found."
                << QwLog::endl;
    }
  }

  fIndependentValues.resize(fIndependentVar.size());
  fDependentValues.resize(fDependentVar.size());

  nP = fIndependentName.size();
  nY = fDependentName.size();

  initHistos(fIndependentName,fDependentName);

  linReg.setDims(nP, nY);
  linReg.init();

  fErrCounts_IV.resize(fIndependentVar.size(),0);
  fErrCounts_DV.resize(fDependentVar.size(),0);

  return 0;
}


void QwCorrelator::addEvent(double *Pvec, double *Yvec)
{
  linReg.accumulate(Pvec, Yvec);

  // .... monitoring

  if (fDisableHistos == false) {
    for (int i = 0; i < nP; i++) {
      h1iv[i]->Fill(Pvec[i]);
      for (int j = i + 1; j < nP; j++)
        h2iv[i*nP+j]->Fill(Pvec[i],Pvec[j]);
    }
    for (int j = 0; j < nY; j++) {
      h1dv[j]->Fill(Yvec[j]);
      for (int i = 0; i < nP; i++)
        h2dv[i*nY+j]->Fill(Pvec[i],Yvec[j]);
    }
  }
}

void QwCorrelator::initHistos(std::vector<std::string> Pname, std::vector<std::string> Yname)
{
  QwMessage << "QwCorrelator::initHistos()" << QwLog::endl;

  //..... 1D,  iv
  h1iv = new TH1 *[nP];
  for (int i = 0; i < nP; i++) {
    TH1* h = h1iv[i] = new TH1D(Form(mCore+"P%d",i),Form("iv P%d=%s, pass=%s ;iv=%s (ppm)",i,Pname[i].c_str(),mCore.Data(),Pname[i].c_str()),128,0.,0.);
    h->GetXaxis()->SetNdivisions(4);
  }

  double x1 = 0;
  //..... 2D,  iv correlations
  h2iv = new TH1 *[nP*nP]; // not all are used
  for (int i = 0; i < nP; i++) {
    for (int j = i + 1; j < nP; j++) {
      TH1* h = h2iv[i*nP+j] = new TH2D(Form(mCore+"P%d_P%d",i,j),Form("iv correlation  P%d_P%d, pass=%s ;P%d=%s (ppm);P%d=%s   (ppm)  ",i,j,mCore.Data(),i,Pname[i].c_str(),j,Pname[j].c_str()),64,-x1,x1,64,-x1,x1);
      h->GetXaxis()->SetTitleColor(kBlue);
      h->GetYaxis()->SetTitleColor(kBlue);
      h->GetXaxis()->SetNdivisions(4);
      h->GetYaxis()->SetNdivisions(4);
    }
  }

  //..... 1D,  dv
  h1dv = new TH1 *[nY];
  for (int i = 0; i < nY; i++) {
    TH1* h = h1dv[i] = new TH1D(Form(mCore+"Y%d",i),Form("dv Y%d=%s, pass=%s ;dv=%s (ppm)",i,Yname[i].c_str(),mCore.Data(),Yname[i].c_str()),128,0.,0.);
    h->GetXaxis()->SetNdivisions(4);
  }

  double y1=0;
  //..... 2D,  dv-iv correlations
  h2dv = new TH1 *[nP*nY]; // not all are used
  for (int i = 0; i < nP; i++) {
    for (int j = 0; j < nY; j++) {
      TH1* h = h2dv[i*nY+j] = new TH2D(Form(mCore+"P%d_Y%d",i,j),Form("iv-dv correlation  P%d_Y%d, pass=%s ;P%d=%s (ppm);Y%d=%s   (ppm)  ",i,j,mCore.Data(),i,Pname[i].c_str(),j,Yname[j].c_str()),64,-x1,x1,64,-y1,y1);
      h->GetXaxis()->SetTitleColor(kBlue);
      h->GetYaxis()->SetTitleColor(kBlue);
      h->GetXaxis()->SetNdivisions(4);
      h->GetYaxis()->SetNdivisions(4);
    }
  }

  // store list of names to be archived
  hA[0] = new TH1D(mCore+"NamesIV",Form("IV name list nIV=%d",nP),nP,0,1);
  for (int i = 0; i < nP; i++)
    hA[0]->Fill(TString(Pname[i]),1.*i);
  hA[1] = new TH1D(mCore+"NamesDV",Form("DV name list nIV=%d",nY),nY,0,1);
  for (int i = 0; i < nY; i++)
    hA[1]->Fill(TString(Yname[i]),i*1.);
}

void QwCorrelator::exportAlphas(
    const std::string& outName,
    const std::vector<std::string>& ivName,
    const std::vector<std::string>& dvName)
{
  QwMessage << "::::::::::::::::QwCorrelator::exportAlphas(" << outName << ") :::::::::::" << QwLog::endl;

  TFile* hFile = new TFile(TString(outName),"RECREATE","correlation coefficents");
  linReg.mA.Write("slopes");
  linReg.mAsig.Write("sigSlopes");
  linReg.mRjk.Write("IV_correlation");
  linReg.mMP.Write("IV_mean");
  linReg.mMY.Write("DV_mean");

  // add processed matrices
  double usedEve = linReg.getUsedEve();
  TMatrixD Mstat(1,1);
  Mstat(0,0) = usedEve;
  Mstat.Write("MyStat");

  //... IVs
  TMatrixD MsigIV(nP,1);
  TH1D hiv("IVname","names of IVs",nP,-0.5,nP-0.5);
  double val;
  for (int i = 0; i < nP; i++) {
    Int_t testval = linReg.getSigmaP(i,val);
    assert(testval==0);
    MsigIV(i,0) = val;
    hiv.Fill(TString(ivName[i]),i);
  }
  MsigIV.Write("IV_sigma"); // of distribution
  hiv.Write();

  //... DVs
  TMatrixD MsigDV(nY,1);
  TH1D hdv("DVname","names of IVs",nY,-0.5,nY-0.5);
  for (int i = 0; i < nY; i++) {
    Int_t testval = linReg.getSigmaY(i,val);
    assert(testval==0);
    MsigDV(i,0) = val;
    hdv.Fill(TString(dvName[i]),i);
  }
  MsigDV.Write("DV_sigma"); // of distribution
  hdv.Write();

  //raw matrices
  linReg.mVPP.Write("IV_rawVariance");
  linReg.mVPY.Write("IV_DV_rawVariance");
  linReg.mVY2.Write("DV_rawVariance");

  hFile->Close();

  QwMessage << "saved " << hFile->GetName() << QwLog::endl;
}

void QwCorrelator::exportAlias(
    const std::string& outPath,
    const std::string& macroName,
    const std::vector<std::string>& Pname,
    const std::vector<std::string>& Yname) const
{
  QwMessage << "::::::::::::::::QwCorrelator::exportAlias(" << macroName << ") :::::::::::" << QwLog::endl;

  std::ofstream file(outPath + macroName + ".C");

  file << "void " << macroName << "() {" << std::endl;
  file << "  TTree* tree = (TTree*) gDirectory->Get(\"mul\");" << std::endl;
  for (int iy = 0; iy < nY; iy++) {
    file << "  tree->SetAlias(\"reg_" << Yname[iy] << "\"," << std::endl << "         \"" << Yname[iy];
    for (int j = 0; j < nP; j++) {
      double val= -linReg.mA(j,iy);
      if (val > 0)
        file << "+";
      file << std::scientific << val << "*" << Pname[j];
    }
    file << "\");" << std::endl;
  }

  file << "}" << std::endl;
  file.close();

  QwMessage << "saved " << macroName << QwLog::endl;
}

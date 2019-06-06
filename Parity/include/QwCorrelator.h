/********************************************************************
File Name: QwCorrelator.h

Created by: Michael Vallee
Email: mv836315@ohio.edu

Description:  This is the header file of the QwCorrelator class,
              which is a child of the VQwDataHandler class.  The
              functionality of this class is derived from
              LinRegBlue.

Last Modified: August 1, 2018 1:43 PM
********************************************************************/

#ifndef QWCORRELATOR_H_
#define QWCORRELATOR_H_

// Parent Class
#include "VQwDataHandler.h"

// LinRegBlue Correlator Class
#include "LinReg_Bevington_Pebay.h"

class QwCorrelator : public VQwDataHandler, public MQwDataHandlerCloneable<QwCorrelator>{
 public:
  /// \brief Constructor with name
  QwCorrelator(const TString& name);
  /// \brief Virtual destructor
  virtual ~QwCorrelator();

  void ParseConfigFile(QwParameterFile& file);

  Int_t LoadChannelMap(const std::string& mapfile);
  	
  void readConfig(const char * configFName);
  	
  /// \brief Connect to Channels (asymmetry/difference only)
  Int_t ConnectChannels(QwSubsystemArrayParity& asym, QwSubsystemArrayParity& diff);
		
  void unpackEvent();

  void AccumulateRunningSum();
  void ProcessData(){};
  void CalcCorrelations();
		
 protected:
  bool fDisableHistos;
  
  std::vector< std::string > fIndependentFull;
  std::vector< std::string > fDependentFull;
    
  //  Using the fDependentType and fDependentName from base class, but override the IV arrays
  std::vector< EQwHandleType > fIndependentType;
  std::vector< std::string > fIndependentName;

  std::vector< const VQwHardwareChannel* > fIndependentVar;
  std::vector< Double_t > fIndependentValues;

  std::string fAlphaOutputPath;
  std::string fAliasOutputPath;		

  Int_t fTotalCount;
  Int_t fGoodCount;
  std::vector< Int_t > fErrCounts_IV;
  std::vector< Int_t > fErrCounts_DV;

 private:

  TString mCore;
  int nP, nY;

  // histograms
  enum {mxHA=4};
  TH1 * hA[mxHA];

  // monitoring histos for iv & dv
  TH1 ** h1iv, **h2iv, ** h1dv, **h2dv;
  void initHistos(std::vector<std::string> ivName, std::vector<std::string> dvName);

  LinRegBevPeb linReg;

 public:

  void addEvent(double *Pvec, double *Yvec);

  void exportAlphas(
      const std::string& outPath,
      const std::vector<std::string>& ivName,
      const std::vector<std::string>& dvName);

  void exportAlias(
      const std::string& outPath,
      const std::string& macroName,
      const std::vector<std::string>& ivName,
      const std::vector<std::string>& dvName) const;

};


#endif //QWCORRELATOR_H_

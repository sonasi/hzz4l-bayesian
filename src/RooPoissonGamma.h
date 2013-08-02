/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/
#ifdef __MAKECINT__ 
#pragma link C++ myclass+; 
#endif

#ifndef ROOPOISSONGAMMA
#define ROOPOISSONGAMMA

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
 
class RooPoissonGamma : public RooAbsPdf {
public:
  RooPoissonGamma() {} ; 
  RooPoissonGamma(const char *name, const char *title,
                RooAbsReal& x, const RooArgList& meanList, const RooArgList& errorList, Bool_t noRounding=kFALSE);
  RooPoissonGamma(const RooPoissonGamma& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPoissonGamma(*this,newname); }
  inline virtual ~RooPoissonGamma() { }

//  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ; 
//  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);


protected:

  RooRealProxy _x ;
  RooListProxy _meanList ;
  RooListProxy _errorList ;
  TIterator* _meanIter ;
  TIterator* _errorIter ;
  Bool_t  _noRounding ;
  Bool_t  _protectNegative ;

  ROOT::Math::Random<ROOT::Math::GSLRngMT>* _gslRan;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooPoissonGamma,0) // Poisson function marginalized over statistical errors
};
 
#endif

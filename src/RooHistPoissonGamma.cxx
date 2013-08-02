/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofit:$Id: RooHistPoissonGamma.cxx 44982 2012-07-10 08:36:13Z moneta $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// RooHistPoissonGamma implements a probablity density function sampled from a 
// multidimensional histogram. The histogram distribution is explicitly
// normalized by RooHistPoissonGamma and can have an arbitrary number of real or 
// discrete dimensions.
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"

#include "RooHistPoissonGamma.h"
#include "RooDataHist.h"

#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooRandom.h"

#include <math.h> 
#include "TMath.h" 
#include "Math/ProbFuncMathCore.h"


using namespace std;

ClassImp(RooHistPoissonGamma)



//_____________________________________________________________________________
RooHistPoissonGamma::RooHistPoissonGamma() : _dataHist(0), _totVolume(0), _unitNorm(kFALSE)
{
  // Default constructor
  // coverity[UNINIT_CTOR]
  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _pdfObsList.createIterator() ;
}




//_____________________________________________________________________________
RooHistPoissonGamma::RooHistPoissonGamma(const char *name, const char *title, const RooArgSet& vars, 
		       const RooDataHist& dataHist, const RooArgList& histList, const RooArgList& muList, Int_t intOrder) :
  RooAbsPdf(name,title), 
  _pdfObsList("pdfObs","List of p.d.f. observables",this),
  _histList("histList","List of mus for sources",this),
  _muList("muList","List of means for sources",this),
  _codeReg(10),
  _gslRan(new ROOT::Math::Random<ROOT::Math::GSLRngMT>()),
  _dataHist((RooDataHist*)&dataHist), 
  _intOrder(intOrder),
  _cdfBoundaries(kFALSE),
  _totVolume(0),
  _unitNorm(kFALSE)
{
  // Constructor from a RooDataHist. RooDataHist dimensions
  // can be either real or discrete. See RooDataHist::RooDataHist for details on the binning.
  // RooHistPoissonGamma neither owns or clone 'dhist' and the user must ensure the input histogram exists
  // for the entire life span of this PDF.


  _histObsList.addClone(vars) ;
  _pdfObsList.add(vars) ;



  // Verify that vars and dataHist.get() have identical contents
  const RooArgSet* dvars = dataHist.get() ;
  if (vars.getSize()!=dvars->getSize()) {
    coutE(InputArguments) << "RooHistPoissonGamma::ctor(" << GetName() 
			  << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
    assert(0) ;
  }
  TIterator* iter = vars.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dvars->find(arg->GetName())) {
      coutE(InputArguments) << "RooHistPoissonGamma::ctor(" << GetName() 
			    << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
      assert(0) ;
    }
  }
  delete iter ;


    _muIter = _muList.createIterator();
    TIterator* muIter = muList.createIterator() ;
    RooAbsArg* mu ;
    while((mu = (RooAbsArg*)muIter->Next())) {
        if (!dynamic_cast<RooAbsReal*>(mu)) {
            cout << "RooPolynomial::ctor(" << GetName() << ") ERROR: muficient " << mu->GetName() 
                 << " is not of type RooAbsReal" << endl ;
            assert(0) ;
        }
        _muList.add(*mu) ;
    }
    delete muIter ;


  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _pdfObsList.createIterator() ;
  
  
    _histIter = _histList.createIterator();
    TIterator* histIter = histList.createIterator() ;
    RooHistPdf* histPdf ;
    while((histPdf = (RooHistPdf*)histIter->Next())) {
        if (!dynamic_cast<RooHistPdf*>(histPdf)) {
            cout << "RooPolynomial::ctor(" << GetName() << ") ERROR: muficient " << histPdf->GetName() 
                 << " is not of type RooAbsReal" << endl ;
            assert(0) ;
        }


            RooDataHist* dhist;
            dhist = &histPdf->dataHist();

            if(_dataHist->numEntries() != dhist->numEntries()) {
            coutE(InputArguments) << "Data histograms (" << _dataHist->GetName() << ") and input histogram (" << dhist->GetName() << ") must contain the same number of bins (" << _dataHist->numEntries()  << " vs. " <<  dhist->numEntries() << "). " << endl ;
            cout  << endl;

            assert(0) ;
            }

            // Verify that vars and dhist.get() have identical contents
            const RooArgSet* dvars = dataHist.get() ;
            if (vars.getSize()!=dvars->getSize()) {
            coutE(InputArguments) << "RooHistPoissonGamma::ctor(" << GetName() 
            << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
            assert(0) ;
            }
            TIterator* iter = vars.createIterator() ;
            RooAbsArg* arg ;
            while((arg=(RooAbsArg*)iter->Next())) {
            if (!dvars->find(arg->GetName())) {
            coutE(InputArguments) << "RooHistPoissonGamma::ctor(" << GetName() 
            << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
            assert(0) ;
            }
            }
            delete iter ;
        
           _histList.add(*histPdf) ;
    }
    delete histIter ;
}






//_____________________________________________________________________________
RooHistPoissonGamma::RooHistPoissonGamma(const RooHistPoissonGamma& other, const char* name) :
  RooAbsPdf(other,name), 
  _pdfObsList("pdfObs",this,other._pdfObsList),
  _dataHist(other._dataHist),
  _codeReg(other._codeReg),
  _intOrder(other._intOrder),
  _cdfBoundaries(other._cdfBoundaries),
  _totVolume(other._totVolume),
  _unitNorm(other._unitNorm)
{
  // Copy constructor

  _histObsList.addClone(other._histObsList) ;
  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _pdfObsList.createIterator() ;
}




//_____________________________________________________________________________
RooHistPoissonGamma::~RooHistPoissonGamma()
{
  // Destructor

  delete _histObsIter ;
  delete _pdfObsIter ;
}





//_____________________________________________________________________________
Double_t RooHistPoissonGamma::evaluate() const
{
    // Return the current value: The value of the bin enclosing the current coordinates
    // of the observables, normalized by the histograms contents. Interpolation
    // is applied if the RooHistPoissonGamma is configured to do that

    // Transfer values from   
    if (_pdfObsList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (harg != parg) {
    parg->syncCache() ;
    harg->copyCache(parg,kTRUE) ;
    }
    }
    }


    double prod = 0.;

    for (int k=0; k < _dataHist->numEntries(); k++){
//        _dataHist->get(k);
        const RooArgSet* row = _dataHist->get(k) ;

        Double_t _x =  _dataHist->weight(*row,0,kFALSE,_cdfBoundaries) ;  
        if (_x<0) {
            _x=0 ;
        }

//        cout << k << "\n\nData Value: " << _x << endl;

        int Di = floor(_x);
        vector<double> A;
        vector<double> dA;

        double value = 0;
        _histIter->Reset();
        RooHistPdf* histPdf ;
        const RooArgSet* nset1 = _histList.nset() ;
        
        _muIter->Reset();
        RooAbsReal* thismu;

        
        while((histPdf = (RooHistPdf*)_histIter->Next())) {
            thismu = (RooAbsReal*)_muIter->Next();
            Double_t mu = thismu->getVal();
//            cout << "mu: " << mu << endl;
            RooDataHist* dhist;
            dhist = &histPdf->dataHist();
            const RooArgSet* row = dhist->get(k);
            Double_t weight = mu * dhist->weight(*row, 0 , kFALSE, _cdfBoundaries);  
            Double_t error = fabs( mu * dhist->weightError(RooAbsData::SumW2) );  
            A.push_back(weight);
            dA.push_back(error);
//            cout << "Weight: " << weight << " +- " << lo << " " << hi << endl;
            if (weight<0) {
                weight=0;
            }
        }

        const int MAXSRC = 8; // Maximum number of sources
        const int MAXBUF = 50000; // Maximum count per bin
    
        double s [MAXSRC];
        double f [MAXSRC];
        double x [MAXSRC];
        long double c[MAXSRC][MAXBUF];
    
        int N = A.size(); // Number of sources (N)
    
        // Check inputs
        if ( A.size() != dA.size() )
        {
            std::cout << "**Error - poissongamma - mis-match in number of sources"
            << endl
            << "size(dA): " << dA.size() << " differs from size(A) = " << A.size()
            << std::endl;
            exit(0);
        }
    
    
        // first do zero...      
        for (int j = 0; j < N; ++j)
        {
            s[j] = A[j];
            f[j] = 1.;
            x[j] = 1.;
        
            if ( dA[j] > 0 )
            {
            f[j] = A[j] / (dA[j]*dA[j]);
            } else {
            f[j] = 10000;
            }


            // This 'if' statement makes sure that we don't overwhelm the TMath::Gamma function.
            // In the case that we would overflow this function, we set the scale parameter 
            // to the highest computable scale.
            // This is conservative in the sense that we the prior on the total 
            // count is more smeared out than it would be if we could calculate the 
            // exact value -- it gives you something very close to a poisson with
            // no prior on the mu, but using an exact poisson would be non-conservative
            // as it assumes we know the "true mu" exactly.  This puts a small as possible error on the true mu.  
            // Another option would be to just fail the function and return an error. 

            if(f[j]*s[j] + Di + 0.5 > 140) {
               f[j] = (139.5 - Di)/s[j];
            }

    //        if( s[j] * f[j] > 100) f[j] = 100 / s[j];

            // For bins with zero bin content
            // c[j][0] = 0;
            c[j][0] = 0;
            if ( fabs(f[j]) > 0. )
            {
            x[j] *= f[j];
            s[j] *= f[j];
            }
            else 
            {
            x[j] = 100000.;
            s[j] = 0.;                    
            }
        
            c[j][0] = exp(  (s[j]+0.5) * log(x[j])  - (s[j]+0.5) * log(x[j]+1.) );
        
    //        cout << "TMath::Gamma(" << s[j] + Di + 0.5 << "): " << TMath::Gamma(s[j] + Di + 0.5) << endl;

        }
    
    
        // ...then 1 to D
        if ( Di > 0 )
        {
            for (int k = 1; k < Di+1; ++k)
            for (int j = 0; j < N; ++j)
            {
                c[j][k] = 
                exp(  (s[j]+0.5) * log(x[j])  - (k + s[j]+0.5) * log(x[j]+1) ) *
                TMath::Gamma(s[j] + k + 0.5) / 
                (  TMath::Gamma(k+1) * TMath::Gamma( s[j] + 0.5 )  ) ;
                if(x[j] == 0) c[j][k] = 0;
            }
        }
    
        // compute sum
        double sum = 0.0;
        switch (N)
        {
            case 1:
                sum += c[0][Di];
                break;
    
            case 2:
                for (int j = 0; j < Di+1; ++j){
                sum += 
                    c[0][j] *
                    c[1][Di-j];
                }        
            break;
        
            case 3:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][Di-j-k];
            break;
        
            case 4:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][Di-j-k-l];
            break;
        
            case 5:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                for (int m = 0; m < Di+1-j-k-l; ++m)
                sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][m] *
                    c[4][Di-j-k-l-m];
            break;
        
            case 6:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                for (int m = 0; m < Di+1-j-k-l; ++m)
                for (int n = 0; n < Di+1-j-k-l-m; ++n)
                sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][m] *
                    c[4][n] *
                    c[5][Di-j-k-l-m-n];
            break;
        
            case 7:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                for (int m = 0; m < Di+1-j-k-l; ++m)
                for (int n = 0; n < Di+1-j-k-l-m; ++n)
                for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][m] *
                    c[4][n] *
                    c[5][jj] *
                    c[6][Di-j-k-l-m-n-jj];
            break;
        
            case 8:
                for (int j = 0; j < Di+1; ++j)
                for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                for (int m = 0; m < Di+1-j-k-l; ++m)
                for (int n = 0; n < Di+1-j-k-l-m; ++n)
                for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                for (int kk = 0; kk < Di+1-j-k-l-m-n-jj; ++kk)
                    sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][m] *
                    c[4][n] *
                    c[5][jj] *
                    c[6][kk] *
                    c[7][Di-j-k-l-m-n-jj-kk];
                break;
        };
//        cout << _x << "\t sum: " << sum << "\tprod " << prod << endl;
        prod -= log(sum);
  }

  //cout << "RooHistPoissonGamma::evaluate(" << GetName() << ") ret = " << ret << " at " << ((RooAbsReal*)_histObsList.first())->getVal() << endl ;
  return prod ;
}


//_____________________________________________________________________________
Double_t RooHistPoissonGamma::totVolume() const
{
  // Return the total volume spanned by the observables of the RooHistPoissonGamma

  // Return previously calculated value, if any
  if (_totVolume>0) {
    return _totVolume ;
  }
  _totVolume = 1. ;
  TIterator* iter = _histObsList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    RooRealVar* real = dynamic_cast<RooRealVar*>(arg) ;
    if (real) {
      _totVolume *= (real->getMax()-real->getMin()) ;
    } else {
      RooCategory* cat = dynamic_cast<RooCategory*>(arg) ;
      if (cat) {
	_totVolume *= cat->numTypes() ;
      }
    }
  }
  delete iter ;
  return _totVolume ;
}

namespace {
    bool fullRange(const RooAbsArg& x ,const char* range)  {
      if (range == 0 || strlen(range) == 0 ) return true;
      const RooAbsRealLValue *_x = dynamic_cast<const RooAbsRealLValue*>(&x);
      if (!_x) return false;
      return ( _x->getMin(range) == _x->getMin() && _x->getMax(range) == _x->getMax() ) ; 
    }
}


//_____________________________________________________________________________
Int_t RooHistPoissonGamma::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  // Determine integration scenario. If no interpolation is used,
  // RooHistPoissonGamma can perform all integrals over its dependents
  // analytically via partial or complete summation of the input
  // histogram. If interpolation is used on the integral over
  // all histogram observables is supported

  // First make list of pdf observables to histogram observables
  // and select only those for which the integral is over the full range
  RooArgList hobsl(_histObsList),pobsl(_pdfObsList) ;
  RooArgSet allVarsHist ;
  TIterator* iter = allVars.createIterator() ;
  RooAbsArg* pdfobs ;
  while((pdfobs=(RooAbsArg*)iter->Next())) {
    Int_t idx = pobsl.index(pdfobs) ;
    if (idx>=0) {
      RooAbsArg* hobs = hobsl.at(idx) ;
      if (hobs && fullRange( *hobs, rangeName ) ) {
	allVarsHist.add(*hobs) ;
      }
    }
  }
  delete iter ;

  // Simplest scenario, integrate over all dependents
  RooAbsCollection *allVarsCommon = allVarsHist.selectCommon(_histObsList) ;  
  Bool_t intAllObs = (allVarsCommon->getSize()==_histObsList.getSize()) ;
  delete allVarsCommon ;
  if (intAllObs) {
    analVars.add(allVars) ;
    return 1000 ;
  }

  // Disable partial analytical integrals if interpolation is used
//   if (_intOrder>0) {
//     return 0 ;
//   }

  // Find subset of _histObsList that integration is requested over
  RooArgSet* allVarsSel = (RooArgSet*) allVarsHist.selectCommon(_histObsList) ;
  if (allVarsSel->getSize()==0) {
    delete allVarsSel ;
    return 0 ;
  }

  // Partial integration scenarios.
  // Build unique code from bit mask of integrated variables in depList
  Int_t code(0),n(0) ;
  iter = _histObsList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (allVars.find(arg->GetName())) {
      code |= (1<<n) ;
      analVars.add(*pobsl.at(n)) ;
    }
    n++ ;
  }
  delete iter ;

  return code ;

}



//_____________________________________________________________________________
Double_t RooHistPoissonGamma::analyticalIntegral(Int_t code, const char* /*rangeName*/) const 
{
  // Return integral identified by 'code'. The actual integration
  // is deferred to RooDataHist::sum() which implements partial
  // or complete summation over the histograms contents

  // WVE needs adaptation for rangeName feature
  // Simplest scenario, integration over all dependents
  if (code==1000) {    
    return _dataHist->sum(kFALSE) ;
  }

  // Partial integration scenario, retrieve set of variables, calculate partial sum
  RooArgSet intSet ;
  TIterator* iter = _histObsList.createIterator() ;
  RooAbsArg* arg ;
  Int_t n(0) ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (code & (1<<n)) {
      intSet.add(*arg) ;
    }
    n++ ;
  }
  delete iter ;

  // WVE must sync hist slice list values to pdf slice list
  // Transfer values from   
  if (_pdfObsList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
      parg = (RooAbsArg*)_pdfObsIter->Next() ;
      if (harg != parg) {
	parg->syncCache() ;
	harg->copyCache(parg,kTRUE) ;
      }
    }
  }  


  Double_t ret =  _dataHist->sum(intSet,_histObsList,kTRUE,kTRUE) ;

//   cout << "intSet = " << intSet << endl ;
//   cout << "slice position = " << endl ;
//   _histObsList.Print("v") ;
//   cout << "RooHistPoissonGamma::ai(" << GetName() << ") code = " << code << " ret = " << ret << endl ;

  return ret ;
}



//_____________________________________________________________________________
list<Double_t>* RooHistPoissonGamma::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  // Return sampling hint for making curves of (projections) of this function
  // as the recursive division strategy of RooCurve cannot deal efficiently
  // with the vertical lines that occur in a non-interpolated histogram

  // No hints are required when interpolation is used
  if (_intOrder>0) {
    return 0 ;
  }

  // Check that observable is in dataset, if not no hint is generated
  RooAbsLValue* lvarg = dynamic_cast<RooAbsLValue*>(_dataHist->get()->find(obs.GetName())) ;
  if (!lvarg) {
    return 0 ;
  }

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  list<Double_t>* hint = new list<Double_t> ;

  // Widen range slighty
  xlo = xlo - 0.01*(xhi-xlo) ;
  xhi = xhi + 0.01*(xhi-xlo) ;

  Double_t delta = (xhi-xlo)*1e-8 ;
 
  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]-delta) ;
      hint->push_back(boundaries[i]+delta) ;
    }
  }

  return hint ;
}



//______________________________________________________________________________
std::list<Double_t>* RooHistPoissonGamma::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  // Return sampling hint for making curves of (projections) of this function
  // as the recursive division strategy of RooCurve cannot deal efficiently
  // with the vertical lines that occur in a non-interpolated histogram

  // No hints are required when interpolation is used
  if (_intOrder>0) {
    return 0 ;
  }

  // Check that observable is in dataset, if not no hint is generated
  RooAbsLValue* lvarg = dynamic_cast<RooAbsLValue*>(_dataHist->get()->find(obs.GetName())) ;
  if (!lvarg) {
    return 0 ;
  }

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  list<Double_t>* hint = new list<Double_t> ;

  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]) ;
    }
  }

  return hint ;
}




//_____________________________________________________________________________
Int_t RooHistPoissonGamma::getMaxVal(const RooArgSet& vars) const 
{
  // Only handle case of maximum in all variables
  RooAbsCollection* common = _pdfObsList.selectCommon(vars) ;
  if (common->getSize()==_pdfObsList.getSize()) {
    delete common ;
    return 1;
  }
  delete common ;
  return 0 ;
}


//_____________________________________________________________________________
Double_t RooHistPoissonGamma::maxVal(Int_t code) const 
{
  assert(code==1) ;

  Double_t max(-1) ;
  for (Int_t i=0 ; i<_dataHist->numEntries() ; i++) {
    _dataHist->get(i) ;
    Double_t wgt = _dataHist->weight() ;
    if (wgt>max) max=wgt ;
  }

  return max*1.05 ;
}


//_____________________________________________________________________________
Bool_t RooHistPoissonGamma::importWorkspaceHook(RooWorkspace& ws) 
{  
  // Check if our datahist is already in the workspace
  std::list<RooAbsData*> allData = ws.allData() ;
  std::list<RooAbsData*>::const_iterator iter ;
  for (iter = allData.begin() ; iter != allData.end() ; ++iter) {
    // If your dataset is already in this workspace nothing needs to be done
    if (*iter == _dataHist) {
      return kFALSE ;
    }
  }

  // Check if dataset with given name already exists
  RooAbsData* wsdata = ws.data(_dataHist->GetName()) ;
  if (wsdata) {
    if (wsdata->InheritsFrom(RooDataHist::Class())) {
      // Exists and is of correct type -- adjust internal pointer
      _dataHist = (RooDataHist*) wsdata ;
      return kFALSE ;
    } else {
      // Exists and is NOT of correct type -- abort
      
    }
  }

  // We need to import our datahist into the workspace
  Bool_t flag = ws.import(*_dataHist) ;
  if (flag) {
    coutE(ObjectHandling) << "RooHistPoissonGamma::importWorkspaceHook(" << GetName() 
			  << ") error importing RooDataHist into workspace: dataset of different type with same name already exists." << endl ;
    return kTRUE ;
  }

  // Redirect our internal pointer to the copy in the workspace
  _dataHist = (RooDataHist*) ws.data(_dataHist->GetName()) ;
  return kFALSE ;
}





//_____________________________________________________________________________
void RooHistPoissonGamma::genHist(Int_t seed) const
{
    _gslRan->SetSeed(seed);
    // Generate a fake data histogram using the Poisson-Gamma model 
    if (_pdfObsList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (harg != parg) {
    parg->syncCache() ;
    harg->copyCache(parg,kTRUE) ;
    }
    }
    }


    double prod = 0.;

    for (int k=0; k < _dataHist->numEntries(); k++){
//        _dataHist->get(k);
        const RooArgSet* row = _dataHist->get(k) ;

        Double_t _x =  _dataHist->weight(*row,0,kFALSE,_cdfBoundaries) ;  
        if (_x<0) {
            _x=0 ;
        }


        int Di = floor(_x);
        vector<double> A;
        vector<double> dA;

        double value = 0;
        _histIter->Reset();
        RooHistPdf* histPdf ;
        const RooArgSet* nset1 = _histList.nset() ;
        
        _muIter->Reset();
        RooAbsReal* thismu;

        
        while((histPdf = (RooHistPdf*)_histIter->Next())) {
            thismu = (RooAbsReal*)_muIter->Next();
            Double_t mu = thismu->getVal();
//            cout << "mu: " << mu << endl;
            RooDataHist* dhist;
            dhist = &histPdf->dataHist();
            const RooArgSet* row = dhist->get(k);
            Double_t weight = mu * dhist->weight(*row, 0 , kFALSE, _cdfBoundaries);  
            Double_t error = mu * dhist->weightError(RooAbsData::SumW2) ;  
            A.push_back(weight);
            dA.push_back(error);
//            cout << "Weight: " << weight << " +- " << error << endl;
            if (weight<0) {
                weight=0;
            }
        }


        const int MAXSRC = 8; // Maximum number of sources
        const int MAXBUF = 50000; // Maximum count per bin
    
        double s [MAXSRC];
        double f [MAXSRC];
        double x [MAXSRC];
        long double c[MAXSRC][MAXBUF];
    
        int N = A.size(); // Number of sources (N)
    
        // Check inputs
        if ( A.size() != dA.size() )
        {
            std::cout << "**Error - poissongamma - mis-match in number of sources"
            << endl
            << "size(dA): " << dA.size() << " differs from size(A) = " << A.size()
            << std::endl;
            exit(0);
        }

      double total = 0;


      // first do zero...      
      for (int j = 0; j < N; ++j)
      {
        s[j] = A[j];
        f[j] = 1.;
        x[j] = 1.;

        if ( dA[j] > 0 )
        {
        f[j] = A[j] / (dA[j]*dA[j]);
        } else {
        f[j] = 10000;
        }

        if ( f[j] > 0 )
        {
            x[j] *= f[j];
            s[j] *= f[j];
        }
        total += RooRandom::randomGenerator()->Poisson(_gslRan->Gamma(s[j]+0.5, 1.0/f[j])); 
      } 


      _dataHist->set(*row, total);

   }       
}





//_____________________________________________________________________________
Double_t RooHistPoissonGamma::evaluatePoisson() const
{
    // Return the current value: The value of the bin enclosing the current coordinates
    // of the observables, normalized by the histograms contents. Interpolation
    // is applied if the RooHistPoissonGamma is configured to do that

    // Transfer values from   
    if (_pdfObsList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (harg != parg) {
    parg->syncCache() ;
    harg->copyCache(parg,kTRUE) ;
    }
    }
    }

    double prod = 0.;

    for (int k=0; k < _dataHist->numEntries(); k++){
//        _dataHist->get(k);
        const RooArgSet* row = _dataHist->get(k) ;

        Double_t _x =  _dataHist->weight(*row,0,kFALSE,_cdfBoundaries) ;  
        if (_x<0) {
            _x=0 ;
        }

//        cout << k << "\n\nData Value: " << _x << endl;

        double Di = _x;
        double weight = 0;

        double value = 0;
        _histIter->Reset();
        RooHistPdf* histPdf ;
        const RooArgSet* nset1 = _histList.nset() ;
        
        _muIter->Reset();
        RooAbsReal* thismu;
        
        while((histPdf = (RooHistPdf*)_histIter->Next())) {
            thismu = (RooAbsReal*)_muIter->Next();
            Double_t mu = thismu->getVal();
//            cout << "mu: " << mu << endl;
            RooDataHist* dhist;
            dhist = &histPdf->dataHist();
            const RooArgSet* row = dhist->get(k);
            weight += mu * dhist->weight(*row, 0 , kFALSE, _cdfBoundaries);  

            

            if (weight<0) {
                weight=0;
            }
        }

        prod -= ( Di * log( weight ) - weight - log( TMath::Gamma( Di + 1 ) )  );
    
  }
  //cout << "RooHistPoissonGamma::evaluate(" << GetName() << ") ret = " << ret << " at " << ((RooAbsReal*)_histObsList.first())->getVal() << endl ;
  return prod ;
}






//_____________________________________________________________________________
void RooHistPoissonGamma::genHistAsimov() const
{
    // Generate a fake data histogram using the Poisson-Gamma model 
    if (_pdfObsList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (harg != parg) {
    parg->syncCache() ;
    harg->copyCache(parg,kTRUE) ;
    }
    }
    }


    double prod = 0.;

    for (int k=0; k < _dataHist->numEntries(); k++){
//        _dataHist->get(k);
        const RooArgSet* row = _dataHist->get(k) ;

        Double_t _x =  _dataHist->weight(*row,0,kFALSE,_cdfBoundaries) ;  
        if (_x<0) {
            _x=0 ;
        }


        int Di = floor(_x);
        vector<double> A;
        vector<double> dA;

        double value = 0;
        _histIter->Reset();
        RooHistPdf* histPdf ;
        const RooArgSet* nset1 = _histList.nset() ;
        
        _muIter->Reset();
        RooAbsReal* thismu;

        
        while((histPdf = (RooHistPdf*)_histIter->Next())) {
            thismu = (RooAbsReal*)_muIter->Next();
            Double_t mu = thismu->getVal();
//            cout << "mu: " << mu << endl;
            RooDataHist* dhist;
            dhist = &histPdf->dataHist();
            const RooArgSet* row = dhist->get(k);
            Double_t weight = mu * dhist->weight(*row, 0 , kFALSE, _cdfBoundaries);  
            Double_t error = mu * dhist->weightError(RooAbsData::SumW2) ;  
            A.push_back(weight);
            dA.push_back(error);
//            cout << "Weight: " << weight << " +- " << error << endl;
            if (weight<0) {
                weight=0;
            }
        }



        int N = A.size(); // Number of sources (N)
    

      double total = 0;


      // first do zero...      
      for (int j = 0; j < N; ++j)
      {
        total += A[j]; 
      } 



      _dataHist->set(*row, total);

   }       
}


#ifndef SYSTEMATICS_HPP
#define SYSTEMATICS_HPP

//----------------------------------------------------------------------------
//  File:    Systematics.cc
//  Purpose: Read systematics from file
//  Created: 05-May-2012 Joseph Bochenek
//$Revision: 1.1 $ 
//----------------------------------------------------------------------------

#include <fstream>
#include <string>
#include <vector>
#include <map>
//#include <boost/regex.hpp>
#include <utility>
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"

#ifdef __WITH_CINT__
#include "TObject.h"
#endif


///
class Systematics
{
 public:

	///
	Systematics();
	
	///
	Systematics(int random_seed, std::vector<std::string> vars);

	///
	double doline(std::vector<double> &data, int channel, std::string sample);

	///
	double setvalues(std::vector<double> &data, int channel, std::string sample);

	///
	void setrandom();
	
	~Systematics();
	
	// double getrandom_lnN();
		
	///
	void   close();


	// ---------------------------------------------------
	// Theoretical Uncertainties
	// ---------------------------------------------------
	
	// pdf+alpha_s
	double pdf_gg();  			// mass dep fit
	double pdf_qqbar();  		// mass dep fit

	// PDF 4l acceptance
	double pdf_H4l_accept();  	// For signal 2% mass ind, for ZZ, mass dep. fit 
		
	// QCD Scale uncertainties (given by fits)
	double 	QCDscale_ggH(); //
	double 	QCDscale_ggVV(); //
	double 	QCDscale_qqH(); //
	double 	QCDscale_ttH(); //
	double 	QCDscale_VH(); //
	double 	QCDscale_VV(); //
	double 	scale_ZX(); //

	// Higgs p_T re-weighting uncertainty

	// Branching Ratio Uncertainty (2%)
	double 	BRhiggs_ZZ4l(); //

	// ---------------------------------------------------
	// Instrumental uncertainites
	// ---------------------------------------------------
	
	// Lepton efficiency
	double 	CMS_eff_e(); //
	double 	CMS_eff_m(); //

	// Lepton p_T/E_T scale
	void 	CMS_lepton_scale(std::vector<double> &data); //
	
	// Lepton Trigger efficiency 
	double 	CMS_trigger_e(); //
	double 	CMS_trigger_m(); //


	double 	lumi(); //

	float getCorrectionError(std::vector<float> errCorr, std::vector<float> eff_weight_vec, std::vector<int> muIdx, std::vector<int> eleIdx) const;

	float eleRecoIDEffCorr(const float pt, const float eta, float &corrError, const bool is2011A) const;

	float muRecoIDEffCorr(const float pt, const float eta, float &corrError, const bool is2011A) const;

	float getMCWeight(); //const;

    float leptonRecoIDEffCorr() const;
	// Other systematics



 private:
 
  bool _ok;
  int _nsysts;
  int _nrow;
  TRandom3 _myRandom;

  // Define container for leptons
  struct Lepton {
  	double pt;
  	double eta;
  	double phi;
  	int id; // 11 for electron 13 for muon
  };

  std::vector<Lepton> lepton;
  double Z1Mass;
  double Z2Mass;
  double m4l;
  int _channel;
  std::string _sample;
  
  int _random_seed;

  double _data_to_mc_random_syst_mu;
  double _data_to_mc_random_stat_mu;
  double _data_to_mc_random_reco_e;
  double _data_to_mc_random_id_e;
  double _muon_pt_scale_random;
  double _electron_Et_scale_random;
  
  
  double _pt_scale_random;

  //  typedef double (Systematics::*fptr)();

  //  std::vector<fptr> _f;
  std::ifstream*             _stream;
  std::vector<double>		 _syst_rand;
  std::vector<double>        _data;
  std::vector<double>        _buffer;
  std::vector<std::string>   _name;
  std::map<std::string, int> _var;

  std::string MC_type;

  TH2F *ele_scale_2012;
  TH2F *mu_scale_2012;
  TH2F *mu_scale_2011;

  TFile* ele_scale_factors;
  TFile* mu_scale_factors;
  
#ifdef __WITH_CINT__
 public:
  ClassDef(Systematics, 1)
#endif
};

#endif


//----------------------------------------------------------------------------
//  File:    Systematics.cc
//  Purpose: Read systematics from file
//  Created: 05-May-2012 Joseph Bochenek
//$Revision: 1.3 $ 
//----------------------------------------------------------------------------

#include "Systematics.h"
#include <TLorentzVector.h>
#include <TMath.h>


#ifdef __WITH_CINT__
ClassImp(Systematics)
#endif

#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TGaxis.h>
#include <TFile.h>
#include <TH2F.h>

//#include <RooRandom.h>

//----------------------------------------------------------------------------

  using namespace std;

// Hide error within an anonymous namespace so that
// it is visible only within this programming unit
namespace {
  void error(string message)
  {
    cerr << "read ** error *** " << message << endl;
  }
  
}

Systematics::Systematics()
  : _ok(true),
    _nsysts(15),
    _nrow(0),
    _myRandom(TRandom3())
{}

Systematics::Systematics(int random_seed, std::vector<std::string> vars)
  : _ok(true),
    _nsysts(15),
    _nrow(0),
    _myRandom(TRandom3(random_seed))
{
  cout << "Systematics: Seed = " << random_seed << endl;
  _random_seed = random_seed;
  setrandom();
  _name = vars;
	
  for(int i=0; i < (int)_name.size(); i++)
    {
      _var[_name[i]] = i;
    }


  // Define a pointer to functions that adjust systematics on the weight
  // so that we can loop over these for each event.
  // 	_f.push_back(&Systematics::BRhiggs_ZZ4l);
  // 	_f.push_back(&Systematics::CMS_eff_e);
  // 	_f.push_back(&Systematics::CMS_eff_m);
  // 	_f.push_back(&Systematics::CMS_trigger_e);
  // 	_f.push_back(&Systematics::CMS_trigger_m);
  // 	_f.push_back(&Systematics::lumi);
  // 	_f.push_back(&Systematics::pdf_gg);
  // 	_f.push_back(&Systematics::pdf_H4l_accept);
  // 	_f.push_back(&Systematics::pdf_qqbar);
  // 	_f.push_back(&Systematics::QCDscale_ggH);
  // 	_f.push_back(&Systematics::QCDscale_ggVV);
  // 	_f.push_back(&Systematics::QCDscale_qqH);
  // 	_f.push_back(&Systematics::QCDscale_ttH);
  // 	_f.push_back(&Systematics::QCDscale_VH);

  Lepton lepton0 = *(new Lepton);
  Lepton lepton1 = *(new Lepton);
  Lepton lepton2 = *(new Lepton);
  Lepton lepton3 = *(new Lepton);
  lepton.push_back(lepton0);
  lepton.push_back(lepton1);
  lepton.push_back(lepton2);
  lepton.push_back(lepton3);

  ele_scale_factors = new TFile("Electron_scale_factors_IDISOSIP_combined.root");
  ele_scale_2012 = (TH2F*)gDirectory->Get("h_electron_scale_factor_RECO_ID_ISO_SIP");

  mu_scale_factors = new TFile("muonid_hcp-05.10.2012-with-tk-v2.root");
  mu_scale_2012 = (TH2F*)gDirectory->Get("TH2D_ALL_2012");
  mu_scale_2011 = (TH2F*)gDirectory->Get("TH2D_ALL_2011B");
}


Systematics::~Systematics() { close(); }

// ------------------------------------------------------------------------------------
// Function: doline()
// ------------------------------------------------------------------------------------

// This is the function which varies the systematics nuisance parameters and recalculates
// weights and variables.  The sample must be set according to the sample argument.  The 
// following samples are implemented right now (The sample names are borrowed from the 
// Higgs combination tool):
//
// sig_ggH - gluon gluon fusion signal
// sig_VBF - VBF signal
// bkg_qqzz - qq->ZZ POWHEG MC
// bkg_ggzz - gg->ZZ
// bkg_ZX - reducible background
//
// Data are passed by reference and alterned according to the systematics.
// ------------------------------------------------------------------------------------

double Systematics::doline(std::vector<double> &data, int channel, std::string sample)
{
  _ok = true;
  _nrow++;
  //  cout << "Sample Name: " << sample << endl;
  //  cout << "channel: " << _channel << endl;

  //  	lepton.push_back(new Lepton);

  //  cout << "from  doline: data[0]: " << data.front() << endl;
  //  cout << "from  doline: data[-1]: " << data.back() << endl;

  // Initiate the leptons depending on the channel
  if ( _var.find("lept1_pt") != _var.end() )
    lepton[0].pt = data[_var["lept1_pt"]];
  else
    {
      cout << "error*** lept1_pt not specified in _var" << endl;
      exit(0);
    }
  lepton[0].eta = data[_var["lept1_eta"]];
  lepton[0].phi = data[_var["lept1_phi"]];
  lepton[1].pt = data[_var["lept2_pt"]];
  lepton[1].eta = data[_var["lept2_eta"]];
  lepton[1].phi = data[_var["lept2_phi"]];
  lepton[2].pt = data[_var["lept3_pt"]];
  lepton[2].eta = data[_var["lept3_eta"]];
  lepton[2].phi = data[_var["lept3_phi"]];
  lepton[3].pt = data[_var["lept4_pt"]];
  lepton[3].eta = data[_var["lept4_eta"]];
  lepton[3].phi = data[_var["lept4_phi"]];

  if(channel == 1) { 
    lepton[0].id = 13;
    lepton[1].id = 13;
    lepton[2].id = 13;
    lepton[3].id = 13;
  }
  else if(channel == 2) { 
    lepton[0].id = 11;
    lepton[1].id = 11;
    lepton[2].id = 11;
    lepton[3].id = 11;
  }
  else if(channel == 3) { 
    lepton[0].id = 11;
    lepton[1].id = 11;
    lepton[2].id = 13;
    lepton[3].id = 13;
  }
  else{
    cout << "WARNING: No channel type defined (systematics.cc)" << endl;
    exit(0);
  }
  
  MC_type = "Summer12";
	
  //  cout << "Systematics\n nsyst: " << _nsysts << "\t _f.size(): " << _f.size() << endl;
  double weight_factor = 1.0;
  m4l = data[_var["mbestH"]];
  Z1Mass = data[_var["mZ"]];
  Z2Mass = data[_var["mZstar"]];

  _channel = channel;
  _sample = sample;
  _data = data;
    
  //  cout << "m4l: " << m4l << endl;
  //  cout << "pt3: " << lepton[2].pt << endl;

  weight_factor *= pow(BRhiggs_ZZ4l(), _syst_rand[0]);
  weight_factor *= pow(CMS_eff_e(), _syst_rand[1]);
  weight_factor *= pow(CMS_eff_m(), _syst_rand[2]);
  weight_factor *= pow(CMS_trigger_e(), _syst_rand[3]);
  weight_factor *= pow(CMS_trigger_m(), _syst_rand[4]);
  weight_factor *= pow(lumi(), _syst_rand[5]);
  weight_factor *= pow(pdf_gg(), _syst_rand[6]);
  weight_factor *= pow(pdf_H4l_accept(), _syst_rand[7]);
  weight_factor *= pow(pdf_qqbar(), _syst_rand[8]);
  weight_factor *= pow(QCDscale_ggH(), _syst_rand[9]);
  weight_factor *= pow(QCDscale_ggVV(), _syst_rand[10]);
  weight_factor *= pow(QCDscale_qqH(), _syst_rand[11]);
  weight_factor *= pow(QCDscale_ttH(), _syst_rand[12]);
  weight_factor *= pow(QCDscale_VH(), _syst_rand[13]);
  weight_factor *= pow(scale_ZX(), _syst_rand[14]);
  // To adjust the weight 
  //   for(int i = 0; i < _nsysts; i++){
  // 	weight_factor *= pow(  (this->*_f[i])()  , _syst_rand[i]);	
  //   }
//cout << weight_factor << endl;
  
  // Data/MC correction
  float datamcweight = leptonRecoIDEffCorr();
  weight_factor *= datamcweight;
  //  cout << "data/mc: " <<  datamcweight << endl;
  
  CMS_lepton_scale(data);
  
  return weight_factor;
  
}

void 
Systematics::setrandom()
{
  for(int i = 0; i < _nsysts; i++) _syst_rand.push_back( _myRandom.Gaus(0.0, 1.0) );
  _data_to_mc_random_syst_mu = _myRandom.Gaus(0.0, 1.0);
  _data_to_mc_random_stat_mu = _myRandom.Gaus(0.0, 1.0);
  _data_to_mc_random_reco_e = _myRandom.Gaus(0.0, 1.0);
  _data_to_mc_random_id_e = _myRandom.Gaus(0.0, 1.0);
  _pt_scale_random = _myRandom.Gaus(0.0, 1.0);
  _electron_Et_scale_random = _myRandom.Gaus(0.0, 1.0);
  _muon_pt_scale_random = _myRandom.Gaus(0.0, 1.0);
}


// ------------------------------------------------------------------------------------
// Function: setvalues()
// ------------------------------------------------------------------------------------
//
// Set the values of the variables so that one can vary individual systematics by calling
// their own functions 
//
// sig_ggH - gluon gluon fusion signal
// sig_VBF - VBF signal
// bkg_qqzz - qq->ZZ POWHEG MC
// bkg_ggzz - gg->ZZ
//
// Data are passed by reference and alterned according to the systematics.
// ------------------------------------------------------------------------------------

double Systematics::setvalues(std::vector<double> &data, int channel, std::string sample)
{
    _ok = true;
    _nrow++;

//  cout << "Sample Name: " << sample << endl;
//  cout << "channel: " << _channel << endl;
//  lepton.push_back(new Lepton);

    MC_type = "Summer12";

	// Initiate the leptons depending on the channel
	lepton[0].pt = data[_var["lept1_pt"]];
	lepton[0].eta = data[_var["lept1_eta"]];
	lepton[0].phi = data[_var["lept1_phi"]];
	lepton[1].pt = data[_var["lept2_pt"]];
	lepton[1].eta = data[_var["lept2_eta"]];
	lepton[1].phi = data[_var["lept2_phi"]];
	lepton[2].pt = data[_var["lept3_pt"]];
	lepton[2].eta = data[_var["lept3_eta"]];
	lepton[2].phi = data[_var["lept3_phi"]];
	lepton[3].pt = data[_var["lept4_pt"]];
	lepton[3].eta = data[_var["lept4_eta"]];
	lepton[3].phi = data[_var["lept4_phi"]];

	if(channel == 1) { 
	lepton[0].id = 13;
	lepton[1].id = 13;
	lepton[2].id = 13;
	lepton[3].id = 13;
	}
	else if(channel == 2) { 
	lepton[0].id = 11;
	lepton[1].id = 11;
	lepton[2].id = 11;
	lepton[3].id = 11;
	}
	else if(channel == 3) { 
	lepton[0].id = 11;
	lepton[1].id = 11;
	lepton[2].id = 13;
	lepton[3].id = 13;
	}
	else{
		cout << "WARNING: No channel type defined (systematics.cc)" << endl;
		exit(0);
	}
	
//  cout << "Systematics\n nsyst: " << _nsysts << "\t _f.size(): " << _f.size() << endl;
  double weight_factor = 1.0;
  m4l = data[_var["mbestH"]];
  Z1Mass = data[_var["mZ"]];
  Z2Mass = data[_var["mZstar"]];
  
  _channel = channel;
  _sample = sample;
  _data = data;

  return 0;
}


void
Systematics::close() 
{ 
  ele_scale_factors->Close();
  mu_scale_factors->Close();
}




double Systematics::pdf_gg()  			// mass dep fit
{
  if(_sample == "bkg_ggzz") return (1 + 0.0035*sqrt(m4l - 30));
  // For now linearly interpolate in the range 110-160
  if(_sample == "sig_ggH") return ( 1.075 + (0.0025)*(m4l-110)/50 );
  return 1.0; 
}

double Systematics::pdf_qqbar()  		// mass dep fit
{
  if(_sample == "bkg_qqzz") return (1 + 0.0066*sqrt(m4l - 10))	;
  if(_sample == "sig_VBF") return (1.022)	;
  return 1.0;
}

// PDF 4l acceptance
double Systematics::pdf_H4l_accept()  	// For signal 2% mass ind, for ZZ, mass dep. fit 
{
  if(_sample == "sig_ggH") return 1.02;
  if(_sample == "sig_VBF") return 1.02;

  return 1.0;
}

// QCD Scale uncertainties (given by fits)
double 	Systematics::QCDscale_ggH() //
{
  // Again, linearly interpolate in the range 110-160
  if(_sample == "sig_ggH") return ( 1.1035 + (1.0905-1.1035)*(m4l-110)/50 );
  return 1.00;
}

// These are left as 1.00 for the following systematics since we only use the qq signal anyways
double 	Systematics::QCDscale_ggVV() //
{
  if(_sample == "bkg_qqzz")return 1.04 + 0.10*sqrt((m4l +  40)/40);
  return 1.00;
}

double 	Systematics::QCDscale_qqH() //
{
  return 1.00;
}

double 	Systematics::QCDscale_ttH() //
{
  return 1.00;
}

double 	Systematics::QCDscale_VH() //
{
  return 1.00;
}

double 	Systematics::QCDscale_VV() //
{
  if(_sample == "bkg_qqzz")return 1.00 + 0.01*sqrt((m4l - 20)/13);
  return 1.00;
}

// Higgs p_T re-weighting uncertainty

// Branching Ratio Uncertainty (2%)
double 	Systematics::BRhiggs_ZZ4l() //
{
  if(_sample == "sig_ggH") return 1.02;
  if(_sample == "sig_VBF") return 1.02;

  return 1.0;
}

// ---------------------------------------------------
// Instrumental uncertainites
// ---------------------------------------------------

// Lepton efficiency
double 	Systematics::CMS_eff_e() //
{
  if(_channel == 2) return 1.0447;
  if(_channel == 3) return 1.0224;
  else return 1.00;
}

double 	Systematics::CMS_eff_m() //
{
  if(_channel == 1) return 1.0224;
  if(_channel == 3) return 1.0224;
  else return 1.00;
}



// Lepton Trigger efficiency 
double 	Systematics::CMS_trigger_e() //
{
  if(_sample == "sig_ZX") return 1.00;
  return 1.015;
}

double 	Systematics::CMS_trigger_m() //
{
  if(_sample == "sig_ZX") return 1.00;
  return 1.015;
}

double 	Systematics::lumi() //
{
  if(_sample == "sig_ZX") return 1.00;
  return 1.045;
}

double 	Systematics::scale_ZX() //
{
  if(_sample == "sig_ZX") return 1.5;
  return 1.00;
}

// Lepton p_T/E_T scale
void 	Systematics::CMS_lepton_scale(std::vector<double> &data) //
{
  if(_sample == "sig_ZX") exit(0);

  const double _muon_pt_scale_sigma = 1.005;
  const double _electron_Et_scale_sigma = 1.004;
  double e_scale_factor = pow(  _electron_Et_scale_sigma, _electron_Et_scale_random );
  double mu_scale_factor = pow(  _muon_pt_scale_sigma, _muon_pt_scale_random );
	
  //	cout << "electron scale " << e_scale_factor << endl;
  //	cout << "muon scale " << mu_scale_factor << endl;
	
  TLorentzVector *m = new TLorentzVector[4];
  for(int nLep=0;nLep<4;nLep++){
    m[nLep].SetPtEtaPhiM(lepton[nLep].pt,lepton[nLep].eta,lepton[nLep].phi,0.1506583);
    //	  cout << "lepton " << nLep << " pt " << m[nLep].Pt() << endl;
    if(lepton[nLep].id == 11){ 
      //		if(fabs(lepton[nLep].eta) < 1.5) m[nLep].SetPerp(lepton[nLep].pt-lepton[nLep].pt*0.005);
      //		if(fabs(lepton[nLep].eta) > 1.5) m[nLep].SetPerp(lepton[nLep].pt-lepton[nLep].pt*0.02);
    }
  }

  TLorentzVector new_m[4];
  for(int nLep=0;nLep<4;nLep++){
    if(lepton[nLep].id == 11){ 
      new_m[nLep].SetPtEtaPhiM(  (e_scale_factor*m[nLep].Pt()),lepton[nLep].eta,lepton[nLep].phi,0.1506583);
      //		if(fabs(lepton[nLep].eta) < 1.5) new_m[nLep].SetPerp(new_m[nLep].Pt()-new_m[nLep].Pt()*0.005);
      //		if(fabs(lepton[nLep].eta) > 1.5) new_m[nLep].SetPerp(new_m[nLep].Pt()-new_m[nLep].Pt()*0.02);
    }
    else if(lepton[nLep].id == 13){ 
      new_m[nLep].SetPtEtaPhiM(  (mu_scale_factor * m[nLep].Pt()),lepton[nLep].eta,lepton[nLep].phi,0.1506583);
    }
    //	  data[_var["mZ"]] = new_m[nLep].Pt();
  }

  //  bool isZ1first = false;
  TLorentzVector m12 = m[0] + m[1];
  TLorentzVector m13 = m[0] + m[2];
  TLorentzVector m14 = m[0] + m[3];
  TLorentzVector m23 = m[1] + m[2];
  TLorentzVector m24 = m[1] + m[3];
  TLorentzVector m34 = m[2] + m[3];

  delete [] m;

  TLorentzVector new_m12 = new_m[0] + new_m[1];
  TLorentzVector new_m13 = new_m[0] + new_m[2];
  TLorentzVector new_m14 = new_m[0] + new_m[3];
  TLorentzVector new_m23 = new_m[1] + new_m[2];
  TLorentzVector new_m24 = new_m[1] + new_m[3];
  TLorentzVector new_m34 = new_m[2] + new_m[3];
    
  double new_Z1_Mass = 0;
  double new_Z2_Mass = 0;

  TLorentzVector new_m4l;
  new_m4l = new_m[0] + new_m[1] + new_m[2] + new_m[3];

  data[_var["mbestH"]]  = new_m4l.M();
  data[_var["lept1_pt"]] = new_m[0].Pt();
  data[_var["lept2_pt"]] = new_m[1].Pt();
  data[_var["lept3_pt"]] = new_m[2].Pt();
  data[_var["lept4_pt"]] = new_m[3].Pt();
    
  if(fabs(m12.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m12.M();
    new_Z2_Mass = new_m34.M();	  
  }
  else if(fabs(m13.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m13.M();
    new_Z2_Mass = new_m24.M();
  }
  else if(fabs(m14.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m14.M();
    new_Z2_Mass = new_m23.M();
  }
  else if(fabs(m23.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m23.M();
    new_Z2_Mass = new_m14.M();
  }
  else if(fabs(m24.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m24.M();
    new_Z2_Mass = new_m13.M();
  }
  else if(fabs(m34.M()-Z1Mass) < 0.01){
    new_Z1_Mass = new_m34.M();
    new_Z2_Mass = new_m12.M();
  } else {
    // cout << "NO Z/Z* MATCH (FSR photon?)" << endl;
    new_Z1_Mass = Z1Mass;
    new_Z2_Mass = Z2Mass;
    data[_var["mbestH"]]  = m4l;
    // assert(0);
  }
    
  data[_var["mZ"]]      = new_Z1_Mass;
  data[_var["mZstar"]]  = new_Z2_Mass;

  //if(fabs(new_Z1_Mass - Z1Mass) > 5 ) {assert(0); cout << "NOT MATCHED" ;}
  //  if(fabs(new_Z1_Mass - Z1Mass) > 5 ) {cout << "NOT MATCHED" << endl;}

  
    //    cout << "old Z1: " << Z1Mass << endl;
    //    cout << "new Z1: " << new_Z1_Mass << endl;
    //    cout << "difference: " << (new_Z1_Mass - Z1Mass) << endl;

    //    cout << "old Z2: " << Z2Mass << endl;
    //    cout << "new Z2: " << new_Z2_Mass << endl;

    //    cout << "old m4l: " << m4l << endl;
    //    cout << "new m4l: " << new_m4l.M() << endl;    
  

}





float Systematics::getMCWeight() //const
{

  bool is2011A = false;
  float MC_weight = 1.;
  float error_corr = 0.;
  float eff_weight = 1.;
  vector<Float_t> eff_weight_vec, errCorr;
  vector<int> muIdx, eleIdx;

  float RunAFraction = 0.465;

  //Double_t whatPeriod = RooRandom::randomGenerator()->Uniform();
  Double_t whatPeriod = _myRandom.Rndm();
  if(whatPeriod < RunAFraction){
    is2011A =true;
  }
  else{
    is2011A = false;
  }


  for(int nLep=0;nLep<4;nLep++){

    eff_weight_vec.push_back(1.);
    errCorr.push_back(0.);

    //	cout << "pt " << lepton[nLep].pt << endl;
    //	cout << "eta " << fabs(lepton[nLep].eta) << endl;
	
    //if(isZ1first && nLep < 2) continue;
    //if(!isZ1first && nLep > 1) continue;	
    if(fabs(lepton[nLep].id) == 13){
      eff_weight *= muRecoIDEffCorr( lepton[nLep].pt,fabs(lepton[nLep].eta),errCorr.at(nLep),is2011A); // * muIsolCorr( lepton[nLep].pt,lepton[nLep].eta,errCorr.at(nLep));
      eff_weight_vec.at(nLep) = muRecoIDEffCorr( lepton[nLep].pt,fabs(lepton[nLep].eta),errCorr.at(nLep),is2011A);
      muIdx.push_back(nLep);
    } else if(fabs(lepton[nLep].id) == 11){ 
      eff_weight *= eleRecoIDEffCorr(lepton[nLep].pt,fabs(lepton[nLep].eta),errCorr.at(nLep),is2011A); // * eleIsolCorr(lepton[nLep].pt,lepton[nLep].eta,errCorr.at(nLep));
      eff_weight_vec.at(nLep) = eleRecoIDEffCorr(lepton[nLep].pt,fabs(lepton[nLep].eta),errCorr.at(nLep),is2011A);   
      eleIdx.push_back(nLep);
    }

    //    cout << "Lepton " << nLep << " efficiency weight: " << eff_weight << endl;
  }

  error_corr = getCorrectionError(errCorr,eff_weight_vec,muIdx,eleIdx);
  //  cout << "Correction: " << error_corr << endl;

  MC_weight = MC_weight*eff_weight + error_corr;

  return MC_weight;
}




float Systematics::leptonRecoIDEffCorr() const
{
        if(_sample == "sig_ZX") return 1.00;

        /////////////Lepton Efficiency Scale Factrons/////////////
        // Load histograms
        //////////////////////////////////////////////////////////

        // Execute Efficiency Reweighting
        Double_t eff_weight = 1.;
        Double_t eff_weight_corr = 1.;

        for(int i = 0; i < 4; ++i){
            Int_t PDGID  = lepton[i].id ;
            Double_t Pt  = lepton[i].pt ;
            Double_t Eta = lepton[i].eta;

            if(fabs(PDGID) == 11)
            {
            if ( (Pt > 7 && Pt < 200) && (fabs(Eta) < 2.5) ) {
            Int_t binx = ele_scale_2012->GetXaxis()->FindBin(Pt);
            Int_t biny = ele_scale_2012->GetYaxis()->FindBin(Eta);
            eff_weight*= ele_scale_2012->GetBinContent(binx,biny);
            eff_weight_corr *= (ele_scale_2012->GetBinContent(binx,biny) + ele_scale_2012->GetBinError(binx,biny) * _data_to_mc_random_id_e);
//            cout << "\tweight: " << ele_scale_2012->GetBinContent(binx,biny) << "\tweight corr (e): " << (ele_scale_2012->GetBinContent(binx,biny) + ele_scale_2012->GetBinError(binx,biny) * _data_to_mc_random_id_e) << "\t rand: " << _data_to_mc_random_id_e << endl;
//            cout << "elec " << i << ", pt = " << Pt << ", Eta= " << Eta << " Data/MC Scale: " << eff_weight << " Corrected: " << eff_weight_corr << endl;
            } else {
//            cout << "electron (pt = " << Pt << ", eta = " << Eta << " ) out of range" << endl;
            }
            }
            else if(fabs(PDGID) == 13)
            {
            if( MC_type == "Fall11")
            {
            if ( (Pt > 5 && Pt < 100) && (fabs(Eta) < 2.4) ) {
            Int_t binx = mu_scale_2011->GetXaxis()->FindBin(Pt);
            Int_t biny = mu_scale_2011->GetYaxis()->FindBin(Eta);
            eff_weight*=mu_scale_2011->GetBinContent(binx,biny);
            eff_weight_corr *= (mu_scale_2011->GetBinContent(binx,biny) + mu_scale_2011->GetBinError(binx,biny) * _data_to_mc_random_syst_mu);
//            cout << "\tweight (mu): " << _data_to_mc_random_syst_mu << endl;
            } else {
	      //            cout << "muon (pt = " << Pt << ", eta = " << Eta << " ) out of range" << endl;
            }
            }
            else if( MC_type == "Summer12")
            {
            if ( (Pt > 5 && Pt < 100) && (fabs(Eta) < 2.4) ) {
            Int_t binx = mu_scale_2012->GetXaxis()->FindBin(Pt);
            Int_t biny = mu_scale_2012->GetYaxis()->FindBin(Eta);
            eff_weight *= mu_scale_2012->GetBinContent(binx,biny);
            eff_weight_corr *= (mu_scale_2012->GetBinContent(binx,biny) + mu_scale_2012->GetBinError(binx,biny) * _data_to_mc_random_syst_mu);
//            cout << "\tweight: " << mu_scale_2012->GetBinContent(binx,biny) <<  "\tweight corr (mu): " << (mu_scale_2012->GetBinContent(binx,biny) + mu_scale_2012->GetBinError(binx,biny) * _data_to_mc_random_syst_mu) << "\t rand: " << _data_to_mc_random_syst_mu << endl;
            } else {
	      //            cout << "muon (pt = " << Pt << ", eta = " << Eta << " ) out of range" << endl;
            }
            }
            }
        }
        
//        cout << "old weight: " << eff_weight << "\t new weight: " << eff_weight_corr << endl; 
        
        double datamc_weight_correction_weight = eff_weight_corr / eff_weight; 
//        cout << "Data/MC: " << datamc_weight_correction_weight  <<  endl; 
        return datamc_weight_correction_weight;

}





// Data/MC eff correction for muons
float Systematics::muRecoIDEffCorr(const float pt, const float eta, float &corrError, const bool is2011A) const
{
  float mc = 0.;
  float data = 0.;
  float datamc = 0.;
  float sigmap = 0.;
  float sigmam = 0.;
  float syst = 0.;

  if (is2011A) { // Period 1-4
    if      (pt>5 &&pt<=7 && eta<1.2)          {mc= 0.980421; data= 0.976988; datamc= 0.996499; sigmap= 0.00713744; sigmam= 0.00712514;  syst= 0.010;}
    else if (pt>5 &&pt<=7 && eta>=1.2&&eta<2.5){mc= 0.978089; data= 0.995969; datamc= 1.01828;  sigmap= 0.00426506; sigmam= 0.0132881;  syst= 0.014;}
    else if (pt>7 &&pt<=20&& eta<1.2)          {mc= 0.985403; data= 0.976819; datamc= 0.991289; sigmap= 0.00247849; sigmam= 0.00250117;  syst= 0.010;}
    else if (pt>7 &&pt<=20&& eta>=1.2&&eta<2.5){mc= 0.981673; data= 0.979046; datamc= 0.997324; sigmap= 0.005072;   sigmam= 0.00513563;  syst= 0.014;}
    else if (pt>20&&pt<=30&& eta<0.9)          {mc= 0.979682; data= 0.976321; datamc= 0.99657;  sigmap= 0.00192501; sigmam= 0.00193656;  syst= 0.002;}
    else if (pt>20&&pt<=30&& eta>=0.9&&eta<1.2){mc= 0.988456; data= 0.988986; datamc= 1.00054;  sigmap= 0.00324797; sigmam= 0.00329139;  syst= 0.002;}
    else if (pt>20&&pt<=30&& eta>=1.2&&eta<1.6){mc= 0.990458; data= 0.989029; datamc= 0.998557; sigmap= 0.00272166; sigmam= 0.0027535;  syst= 0.004;}
    else if (pt>20&&pt<=30&& eta>=1.6&&eta<2.1){mc= 0.987618; data= 0.9859;   datamc= 0.998261; sigmap= 0.00246911; sigmam= 0.0024959;  syst= 0.004;}
    else if (pt>20&&pt<=30&& eta>=2.1&&eta<2.5){mc= 0.988993; data= 0.986328; datamc= 0.997305; sigmap= 0.00424767; sigmam= 0.00427749;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta<0.9)          {mc= 0.986718; data= 0.983071; datamc= 0.996304; sigmap= 0.000355638;sigmam= 0.000357272;  syst= 0.002;}
    else if (pt>30&&pt<=50&& eta>=0.9&&eta<1.2){mc= 0.992839; data= 0.989888; datamc= 0.997028; sigmap= 0.00056499; sigmam= 0.000572236;  syst= 0.002;}
    else if (pt>30&&pt<=50&& eta>=1.2&&eta<1.6){mc= 0.992562; data= 0.989956; datamc= 0.997374; sigmap= 0.000542063;sigmam= 0.000548145;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta>=1.6&&eta<2.1){mc= 0.981146; data= 0.976587; datamc= 0.995353; sigmap= 0.000724634;sigmam= 0.000730356;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta>=2.1&&eta<2.5){mc= 0.986463; data= 0.981917; datamc= 0.995391; sigmap= 0.00128694; sigmam= 0.00129656;  syst= 0.004;}
    else if (pt>50        && eta<0.9)          {mc= 0.986101; data= 0.981975; datamc= 0.995817; sigmap= 0.00110161; sigmam= 0.00111603;  syst= 0.002;}
    else if (pt>50        && eta>=0.9&&eta<1.2){mc= 0.992669; data= 0.991975; datamc= 0.999301; sigmap= 0.00179876; sigmam= 0.00184619;  syst= 0.002;}
    else if (pt>50        && eta>=1.2&&eta<1.6){mc= 0.991115; data= 0.987437; datamc= 0.996288; sigmap= 0.00184449; sigmam= 0.00188304;  syst= 0.004;}
    else if (pt>50        && eta>=1.6&&eta<2.1){mc= 0.975529; data= 0.967573; datamc= 0.991844; sigmap= 0.00262898; sigmam= 0.00266344;  syst= 0.004;}
    else if (pt>50        && eta>=2.1&&eta<2.5){mc= 0.976459; data= 0.96634;  datamc= 0.989637; sigmap= 0.00600008; sigmam= 0.00602842;  syst= 0.004;}
    else assert(0);

  } else { //Period 5 (2011B)
    
    if      (pt>5 &&pt<=7 && eta<1.2)          {mc= 0.980421; data= 0.960302; datamc= 0.97948;  sigmap= 0.0211709;  sigmam= 0.0208153;  syst= 0.010;}
    else if (pt>5 &&pt<=7 && eta>=1.2&&eta<2.5){mc= 0.978089; data= 0.88114;  datamc= 0.90088;  sigmap= 0.0346231;  sigmam= 0.0333162;  syst= 0.014;}
    else if (pt>7 &&pt<=20&& eta<1.2)          {mc= 0.985403; data= 0.978096; datamc= 0.992585; sigmap= 0.00893075; sigmam= 0.0089176;  syst= 0.010;}
    else if (pt>7 &&pt<=20&& eta>=1.2&&eta<2.5){mc= 0.981673; data= 0.938443; datamc= 0.955963; sigmap= 0.0190749;  sigmam= 0.0187154;  syst= 0.014;}
    else if (pt>20&&pt<=30&& eta<0.9)          {mc= 0.979682; data= 0.979173; datamc= 0.999481; sigmap= 0.00326466; sigmam= 0.00329495;  syst= 0.002;}
    else if (pt>20&&pt<=30&& eta>=0.9&&eta<1.2){mc= 0.988456; data= 0.984482; datamc= 0.99598;  sigmap= 0.00553251; sigmam= 0.00563272;  syst= 0.002;}
    else if (pt>20&&pt<=30&& eta>=1.2&&eta<1.6){mc= 0.990458; data= 0.9539;   datamc= 0.963091; sigmap= 0.00545461; sigmam= 0.00552482;  syst= 0.004;}
    else if (pt>20&&pt<=30&& eta>=1.6&&eta<2.1){mc= 0.987618; data= 0.957277; datamc= 0.969279; sigmap= 0.0047797;  sigmam= 0.00483757;  syst= 0.004;}
    else if (pt>20&&pt<=30&& eta>=2.1&&eta<2.5){mc= 0.988993; data= 0.96005;  datamc= 0.970735; sigmap= 0.00895747; sigmam= 0.00896605;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta<0.9)          {mc= 0.986718; data= 0.982805; datamc= 0.996034; sigmap= 0.000583258;sigmam= 0.000588616;  syst= 0.002;}
    else if (pt>30&&pt<=50&& eta>=0.9&&eta<1.2){mc= 0.992839; data= 0.984633; datamc= 0.991735; sigmap= 0.00103789; sigmam= 0.00105713;  syst= 0.002;}
    else if (pt>30&&pt<=50&& eta>=1.2&&eta<1.6){mc= 0.992562; data= 0.960313; datamc= 0.96751;  sigmap= 0.00131662; sigmam= 0.00133161;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta>=1.6&&eta<2.1){mc= 0.981146; data= 0.943761; datamc= 0.961896; sigmap= 0.00155456; sigmam= 0.00156864;  syst= 0.004;}
    else if (pt>30&&pt<=50&& eta>=2.1&&eta<2.5){mc= 0.986463; data= 0.958717; datamc= 0.971873; sigmap= 0.00276439; sigmam= 0.00278628;  syst= 0.004;}
    else if (pt>50        && eta<0.9)          {mc= 0.986101; data= 0.978958; datamc= 0.992756; sigmap= 0.0020995;  sigmam= 0.00213748;  syst= 0.002;}
    else if (pt>50        && eta>=0.9&&eta<1.2){mc= 0.992669; data= 0.988582; datamc= 0.995882; sigmap= 0.00350208; sigmam= 0.0036181;  syst= 0.002;}
    else if (pt>50        && eta>=1.2&&eta<1.6){mc= 0.991115; data= 0.955168; datamc= 0.96373;  sigmap= 0.00456411; sigmam= 0.00464968;  syst= 0.004;}
    else if (pt>50        && eta>=1.6&&eta<2.1){mc= 0.975529; data= 0.931463; datamc= 0.954828; sigmap= 0.00559813; sigmam= 0.00567451;  syst= 0.004;}
    else if (pt>50        && eta>=2.1&&eta<2.5){mc= 0.976459; data= 0.92791;  datamc= 0.95028;  sigmap= 0.0131293;  sigmam= 0.0130737;  syst= 0.004;}
    else assert(0);   
  }

  // Vary systematics according to Normal Distribution;
  syst = _data_to_mc_random_syst_mu * syst;
  
  corrError = max(sigmap, sigmam);
  
  corrError = _data_to_mc_random_stat_mu * corrError;
  
  corrError = sqrt(corrError*corrError + syst*syst);

  return datamc;
}

float Systematics::eleRecoIDEffCorr(const float pt, const float eta, float &corrError, const bool is2011A) const
{

  float reco_mc = 0.;
  float reco_data = 0.;
  float reco_datamc = 0.;
  float reco_sigma = 0.;  

  // RECO  
  if      (pt<15)          {reco_data= 0.9220;  reco_mc= 0.8781;  reco_datamc= 1.049; reco_sigma= 0.027;}
  else if (pt>15&&pt<=20)  {reco_data= 0.9430;  reco_mc= 0.9140;  reco_datamc= 1.032; reco_sigma= 0.028;}
  else if (pt>20&&pt<=30)  {reco_data= 0.9337;  reco_mc= 0.9487;  reco_datamc= 0.984; reco_sigma= 0.008;}
  else if (pt>30&&pt<=40)  {reco_data= 0.9625;  reco_mc= 0.9650;  reco_datamc= 0.997; reco_sigma= 0.008;}
  else if (pt>40&&pt<=50)  {reco_data= 0.9732;  reco_mc= 0.9723;  reco_datamc= 1.001; reco_sigma= 0.002;}
  else if (pt>50)          {reco_data= 0.9796;  reco_mc= 0.9783;  reco_datamc= 1.001; reco_sigma= 0.003;}


  float mc = 0.;
  float data = 0.;
  float datamc = 0.;
  float sigma = 0.;

  // ID
  if      (pt>=7 &&pt<=10&&  eta<1.479)             {data= 0.733; mc= 0.736; datamc= 0.995; sigma= 0.022;} 
  else if (pt>=7 &&pt<=10&&  eta>=1.479&&eta<2.5)   {data= 0.674; mc= 0.66;  datamc= 1.02; sigma= 0.03;}
  else if (pt>10        &&  eta>=1.4442&&eta<1.566){data= 0.977; mc= 0.96;  datamc= 1.017; sigma= 0.01;}
  else if (pt>10 &&pt<=15&& eta<0.78)              {data= 0.871; mc= 0.863; datamc= 1.01; sigma= 0.01;}  
  else if (pt>10 &&pt<=15&& eta>=0.78&&eta<1.4442) {data= 0.912; mc= 0.902; datamc= 1.01; sigma= 0.007;}  
  else if (pt>10 &&pt<=15&& eta>=1.566&&eta<2.0)   {data= 0.739; mc= 0.708; datamc= 1.04; sigma= 0.01;}  
  else if (pt>10 &&pt<=15&& eta>=2.0&&eta<2.5)     {data= 0.858; mc= 0.786; datamc= 1.09; sigma= 0.01;}    
  else if (pt>15 &&pt<=20&& eta<0.78)              {data= 0.929; mc= 0.932; datamc= 0.997; sigma= 0.004;} 
  else if (pt>15 &&pt<=20&& eta>=0.78&&eta<1.4442) {data= 0.935; mc= 0.946; datamc= 0.988; sigma= 0.0041;} 
  else if (pt>15 &&pt<=20&& eta>=1.566&&eta<2.0)   {data= 0.871; mc= 0.841; datamc= 1.036; sigma= 0.006;} 
  else if (pt>15 &&pt<=20&& eta>=2.0&&eta<2.5)     {data= 0.924; mc= 0.899; datamc= 1.028; sigma= 0.005;} 
  else if (pt>20 &&pt<=30&& eta<0.78)              {data= 0.947; mc= 0.959; datamc= 0.987; sigma= 0.001;} 
  else if (pt>20 &&pt<=30&& eta>=0.78&&eta<1.4442) {data= 0.97;  mc= 0.961; datamc= 1.009; sigma= 0.001;}
  else if (pt>20 &&pt<=30&& eta>=1.566&&eta<2.0)   {data= 0.944; mc= 0.932; datamc= 1.012; sigma= 0.002;} 
  else if (pt>20 &&pt<=30&& eta>=2.0&&eta<2.5)     {data= 0.959; mc= 0.952; datamc= 1.007; sigma= 0.002;} 
  else if (pt>30 &&pt<=40&& eta<0.78)              {data= 0.964; mc= 0.969; datamc= 0.995; sigma= 0.001;}
  else if (pt>30 &&pt<=40&& eta>=0.78&&eta<1.4442) {data= 0.97;  mc= 0.973; datamc= 0.997; sigma= 0.001;} 
  else if (pt>30 &&pt<=40&& eta>=1.566&&eta<2.0)   {data= 0.971; mc= 0.962; datamc= 1.009; sigma= 0.001;} 
  else if (pt>30 &&pt<=40&& eta>=2.0&&eta<2.5)     {data= 0.969; mc= 0.966; datamc= 1.004; sigma= 0.001;} 
  else if (pt>40 &&pt<=50&& eta<0.78)              {data= 0.971; mc= 0.976; datamc= 0.994; sigma= 0.001;} 
  else if (pt>40 &&pt<=50&& eta>=0.78&&eta<1.4442) {data= 0.979; mc= 0.982; datamc= 0.997; sigma= 0.001;} 
  else if (pt>40 &&pt<=50&& eta>=1.566&&eta<2.0)   {data= 0.981; mc= 0.978; datamc= 1.003; sigma= 0.001;}
  else if (pt>40 &&pt<=50&& eta>=2.0&&eta<2.5)     {data= 0.978; mc= 0.976; datamc= 1.00241; sigma= 0.001;} 
  else if (pt>50 &&pt<=60&& eta<0.78)              {data= 0.974; mc= 0.98;  datamc= 0.994; sigma= 0.001;} 
  else if (pt>50 &&pt<=60&& eta>=0.78&&eta<1.4442) {data= 0.984; mc= 0.985; datamc= 0.996; sigma= 0.001;} 
  else if (pt>50 &&pt<=60&& eta>=1.566&&eta<2.0)   {data= 0.988; mc= 0.982; datamc= 1.006; sigma= 0.001;}
  else if (pt>50 &&pt<=60&& eta>=2.0&&eta<2.5)     {data= 0.982; mc= 0.979; datamc= 1.002; sigma= 0.002;} 
  else if (pt>60        &&  eta<0.78)              {data= 0.979; mc= 0.982; datamc= 0.997; sigma= 0.001;} 
  else if (pt>60        &&  eta>=0.78&&eta<1.4442) {data= 0.985; mc= 0.987; datamc= 0.998; sigma= 0.001;}
  else if (pt>60        &&  eta>=1.566&&eta<2.0)   {data= 0.991; mc= 0.985; datamc= 1.006; sigma= 0.001;} 
  else if (pt>60        &&  eta>=2.0&&eta<2.5)     {data= 0.983; mc= 0.984; datamc= 0.999; sigma= 0.003;} 
  //  else assert(0);
  else sigma = 0;

  // Vary systematics according to Normal Distribution;

  reco_sigma = _data_to_mc_random_reco_e * reco_sigma;

  sigma = _data_to_mc_random_id_e * sigma;

  corrError = sqrt(reco_sigma*reco_sigma+sigma*sigma);

  return reco_datamc*datamc;
}


float Systematics::getCorrectionError(vector<float> errCorr, vector<float> eff_weight_vec, vector<int> muIdx, vector<int> eleIdx) const
{

  Float_t error_corr = 0.;

  // Calculate the error on the 4l correction factor by propagation of uncertainty
  error_corr = TMath::Sqrt( pow( errCorr.at(0)*eff_weight_vec.at(1)*eff_weight_vec.at(2)*eff_weight_vec.at(3) ,2) +
			    pow( errCorr.at(1)*eff_weight_vec.at(2)*eff_weight_vec.at(3)*eff_weight_vec.at(0) ,2) +
			    pow( errCorr.at(2)*eff_weight_vec.at(3)*eff_weight_vec.at(0)*eff_weight_vec.at(1) ,2) +
			    pow( errCorr.at(3)*eff_weight_vec.at(0)*eff_weight_vec.at(1)*eff_weight_vec.at(2) ,2) );

  return error_corr;
}

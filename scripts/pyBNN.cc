//----------------------------------------------------------------------------
//  File:    BNN.cc
//  Purpose: Calculate BNN value 
//  Created: 16-Dec-2012 Joseph Bochenek
//$Revision: 1.1.1.2 $ 
//----------------------------------------------------------------------------


#include <TLorentzVector.h>
#include <TMath.h>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TGaxis.h>
#include <TRandom.h>
//#include <RooRandom.h>



#include "classes/bnn_4e_8TeV_new.cpp"
#include "classes/bnn_4mu_8TeV_new.cpp"
#include "classes/bnn_2e2mu_8TeV_new.cpp"

#include "classes/bnn_4e_7TeV_new.cpp"
#include "classes/bnn_4mu_7TeV_new.cpp"
#include "classes/bnn_2e2mu_7TeV_new.cpp"

#include "classes/Graviton_2e2mu.cpp"
#include "classes/Graviton_4mu.cpp"
#include "classes/Graviton_4e.cpp"

#include "classes/Higgs0M_2e2mu.cpp"
#include "classes/Higgs0M_4e.cpp"
#include "classes/Higgs0M_4mu.cpp"


#include "classes/bnn_4mu_8TeV_new_pt.cpp"
#include "classes/bnn_4e_8TeV_new_pt.cpp"
#include "classes/bnn_2e2mu_8TeV_new_pt.cpp"

//----------------------------------------------------------------------------



class BayesianNN{
    public: 
        BayesianNN() {}

        double getBNNvalue(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ, 
            float ZZpT = -1., 
            float ZZY = -1.
        );

        double getBNNvalueSPS(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ 
        );

        double getBNNvalueGrav(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ 
        );

        void setChannel(std::string channel) { _channel = channel; };
        void setEra(std::string era) { (era=="7TeV") ? _is7TeV=true : _is7TeV=false; }

        void getEra() {  std::cout << _is7TeV << std::endl; };
        void getChannel() { std::cout << _channel << std::endl; }


    private:
        bool _is7TeV;
        std::string _channel;
};



double BayesianNN::getBNNvalue(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ, 
            float pTZZ, 
            float YZZ
        )
{

    std::vector<double> mvainputvars;

    mvainputvars.push_back(costhetastar);
    mvainputvars.push_back(costheta1);
    mvainputvars.push_back(costheta2);
    mvainputvars.push_back(phi);
    mvainputvars.push_back(phistar1);
    mvainputvars.push_back(mZ1);
    mvainputvars.push_back(mZ2);
    mvainputvars.push_back(mZZ);


    if(_channel=="")   {
        std::cout << "Channel not set.  Set with nn.setChannel(CHANNEL), where CHANNEL = 4mu, 4e or 2e2mu." << std::endl;
    }
    
    double y = -1;

    // BNN trained with pt4l
    if(pTZZ > 0){
        mvainputvars.push_back(pTZZ);
        if(_channel=="4mu") y = bnn_4mu_8TeV_new_pt(mvainputvars, 50, 100);
        if(_channel=="4e") y = bnn_4e_8TeV_new_pt(mvainputvars, 50, 100);
        if(_channel=="2e2mu") y = bnn_2e2mu_8TeV_new_pt(mvainputvars, 50, 100);    
    }    
    // BNN with 8 variables (no pt4l)
    else if(_is7TeV)
    {
 //       std::cout << "7TeV" << std::endl;
        if(_channel=="4mu") y = bnn_4mu_7TeV_new(mvainputvars, 50, 100);
        if(_channel=="4e") y = bnn_4e_7TeV_new(mvainputvars, 50, 100);
        if(_channel=="2e2mu") y = bnn_2e2mu_7TeV_new(mvainputvars, 50, 100);    
    }
    else
    {
//        std::cout << "8TeV" << std::endl;
        if(_channel=="4mu") y = bnn_4mu_8TeV_new(mvainputvars, 50, 100);
        if(_channel=="4e") y = bnn_4e_8TeV_new(mvainputvars, 50, 100);
        if(_channel=="2e2mu") y = bnn_2e2mu_8TeV_new(mvainputvars, 50, 100);    
    }    
    
    return y;
       
}

double BayesianNN::getBNNvalueSPS(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ
        )
{

    std::vector<double> mvainputvars;

    mvainputvars.push_back(costhetastar);
    mvainputvars.push_back(costheta1);
    mvainputvars.push_back(costheta2);
    mvainputvars.push_back(phi);
    mvainputvars.push_back(phistar1);
    mvainputvars.push_back(mZ1);
    mvainputvars.push_back(mZ2);
    mvainputvars.push_back(mZZ);

    if(_channel=="")   {
        std::cout << "Channel not set.  Set with nn.setChannel(CHANNEL), where CHANNEL = 4mu, 4e or 2e2mu." << std::endl;
    }

    double y = -1;
    if(_channel=="4mu") y = Higgs0M_4mu(mvainputvars, 50, 100);
    if(_channel=="4e") y = Higgs0M_4e(mvainputvars, 50, 100);
    if(_channel=="2e2mu") y = Higgs0M_2e2mu(mvainputvars, 50, 100);    
  
    return y;
       
}




double BayesianNN::getBNNvalueGrav(
            float mZ1, 
            float mZ2, 
            float costhetastar,
            float costheta1, 
            float costheta2,
            float phi,
            float phistar1,           
            float mZZ
        )
{

    std::vector<double> mvainputvars;

    mvainputvars.push_back(costhetastar);
    mvainputvars.push_back(costheta1);
    mvainputvars.push_back(costheta2);
    mvainputvars.push_back(phi);
    mvainputvars.push_back(phistar1);
    mvainputvars.push_back(mZ1);
    mvainputvars.push_back(mZ2);
    mvainputvars.push_back(mZZ);

    if(_channel=="")   {
        std::cout << "Channel not set.  Set with nn.setChannel(CHANNEL), where CHANNEL = 4mu, 4e or 2e2mu." << std::endl;
    }

    double y = -1;
    if(_channel=="4mu") y = Graviton_4mu(mvainputvars, 50, 100);
    if(_channel=="4e") y = Graviton_4e(mvainputvars, 50, 100);
    if(_channel=="2e2mu") y = Graviton_2e2mu(mvainputvars, 50, 100);    
  
    return y;
       
}
// Expose classes and methods to Python




/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/ActionWithValue.h"
#include "core/ActionSet.h"
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "tools/Grid.h"
#include "tools/Exception.h"
#include "tools/File.h"
#include "tools/Matrix.h"
#include <memory>
#include <iostream>

using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS GAMBES
/*
*/
//+ENDPLUMEDOC

class GAMBES : public Bias {

private:
  unsigned int n_states_;
  vector< unsigned int > n_gaussians_;
  vector< double> resp;
  double scale_;
  double beta;
  unsigned int pace_;
  bool spline;
  vector<double> factor;
  std::string filename;
  vector<vector<double>> mu_k;
  vector<Matrix<double>> cov_k ;
  void getGaussians(); // called only once; initializes states_grid_pntr_ ; but that grid is not read from file but is computed here from function  
  double update_bias(std::vector<double>&, std::vector<double>&,std::vector<double>&);
  vector<double> get_factor();
  std::string bias_filename_;
  std::string states_filename_;
  vector<double> averages_k;
  double average;
  double bias_cutoff;
  bool cutoff;
  bool static_bias;
  double norm_gaussians;
  double bias_mean;
  unsigned int iter;
  double temp;
  double lambda;
  vector<double> static_factors;
  vector<double> grid_integration_weights;
  double energy_offset;
  vector<ActionWithValue*> biases;
    
public:
  explicit GAMBES(const ActionOptions&);
  void calculate();
  void prepare();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(GAMBES,"GAMBES")

void GAMBES::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","NSTATES","2","number of states present in the system to be studied");
  keys.add("compulsory","PACE","500","the frequency of updating the weight of each state");
  keys.add("compulsory","FILENAME","state","the suffix of the names of the files containing the Gaussians describing each state");
  keys.add("optional","TEMPERATURE","the temperature at which the simulation is being carried out");
  keys.add("optional","CUTOFF","the cutoff for bias");
  keys.add("optional","LAMBDA","the lambda for bias cutoff");
  keys.add("optional","STATIC_FACTORS","the constant factors for the Gaussians at which the simulation should be carried out if the bias is static");
  keys.addOutputComponent("_factors","default","two or more weighing factors for bias"
                          "these quantities will named with  the gaussian number followed by "
                          "the character string _factors. These quantities tell the user the value of the factor ");
  keys.addOutputComponent("_averages","AVERAGES","two or more the averages");

  keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the energy and forces due to the external potential");
  keys.addFlag("BIAS_CUTOFF",false,"to specify if there should be a bias cutoff");
  keys.addFlag("STATIC_BIAS",false,"specifies the bias is static");

}

GAMBES::GAMBES(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao),
  n_states_(0),
  n_gaussians_(0),
  resp(0),
  beta(1.0),
  factor(1.0),
  average(0),
  bias_cutoff(0),
  cutoff(false),
  norm_gaussians(0),
  iter(0),
  temp(0),
  lambda(0),
  energy_offset(0.0)
{ 
 
  unsigned int nargs = getNumberOfArguments();
  parse("NSTATES",n_states_);
  parse("FILENAME",filename);
  unsigned int tot_gaussians=0;

  for(unsigned int n=0; n<n_states_; n++) {
    std::unique_ptr<IFile> ifile(new IFile);
    std::string fname = filename+"."+std::to_string(n);
    if(ifile->FileExist(fname)) {
      ifile->open(fname);
      int component;
      int ngaussian = 0;

      while(ifile->scanField("ID",component)) {
        // get weight of Gaussian
        double w;
        ifile->scanField("WEIGHTS",w);
        resp.push_back(w);

        // get mean of the Gaussians
        vector<double> temp_mu(nargs);
        for(unsigned i=0; i<nargs; i++){
          ifile->scanField("mu_"+std::to_string(i),temp_mu[i]);
        }
        if( temp_mu.size()!=nargs) error(" Improper size of MU for state");
        mu_k.push_back(temp_mu);
         
        // get covariance of the Gaussians
        Matrix<double> temp_cov(nargs,nargs);
        unsigned int k=0;
        for(unsigned i=0; i<nargs; i++){ 
          for(unsigned j=0; j<nargs; j++){
            ifile->scanField("cov_"+std::to_string(i)+"_"+to_string(j),temp_cov(i,j));
            k++;
          }
        }
        if( k!=nargs*nargs) error(" Improper size of COVARIANCE for state ");
        cov_k.push_back(temp_cov);       
        // new line
        ifile->scanField();
        ngaussian++;
      }
      n_gaussians_.push_back(ngaussian);

    }
    else { error("Cannot find file "+fname+"\n"); }
    tot_gaussians+=n_gaussians_[n];
  }


  parse("PACE",pace_);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  spline=!nospline;
  static_bias=false;
  if(keywords.exists("STATIC_BIAS")){parseFlag("STATIC_BIAS",static_bias); }
  if(static_bias){
    parseVector("STATIC_FACTORS",static_factors);
    if(static_factors.size()!=n_states_) error("not enough values for STATIC_FACTORS");
  }

  if(keywords.exists("TEMPERATURE")){ parse("TEMPERATURE",temp);}

  cutoff=false;
  parseFlag("BIAS_CUTOFF",cutoff);
  if(cutoff){
    if(keywords.exists("CUTOFF")){parse("CUTOFF",bias_cutoff);}
    if(keywords.exists("LAMBDA")){parse("LAMBDA",lambda);}
    log.printf(" Bias cutoff applied at %lf with lambda %lf\n", bias_cutoff,lambda);
  }


  checkRead();

  log.printf(" Number of states in the potential %u\n",n_states_);
  for (unsigned n=0; n<n_states_; n++) {
    log.printf(" Number of gaussians in state %u is %u\n",n,n_gaussians_[n]);
  }
// add a component for factor
  for (unsigned n=0; n<n_states_; n++) {
    string s = std::to_string(n)+"_factors";
    addComponent(s);
    componentIsNotPeriodic(s);
    averages_k.push_back(0.0);

  }
  for (unsigned n=0; n<n_states_; n++) {
    string s = std::to_string(n)+"_averages";
    addComponent(s);
    componentIsNotPeriodic(s);
  }

  //construct biases from ActionWithValue with a component named bias
  vector<ActionWithValue*> tmpActions=plumed.getActionSet().select<ActionWithValue*>();
  for(unsigned i=0; i<tmpActions.size(); i++) if(tmpActions[i]->exists(tmpActions[i]->getLabel()+".bias")) biases.push_back(tmpActions[i]);
}
 
void GAMBES::prepare() {
  plumed.getAtoms().setCollectEnergy(true);
}

void GAMBES::calculate()
{
   
  beta = 1/plumed.getAtoms().getKbT();
  if(temp>0.0){ beta = 1/(temp*plumed.getAtoms().getKBoltzmann()) ; }
  bias_filename_="bias.";
  states_filename_="state.";

  unsigned nargs=getNumberOfArguments();
  vector<double> cv(nargs), der(nargs);
  for(unsigned i=0; i<nargs; i++) {cv[i]=getArgument(i);  }
  if(int(iter)==0 ){get_factor();}
  double bias;
  vector<double> u(n_states_,0.0);
  bias=update_bias(cv,der,u); 
  setBias(bias);

  for(unsigned i=0; i<nargs; i++) {
    const double f=-der[i];
    setOutputForce(i,f);
  }
  
  average=iter;
  if(iter==0){energy_offset= plumed.getAtoms().getEnergy();}
  double energy=plumed.getAtoms().getEnergy() - energy_offset;
  double energy1=0.0;
  for(unsigned i=0; i<biases.size(); i++) energy1+=biases[i]->getOutputQuantity("bias");
 // std::cout<<energy1<<" "<<bias<<std::endl;
  energy+=energy1;

  for(unsigned n=0; n<n_states_; n++) {
    double u_val =   -(std::log(u[n]))/beta;
    double diff = (- energy + u_val);
    averages_k[n]+=std::exp(beta*diff);
    string s = std::to_string(n)+"_averages" ;
    getPntrToComponent(s)->set(averages_k[n]/average);
  }

  iter+=1;
}



double GAMBES::update_bias(vector<double>& cv,vector<double>& cvder,vector<double>& u)
{ 
  double bias=0.0;
  unsigned int tot_gaussians=0;
  for(unsigned n=0; n<n_states_; n++) {
    tot_gaussians+=n_gaussians_[n];
  }
//  std::cout<<"1"<<std::endl;
  unsigned int nargs = getNumberOfArguments();
  vector<bool> isperiodic(nargs); //= BiasGrid_->getIsPeriodic();  
  vector<double> periodicity(nargs) ;

  for(unsigned i=0;i<nargs;i++){
    isperiodic[i]=getPntrToArgument(i)->isPeriodic();
    double min,max; getPntrToArgument(i)->getDomain(min,max);
    periodicity[i]=max-min;
  } //periodicity made more general
  vector<Matrix<double>> invcov_k(tot_gaussians, Matrix<double>(nargs,nargs)) ;
  vector<double> P_k(tot_gaussians,0.0);
  vector<vector<double>> der_P_k(tot_gaussians,vector<double>(nargs,0.0));
  for(unsigned int n=0; n<tot_gaussians; n++){
    // get inverse covariance matrix
    Invert(cov_k[n],invcov_k[n]);
    
    // get determinants of each covariance matrix
    Matrix<double> dec_cov(nargs,nargs) ;
    cholesky(cov_k[n],dec_cov);
    double determinant=1.0;
    for(unsigned int i=0; i<nargs; i++){ determinant*=dec_cov(i,i);}
    double normalize = std::pow(2*pi,nargs);
    normalize= 1/(sqrt(normalize)*determinant);
    if(n==0){norm_gaussians=normalize;}
    
    vector<double> diff(nargs) ;
    for(unsigned int i=0; i<nargs; i++){
      diff[i] = cv[i]-mu_k[n][i];
      if(isperiodic[i]){
        if(diff[i]>0.5*periodicity[i]){ diff[i] -= periodicity[i] ; } 
        if(diff[i]<=-0.5*periodicity[i]){ diff[i] += periodicity[i] ; } 
      }
    }
    Matrix<double> diff_(1,nargs);
    diff_.setFromVector(diff);
    Matrix<double> transpose_diff_ ;
    transpose(diff_,transpose_diff_);
    Matrix<double> out1;
    mult(diff_, invcov_k[n], out1);
    Matrix<double> out2;
    mult(invcov_k[n],transpose_diff_, out2);
    Matrix<double> out;
    mult(out1,transpose_diff_, out);
    P_k[n] = std::exp(-0.5*out(0,0));

    for(unsigned int i=0; i<nargs; i++){
        der_P_k[n][i]+=-0.5*P_k[n]*(out1(0,i)+out2(i,0));
    }
    //cutoff
    if(cutoff==true){
      double cut_at = std::exp(-beta*bias_cutoff);
      double sw_val = std::exp(lambda*(P_k[n]-cut_at));
      double bias_sw = 1/(1+sw_val);
      double mult = bias_sw + lambda*std::pow(bias_sw,2)*P_k[n];
      for(unsigned i=0;i<nargs;i++){
        der_P_k[n][i]=der_P_k[n][i]*mult;
      }
      P_k[n]= P_k[n]*bias_sw + cut_at;
    }
  }

  //calculate bias
  if(int(iter)%pace_==0 and !static_bias){ factor = get_factor(); }
  else if(int(iter)==0 and static_bias){ factor = get_factor();}

  vector<double> expu(n_states_,0.0);
  vector<vector<double>> der_in(n_states_,vector<double>(nargs,0.0));
  int count=0;
  for(unsigned n=0; n<n_states_; n++) {
    string s = std::to_string(n)+"_factors" ;
    getPntrToComponent(s)->set(factor[n]);
    for(unsigned k=0; k<n_gaussians_[n]; k++){
      u[n]+=resp[count]*P_k[count];
      for(unsigned i=0;i<nargs;i++){
        der_in[n][i]+=resp[count]*der_P_k[count][i];
      }
      count++;
    }
    double u_val = -(std::log(u[n]))/beta;
    expu[n]=factor[n]*std::exp(-beta*u_val);
    bias+= expu[n];
  }
  double nbias;
  bias = +(std::log(bias))/beta;
  nbias = bias + bias_cutoff ;//- (std::log(norm_gaussians*n_states_))/beta;

  for(unsigned i=0;i<nargs;i++){
    double der=0.0;
    for(unsigned n=0; n<n_states_; n++) {
        der += factor[n]*der_in[n][i] ;
    }
    cvder[i]=std::exp(-beta*bias)*der/(beta) ;
  }

  return nbias;
}

// calculate the factors z_1/z_k
vector<double> GAMBES::get_factor()
{ 
  vector<double> factor_c(n_states_);
  vector<double> expectation(n_states_, 1.0);
  if(!static_bias){
    if(int(iter)==0){
      for(unsigned n=0; n<n_states_; n++) {
        expectation[n]=1.0; 
       }
    }
    else{
     for(unsigned n=0; n<n_states_; n++) {
       expectation[n]=averages_k[n]/average ;
     }
   }
    for(unsigned n=0; n<n_states_; n++) {
      factor_c[n]=expectation[0]/expectation[n];
     }
    }
  else{
     for(unsigned n=0; n<n_states_; n++) {
       factor_c[n]=static_factors[n];
     }
  }
  return factor_c;
}

}
}



#include <TFile.h>
#include <fstream>
#include <RooFit.h>
#include <RooArgSet.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TParameter.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCrystalBall.h>
#include <RooExtendPdf.h>
#include <RooAbsData.h>
int main(){
TFile input("/Users/imb/Downloads/dataset_region2_envelope_2023.root");
auto original=input.Get<RooWorkspace>("w");
RooWorkspace out("dcb_workspace");
out.import(*original->data("data_obs_resonant"));
auto mass=out.var("mass");
RooRealVar mu("mu","mu",9.3842240423278422,6,12); mu.setError(0.012114751194227225);
RooRealVar sigma("sigma","sigma",0.15923875753675959,0.005,2); sigma.setError(0.0050968555709650226);
RooRealVar alphaL("alphaL","alphaL",1.8172703241556882,0.1,10); alphaL.setError(0.35257971803772892);
RooRealVar nL("nL","nL",14.387331988483435,1.01,100); nL.setError(8.2390792126211707);
RooRealVar alphaR("alphaR","alphaR",4.0714474204075941,0.1,10); alphaR.setError(0.17612966911347905);
RooRealVar nR("nR","nR",1.0100000085245078,1.01,100); nR.setError(5.2043460515316244);
RooRealVar yield_("yield_","yield_",915998.73359843483,1,19679025.08158935); yield_.setError(43589.677949345962);
RooCrystalBall dcb("dcb","Double-sided Crystal Ball",*mass,mu,sigma,alphaL,nL,alphaR,nR);
RooExtendPdf model("dcb_extended","Extended DCB",dcb,yield_); out.import(model);
TFile output("outputs/dcb_fit.root","RECREATE"); out.Write();
TMatrixDSym covariance(7);
covariance(0,0)=0.00014676800761023482;
covariance(0,1)=-4.9151715837855461e-05;
covariance(0,2)=-0.0029639554626401977;
covariance(0,3)=0.056206633349231254;
covariance(0,4)=-8.9989085615568246e-05;
covariance(0,5)=-1.0185394211140352e-09;
covariance(0,6)=-3.8817238491574644;
covariance(1,0)=-4.9151715837855461e-05;
covariance(1,1)=2.5978729063477745e-05;
covariance(1,2)=0.0012340506228577033;
covariance(1,3)=-0.022469665021468049;
covariance(1,4)=8.450937993185319e-05;
covariance(1,5)=1.1033732460603678e-09;
covariance(1,6)=-11.347142590481551;
covariance(2,0)=-0.0029639554626401977;
covariance(2,1)=0.0012340506228577033;
covariance(2,2)=0.12468077911580143;
covariance(2,3)=-2.7281271315073048;
covariance(2,4)=0.0038239070845808771;
covariance(2,5)=4.1589711932389729e-08;
covariance(2,6)=422.61684312348939;
covariance(3,0)=0.056206633349231254;
covariance(3,1)=-0.022469665021468049;
covariance(3,2)=-2.7281271315073048;
covariance(3,3)=69.267686434421051;
covariance(3,4)=-0.067197939649385183;
covariance(3,5)=-7.1424176320663841e-07;
covariance(3,6)=-2338.90493403755;
covariance(4,0)=-8.9989085615568246e-05;
covariance(4,1)=8.450937993185319e-05;
covariance(4,2)=0.0038239070845808771;
covariance(4,3)=-0.067197939649385183;
covariance(4,4)=0.031035294133790105;
covariance(4,5)=1.4420013552365514e-06;
covariance(4,6)=525.46717402075194;
covariance(5,0)=-1.0185394211140352e-09;
covariance(5,1)=1.1033732460603678e-09;
covariance(5,2)=4.1589711932389729e-08;
covariance(5,3)=-7.1424176320663841e-07;
covariance(5,4)=1.4420013552365514e-06;
covariance(5,5)=1.8065811476639622e-07;
covariance(5,6)=0.00073232514484804744;
covariance(6,0)=-3.8817238491574644;
covariance(6,1)=-11.347142590481551;
covariance(6,2)=422.61684312348939;
covariance(6,3)=-2338.90493403755;
covariance(6,4)=525.46717402075194;
covariance(6,5)=0.00073232514484804744;
covariance(6,6)=1900130046.8879077;
covariance.Write("covariance_mu_sigma_alphaL_nL_alphaR_nR_yield");
TParameter<double>("chi2",559.96308797469794).Write(); TParameter<int>("ndof",314).Write();
TNamed("method","Chi-square with stored SumW2 errors; 129 zero-variance bins omitted; bin-integrated expectations. nR at lower bound; tail error unreliable.").Write();
output.Close();
TFile check("outputs/dcb_fit.root"); auto cw=check.Get<RooWorkspace>("dcb_workspace");
if(!cw || !cw->pdf("dcb") || cw->data("data_obs_resonant")->numEntries()!=450) return 2;
auto cx=cw->var("mass"); RooArgSet norm(*cx); std::ofstream validation("work/root_integrals.csv");validation.precision(17);
for(int i=0;i<450;i++){ cx->setRange("binCheck",6.+i*6./450,6.+(i+1)*6./450);auto integral=cw->pdf("dcb")->createIntegral(norm,RooFit::NormSet(norm),RooFit::Range("binCheck"));validation<<integral->getVal()*cw->var("yield_")->getVal()<<"\n";delete integral; }
return 0;}

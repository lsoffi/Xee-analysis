#include <TFile.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooMultiPdf.h>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(){TFile f("/Users/imb/Downloads/Xee_ee_4_allYears.root");auto w=(RooWorkspace*)f.Get("w");auto x=w->var("mass");auto d=w->data("data_obs");auto sig=w->pdf("shapeSig_Zd_ch1");auto bg=((RooMultiPdf*)w->pdf("shapeBkg_dy_ch1"))->getPdf(1);RooArgSet obs(*x);std::ofstream o("work/xee_2022_nofit.csv");o<<std::setprecision(16)<<"low,high,data,signal_fraction,background_fraction,upsilon1s_fraction,upsilon2s_fraction\n";double sum=0;for(int i=0;i<d->numEntries();i++){auto row=d->get(i);if(std::string(row->getCatLabel("CMS_channel"))!="ch1")continue;double m=row->getRealValue("mass");double lo=m-6.0/900,hi=m+6.0/900;x->setRange("bin",lo,hi);auto si=sig->createIntegral(obs,RooFit::NormSet(obs),RooFit::Range("bin"));auto bi=bg->createIntegral(obs,RooFit::NormSet(obs),RooFit::Range("bin"));auto u1=w->pdf("shapeBkg_upsilon1s_ch1")->createIntegral(obs,RooFit::NormSet(obs),RooFit::Range("bin"));auto u2=w->pdf("shapeBkg_upsilon2s_ch1")->createIntegral(obs,RooFit::NormSet(obs),RooFit::Range("bin"));o<<lo<<","<<hi<<","<<d->weight()<<","<<si->getVal()<<","<<bi->getVal()<<","<<u1->getVal()<<","<<u2->getVal()<<"\n";sum+=d->weight();delete si;delete bi;delete u1;delete u2;}std::ofstream meta("work/xee_2022_nofit_meta.txt");meta<<std::setprecision(16)<<"Data "<<sum<<"\nSignal nominal "<<w->function("n_exp_binch1_proc_Zd")->getVal()<<"\nContinuum nominal "<<w->function("n_exp_binch1_proc_dy")->getVal()<<"\nr "<<w->var("r")->getVal()<<"\n";for(auto a:*sig->getParameters(obs)){auto v=dynamic_cast<RooRealVar*>(a);if(v&&!v->isConstant())meta<<v->GetName()<<" "<<v->getVal()<<"\n";}for(auto proc:{"Zd","dy","upsilon1s","upsilon2s"})meta<<proc<<" "<<w->function((std::string("n_exp_binch1_proc_")+proc).c_str())->getVal()<<"\n";meta.close();o.close();}

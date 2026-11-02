#include <TFile.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooAbsData.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooFitResult.h>
#include <RooMultiPdf.h>
#include <RooMinimizer.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace RooFit;
int main(){TFile f("/Users/imb/Downloads/Xee_ee_4_allYears.root");auto w=(RooWorkspace*)f.Get("w");auto x=w->var("mass");RooArgSet obs(*x);auto d=w->data("data_obs");RooDataHist dh("data_2022","data_2022",obs);for(int i=0;i<d->numEntries();i++){auto row=d->get(i);if(std::string(row->getCatLabel("CMS_channel"))!="ch1")continue;x->setVal(row->getRealValue("mass"));dh.add(obs,d->weight(),d->weight());}
auto bg=((RooMultiPdf*)w->pdf("shapeBkg_dy_ch1"))->getPdf(1);w->var("CMS_EXO25020_bkgScale_upsilon1s_2022")->setConstant(true);auto r=w->var("r");r->setRange(-30,30);r->setVal(1);r->setError(.1);
RooArgList bp(*bg,*w->pdf("shapeBkg_upsilon1s_ch1"),*w->pdf("shapeBkg_upsilon2s_ch1"));RooArgList bn;for(auto p:{"dy","upsilon1s","upsilon2s"})bn.add(*w->function((std::string("n_exp_binch1_proc_")+p).c_str()));RooAddPdf bmodel("background_2022","",bp,bn);RooArgList sp(*w->pdf("shapeSig_Zd_ch1"));sp.add(bp);RooArgList sn(*w->function("n_exp_binch1_proc_Zd"));sn.add(bn);RooAddPdf smodel("signal_background_2022","",sp,sn);
RooArgSet constraints,globals,nuis;for(auto name:{"CMS_eff_e_id","CMS_eff_e_reco","CMS_eff_trigger_e","lumi_1","CMS_scale_e","CMS_EXO25020_signalModelSigmaNuisance_2022"}){constraints.add(*w->pdf((std::string(name)+"_Pdf").c_str()));globals.add(*w->var((std::string(name)+"_In").c_str()));nuis.add(*w->var(name));}
auto params=smodel.getParameters(obs);std::ofstream table("outputs/xee_2022_frozenY1_parameters.csv");table<<std::setprecision(15)<<"stage,parameter,value,error,error_lo,error_hi,constant\n";auto dump=[&](std::string stage){for(auto a:*params){auto v=dynamic_cast<RooRealVar*>(a);if(v)table<<stage<<","<<v->GetName()<<","<<v->getVal()<<","<<v->getError()<<","<<v->getAsymErrorLo()<<","<<v->getAsymErrorHi()<<","<<v->isConstant()<<"\n";}};dump("prefit");w->saveSnapshot("initial2022",*params);TFile out("outputs/xee_2022_frozenY1_roofit.root","RECREATE");std::ofstream summary("outputs/xee_2022_frozenY1_summary.txt");summary<<std::setprecision(12)<<"Native RooFit diagnostic, NOT executed with Combine. 2022/ch1, mass6-12, Chebyshev only.\nUps1 scale frozen at1, r range[-30,30], Minuit2/migrad strategy0.\nBin-center likelihood, original integer weights interpreted as counts; no SumW2 correction.\nHESSE nuisance errors, MINOS r errors, not Combine robustFit scan.\n";
for(int stage=0;stage<2;stage++){w->loadSnapshot("initial2022");r->setVal(stage?1:0);r->setConstant(!stage);auto& model=stage?smodel:bmodel;auto nll=model.createNLL(dh,Extended(),ExternalConstraints(stage?constraints:RooArgSet()),GlobalObservables(stage?globals:RooArgSet()),EvalBackend("legacy"),Offset(true));RooMinimizer minim(*nll);minim.setMinimizerType("Minuit2");minim.setStrategy(0);minim.setPrintLevel(-1);minim.setMaxFunctionCalls(100000);minim.setEps(.1);int ms=minim.minimize("Minuit2","migrad");int hs=minim.hesse();int mins=stage?minim.minos(RooArgSet(*r)):0;auto fit=minim.save();std::string tag=stage?"fit_s":"fit_b";fit->Write(tag.c_str());fit->Print("v");dump(tag);summary<<tag<<" minimize="<<ms<<" hesse="<<hs<<" minos="<<mins<<" status="<<fit->status()<<" covQual="<<fit->covQual()<<" edm="<<fit->edm()<<" nll_offset="<<fit->minNll()<<"\n";
std::ofstream curves("outputs/xee_2022_frozenY1_"+tag+"_bins.csv");curves<<std::setprecision(15)<<"mass,data,signal,dy,upsilon1s,upsilon2s,total\n";for(int i=0;i<450;i++){dh.get(i);double count=dh.weight();x->setVal(6+(i+.5)*6./450);double y=0;curves<<x->getVal()<<","<<count;for(int j=0;j<4;j++){double v=j==0&&!stage?0:((RooAbsPdf*)sp.at(j))->getVal(obs)*((RooAbsReal*)sn.at(j))->getVal()*6./450;curves<<","<<v;y+=v;}curves<<","<<y<<"\n";}curves.close();for(auto a:sn)summary<<tag<<" "<<a->GetName()<<"="<<((RooAbsReal*)a)->getVal()<<"\n";delete fit;}
summary.close();table.close();out.Close();}

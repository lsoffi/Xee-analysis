#include <TFile.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooAbsData.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooMultiPdf.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace RooFit;
int main(){TFile f("/Users/imb/Downloads/Xee_ee_4_allYears.root");auto w=(RooWorkspace*)f.Get("w");auto x=w->var("mass");RooArgSet obs(*x);auto d=w->data("data_obs");RooDataHist dh("data_2022","data_2022",obs);for(int i=0;i<d->numEntries();i++){auto row=d->get(i);if(std::string(row->getCatLabel("CMS_channel"))!="ch1")continue;x->setVal(row->getRealValue("mass"));dh.add(obs,d->weight(),d->weight());}
auto bg=((RooMultiPdf*)w->pdf("shapeBkg_dy_ch1"))->getPdf(1);RooArgList pdfs(*bg,*w->pdf("shapeBkg_upsilon1s_ch1"),*w->pdf("shapeBkg_upsilon2s_ch1"));RooArgList norms;for(auto p:{"dy","upsilon1s","upsilon2s"})norms.add(*w->function((std::string("n_exp_binch1_proc_")+p).c_str()));RooAddPdf model("background_2022","Chebyshev + Upsilon1S + Upsilon2S",pdfs,norms);auto pars=model.getParameters(obs);pars->Print("v");std::ofstream report("outputs/xee_2022_bonly_fit.txt");report<<std::setprecision(12)<<"2022 only, mass 6-12 GeV, 450 bins, binned extended Poisson likelihood.\nNative RooFit using original workspace components, Chebyshev index 1; no signal, no sideband mask.\nInteger stored weights interpreted as bin counts; SumW2 correction disabled.\nINITIAL PARAMETERS\n";for(auto a:*pars){auto v=dynamic_cast<RooRealVar*>(a);if(v)report<<v->GetName()<<" "<<v->getVal()<<" constant="<<v->isConstant()<<" range="<<v->getMin()<<","<<v->getMax()<<"\n";}
auto fit=model.fitTo(dh,Save(),Extended(),EvalBackend("legacy"),IntegrateBins(1e-4),Strategy(2),Minimizer("Minuit2","migrad"),MaxCalls(100000),Offset(true),PrintLevel(-1),SumW2Error(false));fit->Print("v");report<<"\nRESULT\nstatus "<<fit->status()<<" covQual "<<fit->covQual()<<" edm "<<fit->edm()<<" nll_offset "<<fit->minNll()<<"\n";for(auto a:fit->floatParsFinal()){auto v=(RooRealVar*)a;report<<v->GetName()<<" "<<v->getVal()<<" +/- "<<v->getError()<<"\n";}
TFile out("outputs/xee_2022_bonly_fit.root","RECREATE");fit->Write("fit_b");dh.Write();model.Write();std::ofstream csv("outputs/xee_2022_bonly_fit_bins.csv");csv<<std::setprecision(16)<<"low,high,data,dy,upsilon1s,upsilon2s,total\n";double dev=0;for(int i=0;i<450;i++){dh.get(i);double count=dh.weight(),lo=6+i*6.0/450,hi=6+(i+1)*6.0/450;x->setRange("bin",lo,hi);double total=0;csv<<lo<<","<<hi<<","<<count;for(int j=0;j<3;j++){auto integral=((RooAbsPdf*)pdfs.at(j))->createIntegral(obs,NormSet(obs),Range("bin"));double val=integral->getVal()*((RooAbsReal*)norms.at(j))->getVal();total+=val;csv<<","<<val;delete integral;}csv<<","<<total<<"\n";if(total>0)dev+=2*(total-count+(count>0?count*log(count/total):0));}report<<"Poisson deviance "<<dev<<"\nnominal_dof "<<450-fit->floatParsFinal().getSize()<<"\n";for(auto a:norms)report<<a->GetName()<<" yield "<<((RooAbsReal*)a)->getVal()<<"\n";report.close();csv.close();out.Close();}

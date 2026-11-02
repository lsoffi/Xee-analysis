#include <TFile.h>
#include <RooWorkspace.h>
#include <RooAbsData.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <fstream>
void extract(){
 TFile f("/Users/imb/Downloads/dataset_region2_envelope_2023.root");
 auto w=f.Get<RooWorkspace>("w"); if(!w){throw std::runtime_error("workspace w missing");}
 w->Print();
 auto d=w->data("data_obs_resonant"); if(!d)throw std::runtime_error("dataset missing");
 d->Print("v"); auto obs=d->get(); obs->Print("v");
 std::ofstream o("work/data_obs_resonant.csv");o.precision(17);
 o<<"entry";for(auto a:*obs) o<<","<<a->GetName();o<<",weight,weight_error_lo,weight_error_hi\n";
 for(int i=0;i<d->numEntries();i++){auto r=d->get(i);o<<i;for(auto a:*r){auto v=dynamic_cast<RooRealVar*>(a);if(!v)throw std::runtime_error("non-real observable");o<<","<<v->getVal();}double lo,hi;d->weightError(lo,hi,RooAbsData::SumW2);o<<","<<d->weight()<<","<<lo<<","<<hi<<"\n";}
 std::ofstream m("work/data_metadata.txt");m.precision(17);m<<"class "<<d->ClassName()<<"\nentries "<<d->numEntries()<<"\nsumWeights "<<d->sumEntries()<<"\nweighted "<<d->isWeighted()<<"\n";for(auto a:*obs){auto v=dynamic_cast<RooRealVar*>(a);if(v)m<<"observable "<<v->GetName()<<" min "<<v->getMin()<<" max "<<v->getMax()<<" bins "<<v->getBins()<<" unit "<<v->getUnit()<<"\n";}
}

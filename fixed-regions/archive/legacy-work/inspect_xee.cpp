#include <TFile.h>
#include <TKey.h>
#include <RooWorkspace.h>
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <iostream>
int main(){TFile f("/Users/imb/Downloads/Xee_ee_4_allYears.root");f.ls();for(auto k:*f.GetListOfKeys()){auto key=(TKey*)k;auto o=key->ReadObj();auto w=dynamic_cast<RooWorkspace*>(o);if(w){w->Print();std::cout<<"\nALL VARIABLES\n";w->allVars().Print("v");std::cout<<"\nDATA DETAILS\n";for(auto d:w->allData()){d->Print("v");d->get()->Print("v");}std::cout<<"\nPDF DETAILS\n";for(auto p:w->allPdfs())p->Print("v");}}}

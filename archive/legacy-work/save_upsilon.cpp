#include <TFile.h>
#include <fstream>
#include <TH1D.h>
#include <RooDataHist.h>
#include <RooArgList.h>
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
RooWorkspace out("upsilon_dcb_workspace");
RooRealVar massvar("mass","mass",8.8,9.813333333333334);
massvar.setBins(38);
TH1D histogram("rebinned","rebinned",38,8.8,9.813333333333334);histogram.Sumw2();
histogram.SetBinContent(1,6662.3657377529089);histogram.SetBinError(1,3559.3016368205394);
histogram.SetBinContent(2,19887.943898530888);histogram.SetBinError(2,7407.7731240281682);
histogram.SetBinContent(3,17471.539474746005);histogram.SetBinError(3,6421.4406919242683);
histogram.SetBinContent(4,21438.043209898409);histogram.SetBinError(4,7909.7424806795525);
histogram.SetBinContent(5,21325.880067245987);histogram.SetBinError(5,7210.9218210410618);
histogram.SetBinContent(6,20086.830602783884);histogram.SetBinError(6,8173.7311479727296);
histogram.SetBinContent(7,19871.937492159974);histogram.SetBinError(7,7063.8454005670892);
histogram.SetBinContent(8,17345.569182900137);histogram.SetBinError(8,6457.4114348746243);
histogram.SetBinContent(9,18058.817450702918);histogram.SetBinError(9,6022.3924250557311);
histogram.SetBinContent(10,9473.3074588461022);histogram.SetBinError(10,3402.0789760136167);
histogram.SetBinContent(11,32284.285752414326);histogram.SetBinError(11,8466.5455293236919);
histogram.SetBinContent(12,18422.420035276664);histogram.SetBinError(12,5823.1720227243632);
histogram.SetBinContent(13,35837.149038527088);histogram.SetBinError(13,9165.5555967563214);
histogram.SetBinContent(14,37122.132742387417);histogram.SetBinError(14,8937.434084464343);
histogram.SetBinContent(15,36451.044989982489);histogram.SetBinError(15,9638.2212376820553);
histogram.SetBinContent(16,29886.796811606429);histogram.SetBinError(16,7739.913516579345);
histogram.SetBinContent(17,34512.527034021521);histogram.SetBinError(17,8536.0997928051474);
histogram.SetBinContent(18,46014.418854095034);histogram.SetBinError(18,9694.7457346232695);
histogram.SetBinContent(19,48746.67509855891);histogram.SetBinError(19,11150.114736253861);
histogram.SetBinContent(20,52464.898783741679);histogram.SetBinError(20,11147.311898550919);
histogram.SetBinContent(21,50679.561065304239);histogram.SetBinError(21,10125.08867267381);
histogram.SetBinContent(22,41113.598208363153);histogram.SetBinError(22,9026.4493406870297);
histogram.SetBinContent(23,78154.851763997169);histogram.SetBinError(23,14076.388049167215);
histogram.SetBinContent(24,66971.778322497761);histogram.SetBinError(24,12764.446099990073);
histogram.SetBinContent(25,66571.273693996773);histogram.SetBinError(25,12624.635097467628);
histogram.SetBinContent(26,38983.379254030573);histogram.SetBinError(26,8331.7221298539698);
histogram.SetBinContent(27,62938.455268359321);histogram.SetBinError(27,11916.111817848419);
histogram.SetBinContent(28,47532.521706382839);histogram.SetBinError(28,9375.9085769887242);
histogram.SetBinContent(29,29236.252959417874);histogram.SetBinError(29,7486.1554583881434);
histogram.SetBinContent(30,36877.094394332147);histogram.SetBinError(30,9492.9263586437992);
histogram.SetBinContent(31,30886.706783320824);histogram.SetBinError(31,8494.0884565211345);
histogram.SetBinContent(32,46810.422525526927);histogram.SetBinError(32,11058.900748099859);
histogram.SetBinContent(33,32629.067774032621);histogram.SetBinError(33,8043.4407673868591);
histogram.SetBinContent(34,29250.900606210045);histogram.SetBinError(34,7799.7493170518865);
histogram.SetBinContent(35,19825.527223347261);histogram.SetBinError(35,5809.0888840916923);
histogram.SetBinContent(36,18299.351985283723);histogram.SetBinError(36,6049.7149007398302);
histogram.SetBinContent(37,22383.5684799267);histogram.SetBinError(37,6572.0801970291704);
histogram.SetBinContent(38,14799.858363078423);histogram.SetBinError(38,6669.1631536129225);
RooDataHist data("data_obs_resonant","signal, rebin 2",RooArgList(massvar),RooFit::Import(histogram));out.import(data);auto mass=out.var("mass");
RooRealVar mu("mu","mu",9.4136690501689593,9.1,9.65); mu.setError(0.019966218295381033);
RooRealVar sigma("sigma","sigma",0.02904111826912617,0.01,0.5); sigma.setError(0.0075053316176813007);
RooRealVar alphaL("alphaL","alphaL",0.11659810308683824,0.1,10.0); alphaL.setError(0.055789640472217988);
RooRealVar nL("nL","nL",2.2247210562145345,1.01,100.0); nL.setError(2.6146116094552445);
RooRealVar alphaR("alphaR","alphaR",0.10001359955646637,0.1,10.0); alphaR.setError(7.3886490987944331);
RooRealVar nR("nR","nR",99.999289048766784,1.01,100.0); nR.setError(97.835792233909743);
RooRealVar yield_("yield_","yield_",1246326.4778383498,63865.43770467936,6386543.770467936); yield_.setError(51987.478984249523);
RooCrystalBall dcb("dcb","Double-sided Crystal Ball",*mass,mu,sigma,alphaL,nL,alphaR,nR);
RooExtendPdf model("dcb_extended","Extended DCB",dcb,yield_); out.import(model);
TFile output("outputs/upsilon_dcb_rebin2.root","RECREATE"); out.Write();
TMatrixDSym covariance(7);
covariance(0,0)=0.00039936654079014017;
covariance(0,1)=-3.5622808376085721e-05;
covariance(0,2)=-0.00069887936296135998;
covariance(0,3)=0.024188867324727577;
covariance(0,4)=4.5195579912694264e-06;
covariance(0,5)=6.3985448493784393e-05;
covariance(0,6)=-16.511191312241181;
covariance(1,0)=-3.5622808376085721e-05;
covariance(1,1)=5.6448345329073529e-05;
covariance(1,2)=0.00018193990569897135;
covariance(1,3)=0.0013192255456714065;
covariance(1,4)=0.00010652357808238036;
covariance(1,5)=-0.00010735035561288868;
covariance(1,6)=-3.3505378530895928;
covariance(2,0)=-0.00069887936296135998;
covariance(2,1)=0.00018193990569897135;
covariance(2,2)=0.0031323703049761973;
covariance(2,3)=-0.12122119143165185;
covariance(2,4)=0.00042046128169862001;
covariance(2,5)=-0.0003669705410752342;
covariance(2,6)=65.124343651258812;
covariance(3,0)=0.024188867324727577;
covariance(3,1)=0.0013192255456714065;
covariance(3,2)=-0.12122119143165185;
covariance(3,3)=6.9715356934383541;
covariance(3,4)=9.1815693332684628e-05;
covariance(3,5)=-0.0017428007098845143;
covariance(3,6)=-4250.7444790649552;
covariance(4,0)=4.5195579912694264e-06;
covariance(4,1)=0.00010652357808238036;
covariance(4,2)=0.00042046128169862001;
covariance(4,3)=9.1815693332684628e-05;
covariance(4,4)=0.00032700026348252869;
covariance(4,5)=1.0666097781522482e-05;
covariance(4,6)=-0.12262576460990419;
covariance(5,0)=6.3985448493784393e-05;
covariance(5,1)=-0.00010735035561288868;
covariance(5,2)=-0.0003669705410752342;
covariance(5,3)=-0.0017428007098845143;
covariance(5,4)=1.0666097781522482e-05;
covariance(5,5)=82.534849240143174;
covariance(5,6)=9.1820005317010995;
covariance(6,0)=-16.511191312241181;
covariance(6,1)=-3.3505378530895928;
covariance(6,2)=65.124343651258812;
covariance(6,3)=-4250.7444790649552;
covariance(6,4)=-0.12262576460990419;
covariance(6,5)=9.1820005317010995;
covariance(6,6)=2703098661.1056304;
covariance.Write("covariance_mu_sigma_alphaL_nL_alphaR_nR_yield");
TParameter<double>("weighted_deviance",33.00333611913383).Write(); TParameter<int>("nominal_ndof",31).Write();
TNamed("method","Weighted Bohm-Zech likelihood, signal-only rebin 2; tail parameters at limits; see MINOS errors in JSON").Write();
output.Close();
return 0;}

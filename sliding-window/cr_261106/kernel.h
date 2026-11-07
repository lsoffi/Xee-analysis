#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <string>

namespace ThreeCR {
struct Shape { double mean,sigma,al,ar,nl,nr; };
struct Objective {
  std::vector<double> edges, counts;
  std::vector<double> quadX,quadW;
  std::vector<std::vector<double>> basis;
  std::array<Shape,2> shape;
  int degree=0; // -1 exponential; >=0 Bernstein
  bool floatWidth=true;
  int tailCount=0;
  bool independentAlpha=false;
  bool fractionalN=false;
  int nbase() const {return (floatWidth?7:5)+tailCount;}
  int npar() const {return nbase()+((degree==-3||degree==-4)?2:((degree==-2)?3:((degree<0)?1:degree)));}
  static double quadraticMinimum(double a,double b) {
    double v=std::min(1-a+b,1+a+b);
    if(b>0 && std::abs(a/(2*b))<1)v=std::min(v,1-a*a/(4*b));
    return v;
  }
  static double choose(int n,int k) {
    double r=1; for(int j=1;j<=k;++j) r*=double(n-j+1)/j; return r;
  }
  static double bcdf(double z,int d,int i) {
    double r=0; for(int k=i+1;k<=d+1;++k) r+=choose(d+1,k)*std::pow(z,k)*std::pow(1-z,d+1-k);
    return r/(d+1);
  }
  void prepare() {
    basis.clear(); if(degree<0)return;
    const double a=edges.front(), scale=edges.back()-a;
    for(size_t j=0;j<counts.size();++j) {
      std::vector<double> r;
      for(int i=0;i<=degree;++i)r.push_back(bcdf((edges[j+1]-a)/scale,degree,i)-bcdf((edges[j]-a)/scale,degree,i));
      basis.push_back(r);
    }
  }
  static double cdf(double x,double mu,double sigma,const Shape&s) {
    const double t=(x-mu)/sigma, c=std::sqrt(std::acos(-1.)/2.);
    const double left=std::exp(-.5*s.al*s.al)*s.nl/s.al/(s.nl-1);
    if(t<-s.al)return left*std::pow((s.nl/s.al)/(s.nl/s.al-s.al-t),s.nl-1);
    if(t<=s.ar)return left+c*(std::erf(t/std::sqrt(2.))+std::erf(s.al/std::sqrt(2.)));
    const double mid=c*(std::erf(s.al/std::sqrt(2.))+std::erf(s.ar/std::sqrt(2.)));
    const double right=std::exp(-.5*s.ar*s.ar)*s.nr/s.ar/(s.nr-1);
    return left+mid+right*(1-std::pow((s.nr/s.ar)/(s.nr/s.ar-s.ar+t),s.nr-1));
  }
  std::vector<std::vector<double>> components(const double*p) const {
    const int nb=counts.size(), base=nbase();
    std::vector<std::vector<double>> out(3,std::vector<double>(nb));
    std::vector<double> coefficients;
    double sum=1;
    if(degree>=0) {
      coefficients.push_back(1);
      for(int i=1;i<=degree;++i){coefficients.push_back(std::exp(p[base+i-1]));sum+=coefficients.back();}
    }
    const double a=edges.front(), width=edges.back()-a;
    if(degree==-3 && !(quadraticMinimum(p[base],p[base+1])>1.e-8))return {};
    std::vector<double> erfBins;
    double erfNorm=0;
    if(degree==-2) {
      for(int j=0;j<nb;++j) {
        double v=0;
        const double l=(edges[j]-a)/width,h=(edges[j+1]-a)/width;
        for(size_t q=0;q<quadX.size();++q) {
          double z=.5*(l+h)+.5*(h-l)*quadX[q];
          v+=quadW[q]*std::erf((p[base]-z)/std::exp(p[base+1]))*std::exp(-z/std::exp(p[base+2]));
        }
        v*=.5*(h-l);erfBins.push_back(v);erfNorm+=v;
      }
      if(!(erfNorm>0))return {};
    }
    for(int j=0;j<nb;++j) {
      double frac=0;
      if(degree==-3||degree==-4) {
        double ca=p[base],cb=p[base+1];
        if(degree==-4) {
          const double u=ca,v=cb,center=1+u/4+3*v/4;
          if(u<0||v<0)return {};
          ca=-(u+v)/(2*center);cb=(u-v)/(4*center);
        }
        double l=2*(edges[j]-a)/width-1,h=2*(edges[j+1]-a)/width-1;
        frac=((h-l)+ca*(h*h-l*l)/2+cb*(h*h*h-l*l*l)/3)/(2*(1+cb/3));
      } else if(degree==-2) {
        frac=erfBins[j]/erfNorm;
      } else if(degree<0) {
        const double slope=p[base], l=(edges[j]-a)/width, h=(edges[j+1]-a)/width;
        frac=(std::abs(slope)<1.e-7)?h-l:std::exp(slope*l)*std::expm1(slope*(h-l))/std::expm1(slope);
      } else {
        for(int i=0;i<=degree;++i)frac+=coefficients[i]*basis[j][i];
        frac*=(degree+1)/sum;
      }
      out[0][j]=std::exp(p[0])*frac;
    }
    for(int k=0;k<2;++k) {
      Shape s=shape[k];
      const int ti=floatWidth?7:5;
      const int offset=independentAlpha?2*k:(tailCount==8?4*k:0);
      if(tailCount>=1)s.al*=std::exp(p[ti+offset]);
      if(tailCount>=2)s.ar*=std::exp(p[ti+offset+1]);
      if(tailCount>=4&&!independentAlpha){
        if(fractionalN){s.nl*=std::exp(p[ti+offset+2]);s.nr*=std::exp(p[ti+offset+3]);}
        else{s.nl=1+(s.nl-1)*std::exp(p[ti+offset+2]);s.nr=1+(s.nr-1)*std::exp(p[ti+offset+3]);}
      }
      const double mu=shape[k].mean+shape[k].sigma*p[3+k];
      const double sigma=shape[k].sigma*(floatWidth?std::exp(p[5+k]):1.);
      const double low=cdf(edges.front(),mu,sigma,s);
      const double norm=cdf(edges.back(),mu,sigma,s)-low;
      if(!(norm>0))return {};
      for(int j=0;j<nb;++j)out[k+1][j]=std::exp(p[k+1])*(cdf(edges[j+1],mu,sigma,s)-cdf(edges[j],mu,sigma,s))/norm;
    }
    return out;
  }
  double evaluate(const double*p) const {
    auto c=components(p); if(c.empty())return 1.e100;
    double result=0;
    for(size_t i=0;i<counts.size();++i) {
      const double m=c[0][i]+c[1][i]+c[2][i], n=counts[i];
      if(!(m>0)||!std::isfinite(m))return 1.e100;
      result+=m-n+(n>0?n*std::log(n/m):0);
    }
    return result;
  }
};
}

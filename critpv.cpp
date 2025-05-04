#include <math.h>
#include <fstream>
#include <string.h>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <random> 
#include "distributions.cpp"

using namespace std;

double crit_signe(double stat);

//////////////////////////////////////////////////////////

main () {
 int i,k;
 long int mm,j;
 string st,ff;
 double *cx1,*cx2,*lgN,*x1,*x2,*x3,*x4,*nsigma,alfa1,lgLG1,alfa2,lgLG2,p,wp;
 long int *nplus,*nminus,*nnull;

 ff="critpv";
 ifstream fopen(ff+".inp");
 ofstream fout(ff+".out");

 fopen>>st;
 fopen>>mm;
 fopen>>st;
 fopen>>k;

 cx1=new double[k];
 cx2=new double[k];
 lgN=new double[k];
 x1=new double[k];
 x2=new double[k];
 x3=new double[k];
 x4=new double[k];
 nplus=new long int[k];
 nnull=new long int[k];
 nminus=new long int[k];
 nsigma=new double[k];

for (i=0;i<k;i++) {nplus[i]=0;nnull[i]=0;nminus[i]=0;}

 fopen>>st;
 for (i=0;i<k;i++) fopen>>lgN[i];
 fopen>>st;
 for (i=0;i<k;i++) fopen>>nsigma[i];
 fopen>>st;
 for (i=0;i<k;i++) fopen>>cx1[i];
 fopen>>st;
 fopen>>alfa1;
 fopen>>st;
 fopen>>lgLG1;
 fopen>>st;
 for (i=0;i<k;i++) fopen>>cx2[i];
 fopen>>st;
 fopen>>alfa2;
 fopen>>st;
 fopen>>lgLG2;
 fopen.close();

 double c001,c005,c01;
 long int pcrit;
 c001=1.2879;c005=0.98;c01=0.8224;
 //pcrit=((mm-1.)/2.-c005*sqrt(mm+1.));
 pcrit=normal_inv(0.05,0.5*mm,0.5*sqrt(mm));
 fout<<"pcrit="<<pcrit<<endl;  
 
 random_device rd;
 mt19937 g(rd());
 uniform_real_distribution<> dist(0,1);
 
 for(j=0;j<mm;j++) {
    for(i=0;i<k;i++) {
      p=dist(g);
      wp=log10(log(1./(1.-p)))-log10(log(2.));
      x1[i]=(0.5*cx1[i]/alfa1)*(1.+pow(10,nsigma[i]*(1.923-lgLG1+wp))); 
    }
    for(i=0;i<k;i++) x3[i]=x1[i]/x1[2];
    for(i=0;i<k;i++) {
        p=dist(g);
        wp=log10(log(1./(1.-p)))-log10(log(2.));
        x2[i]=(0.5*cx2[i]/alfa2)*(1.+pow(10,nsigma[i]*(1.923-lgLG2+wp))); 
    }
    for(i=0;i<k;i++) x4[i]=x2[i]/x2[2];
    for(i=0;i<k;i++) {
      if(x3[i]>x4[i]) nplus[i]=nplus[i]+1;
      if(x3[i]==x4[i]) nnull[i]=nnull[i]+1;
      if(x3[i]<x4[i])  nminus[i]=nminus[i]+1;
    }
 }
 
 fout<<"nplus=";    
 for(i=0;i<k;i++) fout<<nplus[i]<<"    ";
 fout<<endl;
 fout<<"nnull=";
 for(i=0;i<k;i++) fout<<nnull[i]<<"    ";
 fout<<endl;
 fout<<"nminus=";
 for(i=0;i<k;i++) fout<<nminus[i]<<"    ";
 fout<<endl;

 for(i=0;i<k;i++) {
    if (fmin(nplus[i],nminus[i])>=pcrit) fout<<"H0+"<<"   ";
    if (fmin(nplus[i],nminus[i])<pcrit) fout<<"H0-"<<"   ";
 }


 fout.close();
 delete [] cx1,cx2,lgN,x1,x2,x3,x4,nsigma,nplus,nminus,nnull;
 return 0;
}

////////////////////////////////////////////////////////////////////////////////

double crit_signe(double stat) {
   long int i,n;
   double *pw,*w;
   double pvalue;

   n=100;//int(stat);
   w=new double[n+1];
   pw=new double[n+1];
   for(i=0;i<=n;i++) {
      w[i]=double(i);
      pw[i]=binomial_cdf(w[i],n,0.5);
      if(w[i]>=stat) {pvalue=pw[i];break;}
   }

   delete [] w,pw;
   return pvalue;
}

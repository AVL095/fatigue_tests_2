#include <cmath>
#include <random>

using namespace std;

double pi=M_PI;

double lowRegGamma(double a, double x);
double gammapinv(double p, double a);
double betaln(double x, double y);
double betacf(double x, double a, double b);
double ibeta(double x, double a, double b);
double ibetainv(double p, double a, double b);
double beta_pdf(double x,double alpha, double beta);
double beta_cdf(double x,double alpha, double beta);
double beta_inv(double x,double alpha,double beta);
double betafn(double x, double y);

double normal_pdf(double x,double mu,double sigma);
double normal_pdf(double x,double mu,double sigma);
double normal_inv(double p, double mu, double sigma);
double studentt_pdf(double x,double dof);
double studentt_cdf(double x,double dof);
double studentt_inv(double p,double dof);
double chisquare_pdf(double x,double dof);
double chisquare_cdf(double x,double dof);
double chisquare_inv(double p,double dof);
double centralF_pdf(double x,double df1,double df2);
double centralF_cdf(double x,double df1,double df2);
double centralF_inv(double x, double df1,double df2);
double  poisson_pdf(long k,long l);
double poisson_cdf(double x,long l);
   

double binomial_pdf(double k, double n, double p);
double binomial_cdf(double x, int n, double p);
double factorialln(double n);
double factorial(double n);
double combination(double n, double m);
double combinationln(double n, double m);
double permutation(double n, double m);
double betinc(double x,double a,double b,double eps);

double invnontap(double,double,double);
double rtwi(double (*fct)(double), double xst, double eps, int iend);
double newton(double fct(double),double dfdx(double),double x, double eps,int &lim);
double noncentralt_cdf(double x, double dof, double ncp);
double noncentralt_pdf(double x, double dof, double ncp);
double noncentralt_inv(double beta,double dof,double ncp);
double lstat_invnontapp(double beta,double f,double d);
double lstat_funtquantile(double x);
double lstat_jennet_welch(double beta,double f,double d);
double lstat_f(double x);
double lstat_dfdx(double x);
void standart(double *x,int k,double p,double &median,double &mean,double &var,double &st,double &cvar,double &xp);

struct nct_tt {
 double beta;
 double f;
 double delta;
};
nct_tt nctt;



//###################################################################

double betacf(double x, double a, double b) {
    const double fpmin = 1e-30;
    int m = 1;
    double qab = a + b;
    double qap = a + 1;
    double qam = a - 1;
    double c = 1.;
    double d = 1.- qab * x / qap;
    double m2, aa, del, h;

    // These q's will be used in factors that occur in the coefficients
    if (abs(d) < fpmin)   d = fpmin;
    d = 1./d;
    h = d;

    for (; m <= 100; m++) {
        m2 = 2.*m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        // One step (the even one) of the recurrence
        d = 1. + aa * d;
        if (abs(d) < fpmin)
            d = fpmin;
        c = 1. + aa / c;
        if (abs(c) < fpmin)
            c = fpmin;
        d = 1./d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        // Next step of the recurrence (the odd one)
        d = 1.+aa * d;
        if (abs(d) < fpmin)
            d = fpmin;
        c = 1.+ aa / c;
        if (abs(c) < fpmin)  c = fpmin;
        d = 1./d;
        del = d * c;
        h *= del;
        if (abs(del - 1.0) < 3e-7)
            break;
    }

    return h;
}

//####################################################################

// Returns the incomplete beta function I_x(a,b)
double ibeta(double x, double a, double b) {
    double bt = (x == 0 || x == 1) ? 0 :  exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.-x));
    if (x < 0 || x > 1)     return false;
    if (x < (a + 1.) / (a+b+2.))  return bt * betacf(x, a, b) / a;
    return 1.-bt*betacf(1.-x,b,a)/b;
}

//#########################################################################

// Returns the inverse of the incomplete beta function
double ibetainv(double p, double a, double b) {
    const double EPS = 1e-8;
    double a1 = a - 1.;
    double b1 = b - 1.;
    int j = 0;
    double lna, lnb, pp, t, u, err, x, al, h, w, afac;

    if (p <= 0.0)
        return 0.0;
    if (p >= 1.0)
        return 1.0;

    if (a >= 1.0 && b >= 1.0) {
        pp = (p < 0.5) ? p : 1.0 - p;
        t = sqrt(-2.0 * log(pp));
        x = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
        if (p < 0.5)
            x = -x;
        al = (x * x - 3.0) / 6.0;
        h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
        w = (x * sqrt(al + h) / h) - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (al + 5.0 / 6.0 - 2.0 / (3.0 * h));
        x = a / (a + b * exp(2.0 * w));
    } else {
        lna = log(a / (a + b));
        lnb = log(b / (a + b));
        t = exp(a * lna) / a;
        u = exp(b * lnb) / b;
        w = t + u;
        if (p < t / w)
            x = pow(a * w * p, 1.0 / a);
        else
            x = 1.0 - pow(b * w * (1.0 - p), 1.0 / b);
    }

    afac = -lgamma(a) - lgamma(b) + lgamma(a + b);

    for (; j < 10; j++) {
        if (x == 0.0 || x == 1.0)
            return x;
        err = ibeta(x, a, b) - p;
        t = exp(a1 * log(x) + b1 * log(1.0 - x) + afac);
        u = err / t;
        x -= (t = u / (1.0 - 0.5 * min(1.0, u * (a1 / x - b1 / (1.0 - x)))));
        if (x <= 0.0)
            x = 0.5 * (x + t);
        if (x >= 1.0)
            x = 0.5 * (x + t + 1.0);
        if (abs(t) < EPS * x && j > 0)
            break;
    }

    return x;
}

//#########################################################

double betafn(double x, double y) {
  if (x <= 0 || y <= 0) return -1;
  return (x + y > 170) ? exp(betaln(x, y)) : tgamma(x) * tgamma(y) /tgamma(x + y);
}


//####################################################################

// extend beta function with static methods
 double beta_pdf(double x,double alpha, double beta) {
    // PDF is zero outside the support
    if (x > 1 || x < 0)   return 0;
    // PDF is one for the uniform case
    if (alpha == 1 && beta == 1)    return 1;
    if (alpha < 512 && beta < 512) {
      return (pow(x, alpha - 1)*pow(1 - x, beta - 1)) /betafn(alpha, beta);
    } else {
      return exp((alpha - 1)*log(x) + (beta - 1)*log(1 - x) -betaln(alpha, beta));
    }
  }

//########################################################################

double beta_cdf(double x,double alpha, double beta) {
    return (x > 1 || x < 0) ? (x > 1) * 1 : ibeta(x, alpha, beta);
  }

//#######################################################################

  double beta_inv(double x,double alpha,double beta) {
    return ibetainv(x, alpha, beta);
  }

//#######################################################################

// natural logarithm of beta function
double betaln(double x,double y) {
  return lgamma(x) + lgamma(y) - lgamma(x + y);
}

//##################################################################

// The lower regularized incomplete gamma function, usually written P(a,x)
double lowRegGamma(double a, double x) {
    double aln = lgamma(a);
    double ap = a;
    double sum = 1 / a;
    double del = sum;
    double b = x + 1 - a;
    double c = 1 / 1.0e-30;
    double d = 1 / b;
    double h = d;
    int i = 1;
    // calculate maximum number of iterations required for a
    int ITMAX = -(log((a >= 1) ? a : 1 / a) * 8.5 + a * 0.4 + 17);
    double an;

    if (x < 0 || a <= 0) {
        return numeric_limits<double>::quiet_NaN();
    } else if (x < a + 1) {
        for (; i <= ITMAX; i++) {
            sum += del *= x / ++ap;
        }
        return (sum * exp(-x + a * log(x) - (aln)));
    }

    for (; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2;
        d = an * d + b;
        c = b + an / c;
        d = 1 / d;
        h *= d * c;
    }
    return (1 - h * exp(-x + a * log(x) - (aln)));
}

//##################################################################

// Returns the inverse of the lower regularized incomplete gamma function
double gammapinv(double p, double a) {
    int j = 0;
    double a1 = a - 1;
    double EPS = 1e-8;
    double gln = lgamma(a);
    double x, err, t, u, pp, lna1, afac;

    if (p >= 1)
        return max(100.0, a + 100.0 * sqrt(a));
    if (p <= 0)
        return 0.0;
    if (a > 1) {
        lna1 = log(a1);
        afac = exp(a1 * (lna1 - 1) - gln);
        pp = (p < 0.5) ? p : 1 - p;
        t = sqrt(-2 * log(pp));
        x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
        if (p < 0.5)
            x = -x;
        x = max(1e-3, a * pow(1 - 1 / (9 * a) - x / (3 * sqrt(a)), 3));
    } else {
        t = 1. - a * (0.253 + a * 0.12);
        if (p < t)
            x = pow(p / t, 1./ a);
        else
            x = 1. - log(1. - (p - t) / (1. - t));
    }

    for (; j < 12; j++) {
        if (x <= 0)
            return 0.0;
        err = lowRegGamma(a, x) - p;
        if (a > 1)
            t = afac * exp(-(x - a1) + a1 * (log(x) - lna1));
        else
            t = exp(-x + a1 * log(x) - gln);
        u = err / t;
        x -= (t = u / (1. - 0.5 * min(1.0, u * ((a - 1.) / x - 1.))));
        if (x <= 0)
            x = 0.5 * (x + t);
        if (abs(t) < EPS * x)
            break;
    }

    return x;
}


//##########################################################

double normal_pdf(double x,double mu,double sigma) {
    double diff;
    diff=(x-mu)/sigma;
    return(exp(-0.5*diff * diff) /(sigma*sqrt(2.*pi)));
}
//############################################

double normal_cdf(double x,double mu,double sigma) {
    return(0.5*(1.0+erf((x-mu)/(sigma *sqrt(2.0)))));
}

//###############################################################

double normal_inv(double p, double mu, double sigma) {

    double q,r,num,x,den;
  
    q = p - 0.5;
    if(fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        num = (((((((2.5090809287301226727e+3 * r +3.3430575583588128105e+4) * r + 6.7265770927008700853e+4) * r + 4.5921953931549871457e+4) * r +
                     1.3731693765509461125e+4) * r +1.9715909503065514427e+3) * r +1.3314166789178437745e+2) * r +3.3871328727963666080e+0) * q;
        den = (((((((5.2264952788528545610e+3 * r +2.8729085735721942674e+4) * r + 3.9307895800092710610e+4) * r + 2.1213794301586595867e+4) * r +
                     5.3941960214247511077e+3) * r + 6.8718700749205790830e+2) * r + 4.2313330701600911252e+1) * r +1.0);
        x = num / den;
        return(mu + (x * sigma));
     }

    if(q <= 0.0)  {
      r=p;
   }
    else {
      r=1.0-p;
    }
    r =sqrt(-log(r));
    if(r <= 5.0) {
        r = r - 1.6;
        num = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r +2.41780725177450611770e-1) * r +1.27045825245236838258e+0) * r +
                     3.64784832476320460504e+0) * r + 5.76949722146069140550e+0) * r +4.63033784615654529590e+0) * r +1.42343711074968357734e+0);
        den = (((((((1.05075007164441684324e-9 * r +5.47593808499534494600e-4) * r +1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r +
                     6.89767334985100004550e-1) * r +1.67638483018380384940e+0) * r +2.05319162663775882187e+0) * r +
                     1.0);
    } else {
        r = r - 5.0;
        num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r +
                     2.96560571828504891230e-1) * r +1.78482653991729133580e+0) * r + 5.46378491116411436990e+0) * r + 6.65790464350110377720e+0);
        den = (((((((2.04426310338993978564e-15 * r +1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +7.86869131145613259100e-4) * r +
                     1.48753612908506148525e-2) * r +1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r +1.0);
    }
    x = num / den;
    if(q < 0.0)  x = -x;
    return(mu + (x * sigma));

}


//#################################################################


 double studentt_pdf(double x,double dof) {
    dof = dof > 1e100 ? 1e100 : dof;
    return (1./(sqrt(dof) * betafn(0.5, dof/2))) * pow(1. + ((x * x) / dof), -((dof + 1.)/2.));
  }

double studentt_cdf(double x,double dof) {
    double dof2 = dof / 2;
    return ibeta((x + sqrt(x * x + dof)) /(2. * sqrt(x * x + dof)), dof2, dof2);
  }

 double studentt_inv(double p,double dof) {
   double x = ibetainv(2.*min(p, 1. - p), 0.5 * dof, 0.5);
    x = sqrt(dof * (1. - x) / x);
    return (p > 0.5) ? x : -x;
  }


//###########################################################################

// extend chisquare function with static methods

 double chisquare_pdf(double x,double dof) {
    if (x < 0)   return 0;
    return (x == 0 && dof == 2) ? 0.5 : exp((dof / 2. - 1.) * log(x) - x / 2 - (dof / 2.) * log(2.) - lgamma(dof / 2.));
  }

double chisquare_cdf(double x,double dof) {
    if (x < 0)   return 0;
    return lowRegGamma(dof / 2., x / 2.);
  }

double chisquare_inv(double p,double dof) {
    return 2*gammapinv(p, 0.5 * dof);
  }

//#################F-Distribution########################################################

// extend F function with static methods
  // This implementation of the pdf function avoids float overflow
  // See the way that R calculates this value:
  // https://svn.r-project.org/R/trunk/src/nmath/df.c
  
double centralF_pdf(double x,double df1,double df2) {
    double p, q, f;

    if (x < 0)   return 0;

    if (df1 <= 2) {
      if (x == 0 && df1 < 2)  return INFINITY;
      if (x == 0 && df1 == 2)   return 1;
      return (1./betafn(df1 /2.,df2/2.))*pow(df1/df2,df1/ 2.)*pow(x, (df1/2.)-1.)*pow((1.+(df1/df2)*x), -(df1+df2) / 2.);
    }

    p = (df1 * x) / (df2 + x * df1);
    q = df2 / (df2 + x * df1);
    f = df1 * q / 2.;
    return(f *binomial_pdf((df1 - 2.)/2.,(df1 + df2-2.)/2., p));
  }

//################################################

 double centralF_cdf(double x,double df1,double df2) {
    if (x < 0)   return 0;
    return ibeta((df1 * x) / (df1 * x + df2), df1 / 2., df2 / 2.);
  }

//#################################################

 double centralF_inv(double x, double df1,double df2) {
    return df2 / (df1 * (1. /ibetainv(x, df1 / 2., df2 / 2.) - 1.));
  }


//###################Binomial,Factorial###################################

double binomial_pdf(double k, double n, double p) {
    return (p == 0 || p == 1) ?   ((n * p) == k ? 1 : 0) :
     combination(n, k) * pow(p, k) * pow(1. - p, n - k);
  }

//########################################################################

double binomial_cdf(double x, int n, double p) {
    double betacdf;
    double eps = 1e-10;

    if (x < 0)    return 0.0;
    if (x >= n)    return 1.0;
    if (p < 0.0 || p > 1.0 || n <= 0)  return numeric_limits<double>::quiet_NaN();

    x = floor(x);
    double z = p;
    double a = x + 1.;
    double b = n - x;
    double s = a + b;
    double bt = exp(lgamma(s) - lgamma(b) - lgamma(a) + a * log(z) + b * log(1 - z));
    if (z < (a + 1.) / (s + 2.))
        betacdf = bt * betinc(z, a, b, eps);
    else
        betacdf = 1.0 - bt * betinc(1.0 - z, b, a, eps);
    return round((1.0 - betacdf) * (1.0 / eps)) / (1.0 / eps);
}

//###############################################

double factorialln(double n) {
  if(n<0) return(-1);
  return lgamma(n + 1);
}

//##################################################

// factorial of n
double factorial(double n) {
  if(n<0) return(-1);
  return tgamma(n + 1);
}

//########################################################

// combinations of n, m
double combination(double n, double m) {
  // make sure n or m don't exceed the upper limit of usable values
  return (n > 170 || m > 170) ? exp(combinationln(n, m)) : (factorial(n) / factorial(m)) / factorial(n - m);
}

//########################################################

double combinationln(double n, double m){
  return factorialln(n) -factorialln(m) - factorialln(n - m);
}

//#########################################################

// permutations of n, m
double permutation(double n, double m) {
  return factorial(n) / factorial(n - m);
}

//###############################################################

// Got this from http://www.math.ucla.edu/~tom/distributions/binomial.html
double betinc(double x,double a,double b,double eps) {
  double a0,b0,a1,b1,m9,a2,c9;
  a0 = 0; b0 = 1; a1 = 1; b1 = 1; m9 = 0; a2 = 0;

  while (abs((a1 - a2) / a1) > eps) {
    a2 = a1;
    c9 = -(a + m9) * (a + b + m9) * x / (a + 2. * m9) / (a + 2. * m9 + 1.);
    a0 = a1 + c9 * a0;
    b0 = b1 + c9 * b0;
    m9 = m9 + 1.;
    c9 = m9 * (b - m9) * x / (a + 2. * m9 - 1) / (a + 2. * m9);
    a1 = a0 + c9 * a1;
    b1 = b0 + c9 * b1;
    a0 = a0 / b1;
    b0 = b0 / b1;
    a1 = a1 / b1;
    b1 = 1.;
  }

  return a1 / a;
}

//################################################################


double noncentralt_pdf(double x, double dof, double ncp) {
    double tol = 1e-14;
    if (abs(ncp) < tol) return studentt_pdf(x, dof);   //ncp approx 0; use student-t
    if (abs(x) < tol) {  // different formula for x == 0
        return exp(lgamma((dof+1.)/2.)-ncp*ncp/2.-0.5*log(M_PI*dof)-lgamma(dof/2.));
    }
    // formula for x != 0
    return (dof/x)*(noncentralt_cdf(x*sqrt(1.+2./dof),dof+2.,ncp)-noncentralt_cdf(x,dof,ncp));
}

//##################################################################

double noncentralt_cdf(double x, double dof, double ncp) {
    double tol,value,prob,lastvalue,y,p,q,a,z,sw; 
    int min_iterations;
    int j;

    tol=1e-14;
    min_iterations = 50;
    if(dof>=15) min_iterations = 100;

   if (dof>5000) return normal_cdf(x,0,1);


   if (dof>80)  {
     a=sqrt(2./dof)*exp(lgamma(0.5*(dof+1.0))-lgamma(0.5*dof));
     sw=1.-a*a; 
     z=abs(ncp-x*a)/sqrt(1.+x*x*sw);
     prob=normal_cdf(z,0,1); 
     return prob;
   }
      
    // turn negative x into positive and flip result afterwards
    bool flip = false;
    if (x < 0) {
        flip = true; ncp = -ncp;
    }


  if(fabs(ncp) < tol) return studentt_cdf(x, dof);    // ncp approx 0; use student-t
  if   (fabs(ncp / (4. *dof)) < tol)    {
               prob =studentt_cdf(x - ncp,dof);return flip ? 1 - prob : prob;
            }


    prob = normal_cdf(-ncp, 0, 1);
    value=tol+1;
    // use value at last two steps to determine convergence
    lastvalue = value;
    y=x*x/(x*x+dof);
    j = 0;
    p = exp(-ncp*ncp/2.);
    q = exp(-ncp*ncp/2.-0.5*log(2.)-lgamma(3./2.)) * ncp;
    
    while (j < min_iterations || lastvalue > tol || value > tol) {
        lastvalue = value;
        if (j > 0) {
            p *= (ncp*ncp)/(2.*j);
            q *= (ncp*ncp)/(2. *(j+0.5));
            //q*=(ncp*ncp)/(2.*tgamma(j+1.5));
        }
        value =p*beta_cdf(y, j+0.5,dof/2.)+q*beta_cdf(y, j+1., dof/2.);
        prob+=0.5*value;
        j++;
    }

    return flip ? (1 - prob) : prob;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double noncentralt_inv_1(double beta, double df, double ncp) {

     double mid,lower_bound,upper_bound,p,eps;
     int i,m;

     eps=1e-15;
     m=100;
     lower_bound = -m;
     upper_bound = m;

    for(i = 0; i < m; ++i) {
        mid= (lower_bound + upper_bound) / 2;
        p=noncentralt_cdf(mid,df,ncp);
        if (p<beta) {
            lower_bound = mid;
        } else {
            upper_bound = mid;
        }
        if(abs(p-beta)<eps) break;
    }
    return (lower_bound + upper_bound) / 2;
}



//###########################################################

double noncentralt_inv(double beta,double dof,double ncp) {

     double eps,xst,x;
     int iend;
     eps=1e-15;iend=15;
    
     if(dof<10) {
       xst=ncp+ normal_inv(beta,0,1);
      }
      else {
       //xst=lstat_invnontapp(beta,dof,ncp);
       xst=invnontap(beta,dof,ncp);
       //xst=lstat_jennet_welch(beta,dof,ncp);
      }

     nctt.beta=beta;nctt.f=dof;nctt.delta=ncp;
     x=rtwi(lstat_funtquantile,xst,eps,iend);
     //x=newton(lstat_f,lstat_dfdx,xst,eps,iend);

     return x;

}


/////////////////////////////////////////////////////////////////////////////////////////
/*
    PURPOSE TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM X=FCT(X) BY MEANS OF WEGSTEIN-S ITERATION METHOD

            USAGE
           double rtwi(double (*fct)(double), double xst, double eps, int iend)

        DESCRIPTION OF PARAMETERS
            x     - RESULTANT ROOT OF EQUATION x=fct(x).
            fct    - NAME OF THE EXTERNAL FUNCTION SUBPROGRAM USED.
            xst    - INPUT VALUE WHICH SPECIFIES THE INITIAL GUESS OF THE ROOT x.
            eps    - INPUT VALUE WHICH SPECIFIES THE UPPER BOUND OF THE ERROR OF RESULT x.
            iend   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.
            ier=0  - NO ERROR,
            ier=1  - NO CONVERGENCE AFTER IEND ITERATION STEPS,
            ier=2  - AT ANY ITERATION STEP THE DENOMINATOR OF ITERATION FORMULA WAS EQUAL TO ZERO.

*/

double rtwi(double (*fct)(double), double xst, double eps, int iend) {
    int i,ier;
    double x,tol,a,b,d,val;

    ier = 0; tol = xst; x = fct(tol);a= x - xst;b = -a;tol = x;val = x - fct(tol);

  for (i = 0; i < iend; i++) {
        if(val==0) return x;
        b = b / val - 1.;
         if (b == 0) {
             ier = 2; return x;
         }
         a = a / b; x = x + a; b = val;tol = x;val = x - fct(tol); tol = eps; d = abs(x);
         if (d > 1.) tol *= d;
         if (abs(a) > tol) continue;
         if (abs(val) >10.*tol) continue;
         return x;
   }

    ier = 1;
    return x;
}

/////////////////////////////////////////////////////////////////////
/*

            PURPOSE
           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM fct(x)=0 BY MEANS OF NEWTON-S ITERATION METHOD.

        USAGE
           double newton(double fct(double),double dfdx(double),double x, double eps,int &lim)
           PARAMETER fct REQUIRES AN EXTERNAL STATEMENT.

        DESCRIPTION OF PARAMETERS
           x      - RESULTANT ROOT OF EQUATION fct(x)=0 and INPUT VALUE WHICH SPECIFIES THE INITIAL GUESS OF THE ROOT x.
           dfdx   - RESULTANT VALUE OF DERIVATIVE AT ROOT x.
           fct    - NAME OF THE EXTERNAL SUBROUTINE USED. IT COMPUTES
                    TO GIVEN ARGUMENT x FUNCTION VALUE fct AND DERIVATIVE
                    dfdx. ITS PARAMETER LIST MUST BE x,fct,dfdx.
           eps    - INPUT VALUE WHICH SPECIFIES THE UPPER BOUND OF THE ERROR OF RESULT x.
           lim   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.
           ier    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS (no output),
           0 - NO ERROR,
           1 - NO CONVERGENCE AFTER IEND ITERATION STEPS,
           2 - AT ANY ITERATION STEP DERIVATIVE dfdx WAS EQUAL TO ZERO.

        REMARKS
           THE PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE ier=2
           IF AT ANY ITERATION STEP DERIVATIVE OF fct(x) IS EQUAL TO 0.
           POSSIBLY THE PROCEDURE WOULD BE SUCCESSFUL IF IT IS STARTED
           ONCE MORE WITH ANOTHER INITIAL GUESS x.

        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
           THE EXTERNAL SUBROUTINE Fct(X) and dfdx(x) MUST BE FURNISHED BY THE USER.

        METHOD
           SOLUTION OF EQUATION F(X)=0 IS DONE BY MEANS OF NEWTON-S
           ITERATION METHOD, WHICH STARTS AT THE INITIAL GUESS XST OF
           A ROOT X. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF
           F(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP
           REQUIRES ONE EVALUATION OF F(X) AND ONE EVALUATION OF THE
           DERIVATIVE OF F(X). FOR TEST ON SATISFACTORY ACCURACY SEE
           FORMULAE (2) OF MATHEMATICAL DESCRIPTION.
           FOR REFERENCE, SEE R. ZURMUEHL, PRAKTISCHE MATHEMATIK FUER
           INGENIEURE UND PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/
           HEIDELBERG, 1963, PP.12-17.
*/
double newton(double fct(double),double dfdx(double),double x, double eps,int &lim) {

    double f,df;
    int iter;

  iter=0;
  do   {
    f=fct(x);
    df=dfdx(x);
    x=x-f/df;
    iter++;
  } while (fabs(f)>eps && iter<lim);
    lim=iter;
    return x;
}


////////////////////////////////////////////////////////

double lstat_funtquantile(double x) {
    return(x-nctt.beta+noncentralt_cdf(x,nctt.f,nctt.delta));
}
///////////////////////////////////////////////////////////

double lstat_f(double x) {
    return(noncentralt_cdf(x,nctt.f,nctt.delta)-nctt.beta);
}

////////////////////////////////////////////////////////

double lstat_dfdx(double x) {
  return noncentralt_pdf(x,nctt.f,nctt.delta);
}


//*****************Приближенные доверительные граница для квантиля**************************
void lmtaprn(double beta,int f,double t1,double t2,double t12,double d,double &tlow,double &tup){

double zb,f1x,f2x,f4x,e11,e2x,e3x,e44;

zb=normal_inv(beta,0,1);
f1x=t2/f;
f2x=2*t12/sqrt(f+1.);
f4x=1.0-f1x/2;
e3x=f4x*f4x-zb*zb*f1x;
e11=f4x*d+zb*zb*f2x/2;
e2x=d*d-zb*zb*t1;
e44=sqrt(abs(e11*e11-e2x*e3x));
tlow=(e11-e44)/e3x;
tup=(e11+e44)/e3x;
}

//************************************************************************
double invnontap(double beta,double f,double d) {
 double zb, z, f4x;

 zb=normal_inv(beta,0,1);
 f4x=1.0-1.0/(4.0*f);
 z=(f4x*d+zb*sqrt(f4x*f4x-zb*zb/(2.*f)+d*d/(2.*f)))/(f4x*f4x-zb*zb/(2.*f));
 return z;
}

/////////////////////////////////////////////////////////////

double lstat_invnontapp(double beta,double f,double d) {

 double z,zb,zx,f4x;
 zb=normal_inv(beta,0,1);
 f4x=1.-1./(4.*f);
 z=f4x*f4x-zb*zb/(2.*f);
 if(z<0) {
   z=1.;f4x=1.;
 }
 zx=(f4x*d+zb*sqrt(z+d*d/(2.*f)))/z;
 return(zx);
}

//########################Poisson Distribution########################################

// extend uniform function with static methods
double  poisson_pdf(long k,long l) {
    if (l < 0 || (k % 1) != 0 || k < 0) return 0;
    return pow(l, k)*exp(-l)/factorial(k);
  }
//###################################################
 double poisson_cdf(double x,long l) {
    long k;
    double sum;
    sum=0;
    if (x < 0) return 0;
    for (k=0; k <= x; k++) {
      sum+=poisson_pdf(k, l);
    }
    return sum;
  }

//###################################################

void standart(double *x,int k,double p,double &median,double &mean,double &var,double &st,double &cvar,double &xp) {
  double s1,s2,alphap;
  int i;
 
  s1=0;s2=0; 
  for(i= 0;i<k;i++) {
    s1+= x[i];s2+= x[i]*x[i];
  }        
  mean=s1/k;
  median=!(k & 1) ? (x[(k/2)-1]+x[(k/2)])/2 : x[(k/2) | 0];
  var=(s2-mean*mean*k)/(k-1.);
  st=sqrt(var);
  cvar=st/mean; 

  i=trunc(p*(k+1.));
  if(i<=1) {
    xp=x[0];
  }
  else {
   alphap=p*(k+1.)-i;
   xp=(1.-alphap)*x[i-1]+alphap*x[i];
  }
}


/////////////////////////////////////////////////////////////

double lstat_jennet_welch(double beta,double f,double d) {

 double b,zb,z,f4;

 zb=normal_inv(beta,0,1);
 b=sqrt(2./f)*tgamma(0.5*(f+1.0))/tgamma(0.5*f);

 f4=b*b-zb*zb*(1.-b*b);
 z=(b*d+zb*sqrt(b*b+(1.-b*b)*(d*d-zb*zb)))/f4;
 return(z);
}

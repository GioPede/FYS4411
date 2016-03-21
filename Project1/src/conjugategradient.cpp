#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vectormatrix.h"
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "sampler.h"



using namespace  std;

//double E_function(Vector &x);

void E_function(Vector &x, double &fp, Vector &g);

void lnsrch(int n, Vector &xold, double fold, Vector &g, Vector &p, Vector &x,
         double *f, double stpmax, int *check, void (*func)(Vector &p, double &fp, Vector &g));

void dfpmin(Vector &p, int n, double gtol, int *iter, double *fret,
        void (*func)(Vector &p,  double &fp, Vector &g));

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
/*
//  this function defines the Energy function
double E_function(Vector  &x)
{
    int numberOfDimensions  = 3;
    int numberOfParticles   = 3;
    int numberOfSteps       = (int) 1e4;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = x(0);          // Variational parameter (fixed or computed via dfpmin)
    double beta             = 1.0;          // Ellyptic factor
    double interactionRange = 1e-3;         // Typical distance of interaction
    double stepLength       = 0.001;        // Metropolis step length.
    double derivationStep   = 0.01;         // Derivation step length
    double equilibration    = 0.1;          // Amount of the total steps used
                                            // for equilibration.

//    conjugateGradient(numberOfParameters, minimumTolerance, minimumEstimate,
//                      numberOfIterations, parameterGuess,
//                      numberOfDimensions, numberOfParticles, numberOfSteps, omega,
//                      alpha, beta, interactionRange, stepLength,
//                      derivationStep, equilibration);

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, derivationStep, alpha, beta, interactionRange));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps, 0);

    return system->getSampler()->getCumulativeEnergy();
} // end of function to evaluate*/

//  this function defines the derivative of the energy
void E_function(Vector &x, double &fp, Vector &g)
{
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e4;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = x(0);          // Variational parameter (fixed or computed via dfpmin)
    double beta             = 1.0;//2.82843;          // Ellyptic factor
    double interactionRange = 0.0043;         // Typical distance of interaction
    double stepLength       = 0.01;        // Metropolis step length.
    double derivationStep   = 0.01;         // Derivation step length
    double equilibration    = 0.1;          // Amount of the total steps used
                                            // for equilibration.

//    conjugateGradient(numberOfParameters, minimumTolerance, minimumEstimate,
//                      numberOfIterations, parameterGuess,
//                      numberOfDimensions, numberOfParticles, numberOfSteps, omega,
//                      alpha, beta, interactionRange, stepLength,
//                      derivationStep, equilibration);

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, derivationStep, alpha, beta, interactionRange));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps, 0);
    fp = system->getSampler()->getCumulativeEnergy();
    g(0) = (2*(system->getSampler()->getCumulativeDePsideAlphaELocal()-
              system->getSampler()->getCumulativeDePsideAlpha()* system->getSampler()->getCumulativeEnergy()));
    //std::cout << system->getSampler()->getCumulativeEnergy() << g(0) << std::endl;
} // end of function to evaluate




#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 1000.0


void dfpmin(Vector &p, int n, double gtol, int *iter, double *fret,
        void (*func)(Vector &p, double &fp, Vector &g))
{

  int check,i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  Vector dg(n), g(n), hdg(n), pnew(n), xi(n);
  Matrix hessian(n,n);

  (*func)(p,fp, g);
  for (i = 0;i < n;i++) {
    for (j = 0; j< n;j++) hessian(i,j)=0.0;
    hessian(i,i)=1.0;
    xi(i) = -g(i);
    sum += p(i)*p(i);
  }
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=1;its<=ITMAX;its++) {

      std::cout<<its<<std::endl;
    *iter=its;
    lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
    fp = *fret;
    for (i = 0; i< n;i++) {
      xi(i)=pnew(i)-p(i);
      p(i)=pnew(i);
    }
    test=0.0;
    for (i = 0;i< n;i++) {
      temp=fabs(xi(i))/FMAX(fabs(p(i)),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i);
    double tempa;
    (*func)(p,tempa,g);
    test=0.0;
    den=FMAX(*fret,1.0);
    for (i=0;i<n;i++) {
      temp=fabs(g(i))*FMAX(fabs(p(i)),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i)-dg(i);
    for (i=0;i<n;i++) {
      hdg(i)=0.0;
      for (j=0;j<n;j++) hdg(i) += hessian(i,j)*dg(j);
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) {
      fac += dg(i)*xi(i);
      fae += dg(i)*hdg(i);
      sumdg += SQR(dg(i));
      sumxi += SQR(xi(i));
    }
    if (fac*fac > EPS*sumdg*sumxi) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=0;i<n;i++) dg(i)=fac*xi(i)-fad*hdg(i);
      for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      hessian(i,j) += fac*xi(i)*xi(j)
        -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
    }
      }
    }
    for (i=0;i<n;i++) {
      xi(i)=0.0;
      for (j=0;j<n;j++) xi(i) -= hessian(i,j)*g(j);
    }
  }
  cout << "too many iterations in dfpmin" << endl;
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n, Vector &xold, double fold, Vector &g, Vector &p, Vector &x,
        double *f, double stpmax, int *check, void (*func)(Vector &p, double &fp, Vector &g))
{
  int i;
  double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  *check=0;
  for (sum=0.0,i=0;i<n;i++) sum += p(i)*p(i);
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=0;i<n;i++) p(i) *= stpmax/sum;
  for (slope=0.0,i=0;i<n;i++)
    slope += g(i)*p(i);
  test=0.0;
  for (i=0;i<n;i++) {
    temp=fabs(p(i))/FMAX(fabs(xold(i)),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=0;i<n;i++) x(i)=xold(i)+alam*p(i);
    Vector tempa(n);
    (*func)(x, *f, tempa);
    if (alam < alamin) {
      for (i=0;i<n;i++) x(i)=xold(i);
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
    tmplam = -slope/(2.0*(*f-fold-slope));
      else {
    rhs1 = *f-fold-alam*slope;
    rhs2=f2-fold2-alam2*slope;
    a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
    b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
    if (a == 0.0) tmplam = -slope/(2.0*b);
    else {
      disc=b*b-3.0*a*slope;
      if (disc<0.0) cout << "Roundoff problem in lnsrch." << endl;
      else tmplam=(-b+sqrt(disc))/(3.0*a);
    }
    if (tmplam>0.5*alam)
      tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX



//   Main function begins here
void conjugateGradient()
{
     int n, iter;
     double gtol, fret;
     double alpha;
     n = 1;
//   reserve space in memory for vectors containing the variational
//   parameters
     Vector g(n), p(n);
     cout << "Read in guess for alpha" << endl;
     cin >> alpha;
     gtol = 1.0e-5;
//   now call dfmin and compute the minimum
     p(0) = alpha;
     dfpmin(p, n, gtol, &iter, &fret, E_function);
     cout << "Value of energy minimum = " << fret << endl;
     cout << "Number of iterations = " << iter << endl;
     cout << "Value of alpha at minimum = " << p(0) << endl;
}  // end of main program


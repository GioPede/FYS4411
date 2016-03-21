#include <iostream>
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
#include "conjugategradient.h"

using namespace std;

int main(int argn, char *argv[]) {

    //conjugateGradient();  // Uncomment for quick usage (very unstable)

    // SYSTEM RELATED VARIABLES
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;                      // Oscillator frequency.
    double alpha            = 0.5; //atof(argv[1]);     // Variational parameter (fixed or computed via dfpmin)
    double beta             = 2.82843;                  // Elliptic factor
    double interactionRange = 0.0043;                   // Typical distance of interaction
    double stepLength       = 0.01;                     // Metropolis step length.
    double derivationStep   = 0.01;                     // Derivation step length
    double equilibration    = 0.1;                      // Amount of the total steps used for equilibration.

    // INITIALIZE THE SYSTEM
    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, derivationStep, alpha, beta, interactionRange));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    // RUN SIMULATION
    system->runMetropolisSteps          (numberOfSteps, 1);  // bool for printing output to file (for blocking)

    // OPTIONAL TERMINAL OUTPUT
    system->getSampler()->printOutputToTerminal();
    return 0;
}

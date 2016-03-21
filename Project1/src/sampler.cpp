#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system, bool print) {
    m_system = system;
    m_stepNumber = 0;
    m_output = fopen("energy.txt", "w+");
    m_print = print;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

// Make sure the sampling variable(s) are initialized at the first step.
void Sampler::initializeFirstStep(){
    m_cumulativeEnergy = 0;
    m_cumulativeEnergySquared = 0;
    m_localEnergy = m_system->getHamiltonian()->
            computeLocalEnergy(m_system->getParticles());
    m_dePsideAlpha = 0.0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
            m_dePsideAlpha -= m_system->getParticles().at(i)->getPosition().at(j)*m_system->getParticles().at(i)->getPosition().at(j);
        }
        m_dePsideAlpha -= m_system->getParticles().at(i)->getPosition().back() * m_system->getParticles().at(i)->getPosition().back() *
                m_system->getWaveFunction()->getParameters().at(1);
    }
    m_dePsideAlphaELocal = m_dePsideAlpha * m_localEnergy;
}

// PERFORM THE ACTUAL SAMPLING
void Sampler::sample(bool acceptedStep) {
    if(acceptedStep){

        m_localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
        m_dePsideAlpha = 0.0;
        for (int i=0; i < m_system->getNumberOfParticles(); i++) {
            for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
              m_dePsideAlpha -= m_system->getParticles().at(i)->getPosition().at(j)
                               *m_system->getParticles().at(i)->getPosition().at(j);
            }
            m_dePsideAlpha -= m_system->getParticles().at(i)->getPosition().back()
                             *m_system->getParticles().at(i)->getPosition().back()
                             *m_system->getWaveFunction()->getParameters().at(1);
        }
        m_dePsideAlphaELocal = m_dePsideAlpha * m_localEnergy;
        m_numberAccepted++;
    }
    if(m_print) fprintf(m_output, "%.15lf\n", m_localEnergy);
    m_cumulativeEnergy  += m_localEnergy;
    m_cumulativeDePsideAlpha  += m_dePsideAlpha;
    m_cumulativeDePsideAlphaELocal  += m_dePsideAlphaELocal;
    m_cumulativeEnergySquared  += m_localEnergy*m_localEnergy;

    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " << m_variance << endl;
    cout << " Acceptance Ratio  : " << m_numberAccepted/m_stepNumber<< endl;
    cout << endl;

    fclose(m_output);
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */

    m_energy = m_cumulativeEnergy / m_stepNumber ;/// m_system->getNumberOfParticles();
    m_variance = m_cumulativeEnergySquared / m_stepNumber - m_energy*m_energy;/// m_system->getNumberOfParticles() / m_system->getNumberOfParticles() - m_energy*m_energy;

}

double Sampler::getCumulativeEnergy() const
{
    return m_cumulativeEnergy / m_stepNumber;
}

double Sampler::getCumulativeDePsideAlpha() const
{
    return m_cumulativeDePsideAlpha / m_stepNumber;
}

double Sampler::getCumulativeDePsideAlphaELocal() const
{
    return m_cumulativeDePsideAlphaELocal / m_stepNumber;
}

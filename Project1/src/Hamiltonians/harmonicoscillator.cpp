#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

// COMPUTE THE LOCAL ENERGY ANALITICALLY
double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    double alpha = m_system->getWaveFunction()->getParameters().at(0);
    double beta = m_system->getWaveFunction()->getParameters().at(1);
    double interactionRange =  m_system->getWaveFunction()->getParameters().at(2);

    // compute the harmonic elliptic potential
    double r2 = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
          r2 += m_system->getParticles().at(i)->getPosition().at(j)*m_system->getParticles().at(i)->getPosition().at(j);
        }
        r2 += particles.at(i)->getPosition().back() * particles.at(i)->getPosition().back() * beta*beta;
    }
    potentialEnergy = r2*0.5;

    // one body kinetic energy
    double oneBodyLaplacian = r2 * 4.0 * alpha * alpha
                              - 2 * ((m_system->getNumberOfDimensions() - 1) * alpha + alpha * beta)
                              * m_system->getNumberOfParticles();

    // interaction laplacian and cross term
    double interactionLaplacian = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        double r_i2 = 0;
        for (int k=0; k < m_system->getNumberOfDimensions(); k++)
            r_i2 += m_system->getParticles().at(i)->getPosition().at(k)*
                    m_system->getParticles().at(i)->getPosition().at(k);

        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, temp;
        for (int j=0; j < i; j++) {
            double r_ij = m_system->getInterparticleDistance(i,j);
            temp = interactionRange / ((r_ij - interactionRange) * r_ij);
            sum1 += temp;
            sum2 += temp*temp;
            sum3 += (interactionRange -2*r_ij)  * interactionRange /
                    ( (r_ij - interactionRange) * (r_ij - interactionRange) * r_ij * r_ij);
            double r_iDotr_j = 0;
            for (int k=0; k < m_system->getNumberOfDimensions()-1; k++)
                r_iDotr_j += m_system->getParticles().at(i)->getPosition().at(k)*m_system->getParticles().at(j)->getPosition().at(k);
            int k = m_system->getNumberOfDimensions()-1;
            r_iDotr_j += beta*m_system->getParticles().at(i)->getPosition().at(k)*m_system->getParticles().at(j)->getPosition().at(k);

            sum4 -= 4 * alpha * interactionRange * r_i2 * r_iDotr_j /
                    ((r_ij - interactionRange) * r_ij * r_ij);
        }
        for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
            double r_ij = m_system->getInterparticleDistance(j,i);
            temp = interactionRange / ((r_ij - interactionRange) * r_ij);
            sum1 += temp;
            sum3 += (interactionRange -2*r_ij) * interactionRange /
                    ( (r_ij - interactionRange) * (r_ij - interactionRange) * r_ij * r_ij);
            double r_iDotr_j = 0;
            for (int k=0; k < m_system->getNumberOfDimensions() - 1; k++) {
                r_iDotr_j += m_system->getParticles().at(i)->getPosition().at(k)*m_system->getParticles().at(j)->getPosition().at(k);
            }
            r_iDotr_j += beta * m_system->getParticles().at(i)->getPosition().back()*
                                    m_system->getParticles().at(j)->getPosition().back();
            sum4 -= 4 * alpha * interactionRange * r_i2 * r_iDotr_j /
                    ((r_ij - interactionRange) * r_ij * r_ij);
        }

        sum2 = sum1*sum1;
        sum1 *= (m_system->getNumberOfDimensions()-1)/sqrt(r_i2);
        interactionLaplacian += sum1 + sum2 + sum3 + sum4;
    }

    kineticEnergy = -0.5 * (oneBodyLaplacian + interactionLaplacian);

    // RETURN THE TOTAL ENERGY
    return kineticEnergy + potentialEnergy ;
}

// COMPUTE THE LOCAL ENERGY NUMERICALLY
double HarmonicOscillator::computeLocalEnergyNumeric(std::vector<Particle*> particles) {
    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    double interactionPotential = 0.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (unsigned int j=0; j < i; j++) {
            if(m_system->getInterparticleDistance(i,j) < m_system->getWaveFunction()->getParameters().at(2)){
                interactionPotential = 1e100;
                break;
            }
        }
    }
    if(interactionPotential != 0) return interactionPotential;
    double r2 = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
          r2 += particles.at(i)->getPosition().at(j)*
                particles.at(i)->getPosition().at(j);
        }
        r2 += particles.at(i)->getPosition().back() *
              particles.at(i)->getPosition().back() *
              m_system->getWaveFunction()->getParameters().at(1);
    }

    potentialEnergy = r2*0.5;
    kineticEnergy = -0.5 * m_system->getWaveFunction()->computeDoubleDerivative(particles)
                         / m_system->getWaveFunction()->evaluate(particles);
    return kineticEnergy + potentialEnergy;
}

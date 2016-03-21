#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double deriveStep,
                               double alpha, double beta, double interactionRange)
               : WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(interactionRange);
    m_interactionRange = interactionRange;
    m_deriveStep = deriveStep;
}


// COMPUTE THE VALUE OF THE WAVE FUNCTION IN THE CURRENT STATE
double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    // INTERACTION
    double InteractingPart = 1.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (unsigned int j=0; j < i; j++) {
                if(m_system->getInterparticleDistance(i,j) < m_interactionRange){ InteractingPart = 0.0; break;}
                InteractingPart *= 1 - m_interactionRange / m_system->getInterparticleDistance(i,j);
        }
    }

    // ONE BODY PART
    double r2 = 0.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
            r2 += particles.at(i)->getPosition().at(j) *
                  particles.at(i)->getPosition().at(j);
        }
        r2 += particles.at(i)->getPosition().back() *
              particles.at(i)->getPosition().back() *
              m_system->getWaveFunction()->getParameters().at(1);
    }
    double NonInteractingPart = exp(-m_system->getWaveFunction()->getParameters().at(0) * r2);

    // RETURN THE TOTAL WAVE FUNCTION
    return  NonInteractingPart * InteractingPart;
}



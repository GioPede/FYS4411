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
    m_deriveStep = deriveStep;
}


// COMPUTE THE VALUE OF THE WAVE FUNCTION IN THE CURRENT STATE
double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    // INTERACTION
    double InteractingPart = 1.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (unsigned int j=0; j < i; j++) {
                if(m_system->getDistance()[i][j] < m_parameters[2]){ InteractingPart = 0.0; break;}
                InteractingPart *= 1 - m_parameters[2] / m_system->getDistance()[i][j];;
        }
    }

    // ONE BODY PART
    double r2 = 0.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions() - 1; j++) {
            r2 += particles[i]->getPosition()[j] *
                  particles[i]->getPosition()[j];
        }
        int j =  m_system->getNumberOfDimensions() - 1;
        r2 += particles[i]->getPosition()[j] *
              particles[i]->getPosition()[j] *
              m_system->getWaveFunction()->getParameters()[1];
    }
    double NonInteractingPart = exp(-m_system->getWaveFunction()->getParameters()[0] * r2);

    // RETURN THE TOTAL WAVE FUNCTION
    return  NonInteractingPart * InteractingPart;
}

void SimpleGaussian::setQuantumForce(std::vector<class Particle*> particles){
    // CREATE THE MATRIX
    m_quantumForce.resize(m_system->getNumberOfDimensions());
    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
        m_quantumForce[i].resize(2*m_system->getNumberOfParticles());

    // ASSIGN THE VALUES
    double a = m_parameters[2];
    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
        for(int j = 0; j < m_system->getNumberOfParticles(); j++){
            m_quantumForce[i][j] = -2*particles[j]->getPosition()[i]*m_parameters[0];
            for(int k = 0; k < j; k++){
                double rjk= m_system->getDistance()[j][k];
                m_quantumForce[i][j]+= a* (particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
            }
            for(int k = j+1; k < m_system->getNumberOfParticles(); k++){
                double rjk= m_system->getDistance()[k][j];
                m_quantumForce[i][j]+= a*(particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
            }
            m_quantumForce[i][j]*=2;
        }
}

// COMPUTE QUANTUM FORCE FOR THE SUGGESTED MOVE
void SimpleGaussian::computeNewQuantumForce(std::vector<class Particle*> particles){
    double a = m_parameters[2];
    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
        for(int j = 0; j < m_system->getNumberOfParticles(); j++){

            m_quantumForce[i][j+m_system->getNumberOfParticles()] = -2*particles[j]->getPosition()[i]*m_parameters[0];

            for(int k = 0; k < j; k++){
                double rjk= m_system->getDistance()[j][k];
                m_quantumForce[i][j+m_system->getNumberOfParticles()]+= a*(particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
            }
            for(int k = j+1; k < m_system->getNumberOfParticles(); k++){
                double rjk= m_system->getDistance()[k][j];
                m_quantumForce[i][j+m_system->getNumberOfParticles()]+= a*(particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
            }
            m_quantumForce[i][j+m_system->getNumberOfParticles()]*=2;
        }
}

void SimpleGaussian::substituteQuantumForce(){
    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
        for(int j = 0; j < m_system->getNumberOfParticles(); j++){
            m_quantumForce[i][j] = m_quantumForce[i][j+m_system->getNumberOfParticles()];
        }
}

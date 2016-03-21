#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <cstdio>

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}


// COMPUTE THE DERIVATIVE NUMERICALLY
double WaveFunction::computeDerivative(int particleIndex, int Dimension) {
    double psi1 = 0.0;
    m_system->getParticles().at(particleIndex)->adjustPosition(m_deriveStep, Dimension);
    psi1 += m_system->getWaveFunction()->evaluate(m_system->getParticles());

    m_system->getParticles().at(particleIndex)->adjustPosition(-2.0*m_deriveStep, Dimension);
    psi1 -= m_system->getWaveFunction()->evaluate(m_system->getParticles());

    m_system->getParticles().at(particleIndex)->adjustPosition(m_deriveStep, Dimension);
    return psi1 / m_deriveStep / 2.0;
}

// COMPUTE THE SECOND DERIVATIVE NUMERICALLY
double WaveFunction::computeDoubleDerivative(std::vector<class Particle*> particles) {
    double psi2 = 0.0;
    double sum = 0.0;
    double diff = 0.0;
    for (unsigned int i=0; i < particles.size(); i++) {
        for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
            diff -= 2.0 * m_system->getWaveFunction()->evaluate(particles);

            particles.at(i)->adjustPosition(m_deriveStep, j);
            sum += m_system->getWaveFunction()->evaluate(particles);

            particles.at(i)->adjustPosition(-2.0*m_deriveStep, j);
            sum += m_system->getWaveFunction()->evaluate(particles);

            particles.at(i)->adjustPosition(m_deriveStep, j);

        }
    }
    psi2 = (sum + diff) / (m_deriveStep * m_deriveStep);
    return psi2;
}

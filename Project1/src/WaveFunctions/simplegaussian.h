#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double deriveStep, double alpha, double beta, double interactionRange);
    double evaluate(std::vector<class Particle*> particles);

private:
    double m_interactionRange = 0;
};

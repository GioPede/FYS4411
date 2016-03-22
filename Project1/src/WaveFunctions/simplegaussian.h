#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double deriveStep, double alpha, double beta, double interactionRange);
    double evaluate(std::vector<class Particle*> particles);
    void setQuantumForce(std::vector<class Particle*> particles);
    void computeNewQuantumForce(std::vector<class Particle*> particles);
    void substituteQuantumForce();
    void clearQuantumForce();
};

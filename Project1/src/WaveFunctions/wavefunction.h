#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double>& getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    double computeDerivative(int particleIndex, int Dimension);
    double computeDoubleDerivative(std::vector<class Particle*> particles);


protected:
    int     m_numberOfParameters = 0;
    double  m_deriveStep = 0.0001;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};


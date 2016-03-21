#pragma once
#include <stdio.h>
#include <stdlib.h>

class Sampler {
public:
    Sampler(class System* system, bool print);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }

    double getCumulativeEnergy() const;
    double getCumulativeDePsideAlpha() const;
    double getCumulativeDePsideAlphaELocal() const;

    void initializeFirstStep();
private:
    bool m_print;
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_variance = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeDePsideAlpha= 0;
    double  m_cumulativeDePsideAlphaELocal = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_localEnergy = 0;
    double  m_dePsideAlpha = 0;
    double  m_dePsideAlphaELocal = 0;
    double  m_numberAccepted = 0;
    class System* m_system = nullptr;
    FILE *  m_output = nullptr;
};

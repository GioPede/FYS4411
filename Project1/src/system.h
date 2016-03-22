#pragma once
#include <vector>

class System {
public:
    bool metropolisStep             ();
    bool metropolisStepImportanceSamplingNumeric();
    void runMetropolisSteps         (int numberOfMetropolisSteps, bool print);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>&   getParticles()      { return m_particles; }
    std::vector<std::vector<double>>&   getDistance()      { return m_distanceMatrix; }

    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }


    void setDistanceMatrix();
    void updateDistanceMatrix();
    void updateDistanceMatrix(int particleIndex);
    bool metropolisStepImportanceSampling();
    void setDeltaR();
    void computeDeltaR();

    void setSampler(class Sampler *sampler);

private:
    std::vector<std::vector<double>> m_distanceMatrix;
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    std::vector<std::vector<double>>m_deltaR;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_waveOld= 0.0;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    class Random*                   m_random = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};


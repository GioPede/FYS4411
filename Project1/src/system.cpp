#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>


// BASIC METROPOLIS STEP, NO IMPORTANCE SAMPLING
bool System::metropolisStep() {
    // SUGGEST A MOVE
    double deltaR = (m_random->nextDouble() - 0.5) * m_stepLength;
    double dimensionOfChange = m_random->nextInt(m_numberOfDimensions);
    int selectedParticleIndex = m_random->nextInt(m_numberOfParticles);
    Particle* selectedParticle = m_particles[selectedParticleIndex];
    selectedParticle->adjustPosition(deltaR, dimensionOfChange);
    updateDistanceMatrix(selectedParticleIndex);

    // COMPUTE NEW WAVE FUNCTION
    double newWaveValue = m_waveFunction->evaluate(m_particles);

    // ACCEPT OR REJECT THE MOVE
    if (newWaveValue*newWaveValue/(m_waveOld*m_waveOld) > m_random->nextDouble()){
        m_waveOld = newWaveValue;
        return true;
    }
    else {
        selectedParticle->adjustPosition(-deltaR, dimensionOfChange);
        updateDistanceMatrix(selectedParticleIndex);
        return false;
    }
}

// ADVANCED METROPOLIS STEP WITH IMPORTANCE SAMPLING
bool System::metropolisStepImportanceSampling() {

    // moving the particles
    computeDeltaR();
    for(int i = 0; i < m_numberOfDimensions; i++){
        for (int j=0;j<m_numberOfParticles;j++){
            m_particles[j]->adjustPosition(m_deltaR[i][j], i);
        }
    }

    // COMPUTE NEW VALUES
    updateDistanceMatrix();
    m_waveFunction->computeNewQuantumForce(m_particles);
    double newWaveValue = m_waveFunction->evaluate(m_particles);

    // compute the green function
    double greenFunction = 0.0;
    for(int i = 0; i < m_numberOfDimensions; i++)
        for(int j= 0;j< m_numberOfParticles;j++){
        greenFunction += 0.5 * (m_waveFunction->getQuantumForce()[i][j+m_numberOfParticles] + m_waveFunction->getQuantumForce()[i][j])
                * (0.5 * 0.5 * m_stepLength *(m_waveFunction->getQuantumForce()[i][j] -m_waveFunction->getQuantumForce()[i][j+m_numberOfParticles])
                   -m_deltaR[i][j]);
    }
    greenFunction = exp(greenFunction);

    // accept or reject the move
    if (greenFunction * newWaveValue*newWaveValue/(m_waveOld*m_waveOld) > m_random->nextDouble()){
        m_waveOld = newWaveValue;
        m_waveFunction->substituteQuantumForce();
        return true;
    }
    else {
        for(int i = 0; i < m_numberOfDimensions; i++){
            for (int j=0;j<m_numberOfParticles;j++){
                m_particles[j]->adjustPosition(-m_deltaR[i][j], i);
            }
        }
        updateDistanceMatrix();
        return false;
    }
}

// MAIN FUNCTION OF THE PROGRAM, LOOP OVER THE MC CYCLES
void System::runMetropolisSteps(int numberOfMetropolisSteps, bool print) {
    // INITIALIZE THE WHOLE SYSTEM
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this, print);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_random = new Random();
    m_random->setSeed(time(NULL));

std::cout<< "hi!" << std::endl;
    setDistanceMatrix();
    std::cout<< "hi!" << std::endl;
    setDeltaR();
    m_waveFunction->setQuantumForce(m_particles);
    m_sampler->initializeFirstStep();
    m_waveOld = m_waveFunction->evaluate(m_particles);


    // START THE METROPOLIS ALGORITHM
    // thermalization
    for (int i=0; i < numberOfMetropolisSteps * m_equilibrationFraction; i++)
        metropolisStepImportanceSampling();
    // sampling
    for (int i = numberOfMetropolisSteps * m_equilibrationFraction; i < numberOfMetropolisSteps; i++){
        bool acceptedStep = metropolisStepImportanceSampling();
        m_sampler->sample(acceptedStep);
    }

    // COMPUTE AVERAGES AND EXIT
    m_sampler->computeAverages();
}

// COMUTE A RANDOM MOVE BASED ON THE CURRENT QUANTUM FORCE
void System::computeDeltaR(){
    for(int i = 0; i < m_numberOfDimensions; i++)
        for(int j = 0; j < m_numberOfParticles; j++){
            m_deltaR[i][j] = m_random->nextGaussian(0, 1) * sqrt(m_stepLength) + 0.5 * m_stepLength * m_waveFunction->getQuantumForce()[i][j];
        }
}

// UPDATE THE CURRENT DISTANCE MATRIX
void System::updateDistanceMatrix() {
    for(int i=0;i< m_numberOfParticles;i++){
        for(int j = 0; j < i; j++){
            m_distanceMatrix[i][j]=0.0;
            for(int k = 0; k < m_numberOfDimensions; k++){
                m_distanceMatrix[i][j] += (m_particles[i]->getPosition()[k] - m_particles[j]->getPosition()[k]) *
                                          (m_particles[i]->getPosition()[k] - m_particles[j]->getPosition()[k]);}
             m_distanceMatrix[i][j] = sqrt(m_distanceMatrix[i][j]);
        }
    }

}

// UPDATE THE CURRENT DISTANCE MATRIX FOR JUST ONE PARTICLE MOVEMENT
void System::updateDistanceMatrix(int particleIndex) {
    for(int j = 0; j < particleIndex; j++){
        m_distanceMatrix[particleIndex][j] = 0.0;
        for(int k = 0; k < m_numberOfDimensions; k++)
            m_distanceMatrix[particleIndex][j] += (m_particles[particleIndex]->getPosition()[k] - m_particles[j]->getPosition()[k]) *
                                                  (m_particles[particleIndex]->getPosition()[k] - m_particles[j]->getPosition()[k]);
        m_distanceMatrix[particleIndex][j] = sqrt(m_distanceMatrix[particleIndex][j]);
    }
    for(int i = particleIndex; i < m_numberOfParticles; i++){
        m_distanceMatrix[i][particleIndex] = 0.0;
        for(int k = 0; k < m_numberOfDimensions; k++)
            m_distanceMatrix[i][particleIndex] += (m_particles[i]->getPosition()[k] - m_particles[particleIndex]->getPosition()[k]) *
                                                  (m_particles[i]->getPosition()[k] - m_particles[particleIndex]->getPosition()[k]);
        m_distanceMatrix[i][particleIndex] = sqrt(m_distanceMatrix[i][particleIndex]);
    }
}

void System::setDeltaR(){
    // CREATE THE MATRIX
    m_deltaR.resize(m_numberOfDimensions);
    for(int i = 0; i < m_numberOfDimensions; i++)
        m_deltaR[i].resize(m_numberOfParticles);
}



void System::setDistanceMatrix() {
    // CREATE THE MATRIX
    m_distanceMatrix.resize(m_numberOfParticles);
    for(int i = 0; i < m_numberOfParticles; i++)
        m_distanceMatrix[i].resize(m_numberOfParticles);

    for(int i = 0; i < m_numberOfParticles; i++)
        for(int j = 0; j < i; j++){
            m_distanceMatrix[i][j] = 0.0;
            for(int k = 0; k < m_numberOfDimensions; k++){
                m_distanceMatrix[i][j] += (m_particles[i]->getPosition()[k] - m_particles[j]->getPosition()[k]) *
                                          (m_particles[i]->getPosition()[k] - m_particles[j]->getPosition()[k]);
            }
            m_distanceMatrix[i][j] = sqrt(m_distanceMatrix[i][j]);
        }
}

void System::setSampler(Sampler *sampler){
    m_sampler = sampler;
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}





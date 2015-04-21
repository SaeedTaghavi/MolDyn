/*
 * MixedLJPotential.hpp
 *
 * A multi-component Lennard-Jones Potential
 *
 *  Created on: Apr 17, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef MIXEDLJPOTENTIAL_HPP_
#define MIXEDLJPOTENTIAL_HPP_

#include "Potential.hpp"
#include <vector>
#include "SymmetricMatrix.hpp"
#include "Component.hpp"
#include <cmath>
#include <stdexcept>

namespace MolDyn {

template<class T>
class MixedLJPotential: public Potential<T> {
private:
	unsigned int maxComponents;
	unsigned int componentCounter;
	std::vector<T> epsilon;
	std::vector<T> sigma;
	std::vector<unsigned int> n_particles;
	std::vector<T> densities;
	/*
	 * A symmetric matrix with the mixed sigma parameters
	 */
	SymmetricMatrix<T> mixedSigma;

	/*
	 * A Symmetric matrix with the mixed epsilon parameters
	 */
	SymmetricMatrix<T> mixedEpsilon;
public:
	/*
	 * Class constructor
	 *
	 * @param[in]	_maxComponents: The number of components in the system
	 * @param[in]	_cutOff: The LJ cutoff distance
	 */
	MixedLJPotential(const unsigned int& _maxComponents, const T& _cutOff);
	/*
	 * evaluate(const T& _distance, T& _ePot, T& _force,
	 * const unsigned int _species1, const unsigned int _species2)
	 *
	 * Evaluates the energy and magnitude of the force according to the
	 * Lennard-Jones pairwise potential for two (not necessarily equal)
	 * particles.
	 *
	 * @param[in]		_distance: The distance between particles
	 * @param[in,out]	_ePot: The potential energy between particles
	 * @param[in,out]	_force: The magnitude of the force
	 * @param[in]		_species1: The species of the first particle (if different from)
	 * @param[in]		_species2: The species of the second particle
	 */
	void evaluate(const T& _distance, T& _ePot, T& _force,
			const unsigned int _species1, const unsigned int _species2);
	/*
	 * addComponent(const T& _sigma, const T& _epsilon, const int& _nParticles,
	 * 		const T& _density)
	 * Adds a component to the class.
	 * The limit of components added is set in the constructor class. Upon reaching
	 * the _maxComponents, next method calls will be ignored
	 *
	 * @param[in]	_sigma: The hard shell diameter
	 * @param[in]	_epsilon: The epsilon LJ parameter
	 * @param[in]	_nParticles: The number of particles holding these parameters
	 * 				(used when computing Long Range Corrections)
	 * @param[in]	_density: The density of the particles in the component
	 */
	void addComponent(const T& _sigma, const T& _epsilon,
			const int& _nParticles, const T& _density);
	/*
	 * addComponent(const MolDyn::Component<T>& component)
	 * Adds a component using the Component envelope
	 * @see addComponent(const T& _sigma, const T& _epsilon, const int& _nParticles,
	 * const T& _density)
	 */
	void addComponent(const MolDyn::Component<T>& component);
	/*
	 * @param[in]	type: The index of the component
	 * @return 		sigma of the selected component
	 */
	T getSigma(const unsigned int& type = 0);
	/*
	 * @param[in]	type: The index of the component
	 * @return 		Epsilon of the selected component
	 */
	T getEpsilon(const unsigned int& type = 0);
	/*
	 * @param[in]	type: The index of the component
	 * @return 		The cutoff distance of the selected component
	 */
	T getCutOff(const unsigned char type1, const unsigned char type2);
	/*
	 * validateCutoff(const T& length)
	 *
	 * Checks if the cutoff length exceeds the size of (half) the length
	 * of the simulation box
	 *
	 * @param[in]	length: half the length of the simulation box
	 * @return		true if the length is greater than LJ cutoff
	 * 				and false otherwise
	 */
	bool validateCutoff(const T& length);

	T getELRC() {
		return Potential<T>::ELRC;
	}
	;
	T getPLRC() {
		return Potential<T>::PLRC;
	}
	;
	virtual ~MixedLJPotential();
private:
	/*
	 * Calculate the Long range corrections for Energy and Pressure
	 * and store them in ELRC and PLRC respectively
	 */
	void setLongRangeCorrections();
};

template<class T>
inline MolDyn::MixedLJPotential<T>::MixedLJPotential(
		const unsigned int& _maxComponents, const T& cutoff) :
		Potential<T>(cutoff), maxComponents(_maxComponents), epsilon(
				_maxComponents), sigma(_maxComponents), mixedEpsilon(
				_maxComponents), mixedSigma(_maxComponents), n_particles(
				_maxComponents), densities(_maxComponents), componentCounter(0) {
}

template<class T>
inline void MolDyn::MixedLJPotential<T>::addComponent(const T& _sigma,
		const T& _epsilon, const int& _nParticles, const T& _density) {
	if (componentCounter < maxComponents) {
		sigma[componentCounter] = _sigma;
		epsilon[componentCounter] = _epsilon;
		n_particles[componentCounter] = _nParticles;
		densities[componentCounter] = _density;
		++componentCounter;
	}

	// Update the mixed parameters matrices
	if (componentCounter == maxComponents) {
		for (unsigned int i = 0; i < maxComponents; ++i) {
			for (unsigned int j = i; j < maxComponents; ++j) {
				mixedEpsilon(i, j) = mixEpsilon(epsilon[i], epsilon[j]);
				mixedSigma(i, j) = mixSigma(sigma[i], sigma[j]);
			}
		}

		std::cout << "Mixing Sigma:";
		std::cout << mixedSigma;

		std::cout << "Mixing Epsilon:";
		std::cout << mixedEpsilon;

		setLongRangeCorrections();
	}
}

template<class T>
inline void MolDyn::MixedLJPotential<T>::evaluate(const T& distance, T& epot,
		T& force, const unsigned int s1, const unsigned int s2) {
	T sr6 = pow(mixedSigma(s1, s2) / distance, 6.0);
	T sr12 = sr6 * sr6;
	force = (24.0 * mixedEpsilon(s1, s2) / distance) * (2 * sr12 - sr6);
	epot = 4.0 * mixedEpsilon(s1, s2) * (sr12 - sr6);
}

template<class T>
inline T MolDyn::MixedLJPotential<T>::getSigma(const unsigned int& type) {
	if (type > maxComponents) {
		throw std::out_of_range("Sigma index is out of bounds");
	}
	return sigma[type];
}

template<class T>
inline T MolDyn::MixedLJPotential<T>::getEpsilon(const unsigned int& type) {
	if (type > maxComponents) {
		throw std::out_of_range("Epsilon index is out of bounds");
	}
	return epsilon[type];
}

template<class T>
inline T MolDyn::MixedLJPotential<T>::getCutOff(const unsigned char type1 = 0,
		const unsigned char type2 = 0) {
	if (type1 >= maxComponents || type2 >= maxComponents) {
		throw std::out_of_range("Cutoff index is out of bounds");
	}
	return mixSigma(sigma[type1], sigma[type2]) * Potential<T>::cutOff;
}

template<class T>
inline MolDyn::MixedLJPotential<T>::~MixedLJPotential() {
	std::vector<T>().swap(sigma);
	std::vector<T>().swap(epsilon);
	std::vector<unsigned int>().swap(n_particles);
	std::vector<T>().swap(densities);

}

template<class T>
inline bool MixedLJPotential<T>::validateCutoff(const T& length) {
	for (unsigned int i = 0; i < maxComponents; ++i) {
		if (Potential<T>::cutOff * sigma[i] >= length) {
			return false;
		}
	}
	return true;
}

template<class T>
inline void MixedLJPotential<T>::addComponent(
		const MolDyn::Component<T>& component) {
	addComponent(component.sigma, component.epsilon, component.nParticles,
			component.density);
}

template<class T>
inline void MolDyn::MixedLJPotential<T>::setLongRangeCorrections() {
	/* CALCULATE LONG RANGE CORRECTIONS FOR ENERGY AND PRESSURE APRIORI
	 * FOR EACH COMPONENT, THEN SUM
	 */
	Potential<T>::ELRC = 0;
	Potential<T>::PLRC = 0;
	T mSigma, mEpsilon, totDens, totParticles;
	T sc3, sc9, s3, ec, pc;

	totDens = 0;
	totParticles = 0;
	mSigma = 0;
	mEpsilon = 1;

	// Average sigma, epsilon and sum densities
	for (unsigned int i = 0; i < maxComponents; ++i) {
		mSigma += sigma[i];
		mEpsilon *= epsilon[i];
		totDens += densities[i];
		totParticles += n_particles[i];
	}

	mSigma /= maxComponents;
	mEpsilon = pow(mEpsilon, 1.0 / maxComponents);

	sc3 = pow(mSigma / (Potential<T>::cutOff * mSigma), 3.0);
	sc9 = pow(sc3, 3.0);
	s3 = pow(mSigma, 3.0);
	ec = (1.0 / 3.0) * sc9 - sc3;
	ec *= (8.0 / 3.0) * PI * mEpsilon * s3 * totParticles * totDens;
	pc = (2.0 / 3.0) * sc9 - sc3;
	pc *= (16.0 / 3.0) * PI * mEpsilon * s3 * totDens * totDens;
	Potential<T>::ELRC += ec;
	Potential<T>::PLRC += pc;
}

} /* namespace MolDyn */

#endif /* MIXEDLJPOTENTIAL_HPP_ */

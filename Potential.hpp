/*
 * Potential.h
 *
 *  Created on: Mar 30, 2015
 *      Author: erick
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_
#include "constants.hpp"

namespace MolDyn {

template<class T>
class Potential {
public:
	T cutoffBlowsAt;
protected:
	T cutOff;
	T ELRC, PLRC;
public:
	Potential(const T& cutoff);
	virtual T getSigma(const unsigned int& type = 0) = 0;
	virtual T getEpsilon(const unsigned int& type = 0) = 0;
	virtual T getCutOff(const unsigned char type1 = 0,
			const unsigned char type2 = 0) = 0;
	virtual bool validateCutoff(const T& length) = 0;
	virtual void evaluate(const T& distance, T& epot, T& force,
			const unsigned int s1, const unsigned int s2) = 0; // Interface method
	virtual T getELRC() = 0;
	virtual T getPLRC() = 0;
	virtual ~Potential();
};

template<class T>
inline MolDyn::Potential<T>::Potential(const T& cutoff) :
		cutOff(cutoff) {
}

template<class T>
inline MolDyn::Potential<T>::~Potential() {
}

/*
 * Pairwise mix sigma according to the Lorentz-Berthelot
 * mixing rules:
 *
 * sigma_{12} = (sigma_1 + sigma_2) / 2
 *
 * @param[in]	sigma1
 * @param[in]	sigma2
 * @return		The mixed sigma
 */
template<typename T>
T mixSigma(const T& sigma1, const T& sigma2) {
	return 0.5 * (sigma1 + sigma2);
}

/*
 * Pairwise mix epsilon according to the Lorentz-Berthelot
 * mixing rules:
 *
 * epsilon_{12} = sqrt(epsilon_1 * epsilon_2)
 *
 * @param[in]	epsilon1
 * @param[in]	epsilon2
 * @return		The mixed epsilon
 */
template<class T>
T mixEpsilon(const T& eps1, const T& eps2) {
	return sqrt(eps1 * eps2);
}

} /* namespace MolDyn */

#endif /* POTENTIAL_HPP_ */

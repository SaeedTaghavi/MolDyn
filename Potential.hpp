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
protected:
	T epsilon;
	T sigma;
	T cutOff;
	T ELRC, PLRC;
public:
	virtual void evaulate(const T& pos, T& E, T& F) = 0; // Interface method
	virtual T getELRC(const T& natoms, const T& density) = 0; // Interface method
	virtual T getPLRC(const T& density) = 0; // Interface method
	virtual ~Potential();
	T getSigma() const {
		return sigma;
	}

	T getEpsilon() const {
		return epsilon;
	}
	T getCutOff() const {
		return cutOff;
	}
	bool validateCutoff(const T& length) const {
		return cutOff <= length;
	}
};

template<class T>
inline MolDyn::Potential<T>::~Potential() {
}

} /* namespace MolDyn */

#endif /* POTENTIAL_HPP_ */

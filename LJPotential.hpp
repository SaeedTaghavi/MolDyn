/*
 * LJPotential.hpp
 *
 *  Created on: Mar 30, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef LJPOTENTIAL_HPP_
#define LJPOTENTIAL_HPP_

#include "Potential.hpp"
#include <cmath>

namespace MolDyn {

template<class T>
class LJPotential: public Potential<T> {
private:
public:
	LJPotential(const T& s, const T& e);
	void evaulate(const T& pos, T& E, T& F);
	T getELRC(const T& natoms, const T& density);
	T getPLRC(const T& density);
	T energyMin() const;
	T rmin() const;
	void setCutOff(T cut);
	virtual ~LJPotential();
private:
	void setLongRangeCorrections();
};

template<class T>
inline MolDyn::LJPotential<T>::LJPotential(const T& s, const T& e) {
	Potential<T>::sigma = s;
	Potential<T>::epsilon = e;
	Potential<T>::cutOff = 2.5 * Potential<T>::sigma;
	setLongRangeCorrections();
}

template<class T>
inline void MolDyn::LJPotential<T>::evaulate(const T& pos, T& E, T& F) {
	T sr6 = pow(Potential<T>::sigma / pos, 6.0);
	T sr12 = sr6 * sr6;
	F = (24.0 * Potential<T>::epsilon / pos) * (2 * sr12 - sr6);
	E = 4.0 * Potential<T>::epsilon * (sr12 - sr6);
}

template<class T>
inline T MolDyn::LJPotential<T>::energyMin() const {
	return -Potential<T>::epsilon;
}

template<class T>
inline T MolDyn::LJPotential<T>::rmin() const {
	return Potential<T>::sigma * pow(2.0, 1.0 / 6.0);
}

template<class T>
inline T MolDyn::LJPotential<T>::getELRC(const T& natoms, const T& density) {
	return Potential<T>::ELRC * natoms * density;
}

template<class T>
inline T MolDyn::LJPotential<T>::getPLRC(const T& density) {
	return Potential<T>::PLRC * density * density;
}

template<class T>
inline void MolDyn::LJPotential<T>::setCutOff(T cut) {
	Potential<T>::cutOff = Potential<T>::sigma * cut;
	setLongRangeCorrections();
}

template<class T>
inline LJPotential<T>::~LJPotential() {

}

template<class T>
inline void MolDyn::LJPotential<T>::setLongRangeCorrections() {
	// CALCULATE LONG RANGE CORRECTIONS FOR ENERGY AND PRESSURE APRIORI
	T sl3 = pow(Potential<T>::sigma / Potential<T>::cutOff, 3.0);
	T sl9 = pow(sl3, 3.0);
	T s3 = pow(Potential<T>::sigma, 3.0);
	Potential<T>::ELRC = (1.0 / 3.0) * sl9 - sl3;
	Potential<T>::ELRC *= (8.0 / 3.0) * PI * Potential<T>::epsilon * s3;
	Potential<T>::PLRC = (2.0 / 3.0) * sl9 - sl3;
	Potential<T>::PLRC *= (16.0 / 3.0) * PI * Potential<T>::epsilon * s3;
}

} /* namespace MolDyn */

#endif /* LJPOTENTIAL_HPP_ */

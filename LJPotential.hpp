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
	unsigned int n_atoms;
	T density;
	T epsilon;
	T sigma;
public:
	LJPotential(const T& s, const T& e, const unsigned int& n, const T& dens,
			const T& cutoff);
	void evaulate(const T& distance, T& epot, T& force,
			const unsigned int type1, const unsigned int type2);
//	T getELRC(const T& natoms, const T& density);
//	T getPLRC(const T& density);
	T energyMin() const;
	T rmin() const;
	void setCutOff(T cut);
	T getSigma(const unsigned int& type) const;
	T getEpsilon(const unsigned int& type = 0) const;
	bool validateCutoff(const T& length) const;
	virtual ~LJPotential();
private:
	void setLongRangeCorrections();
};

template<class T>
inline MolDyn::LJPotential<T>::LJPotential(const T& s, const T& e,
		const unsigned int& n, const T& dens, const T& cutoff = 2.5) :
		n_atoms(n), density(dens), sigma(s), epsilon(e), Potential<T>(
				cutoff * sigma) {
	setLongRangeCorrections();
}

template<class T>
inline void MolDyn::LJPotential<T>::evaulate(const T& distance, T& epot,
		T& force, const unsigned int type1 = 0,
		const unsigned int type2 = 0) {
	T sr6 = pow(sigma / distance, 6.0);
	T sr12 = sr6 * sr6;
	force = (24.0 * epsilon / distance) * (2 * sr12 - sr6);
	epot = 4.0 * epsilon * (sr12 - sr6);
}

template<class T>
inline T MolDyn::LJPotential<T>::energyMin() const {
	return -epsilon;
}

template<class T>
inline T MolDyn::LJPotential<T>::rmin() const {
	return sigma * pow(2.0, 1.0 / 6.0);
}

//template<class T>
//inline T MolDyn::LJPotential<T>::getELRC(const T& natoms, const T& density) {
//	return Potential<T>::ELRC * natoms * density;
//}
//
//template<class T>
//inline T MolDyn::LJPotential<T>::getPLRC(const T& density) {
//	return Potential<T>::PLRC * density * density;
//}

template<class T>
inline void MolDyn::LJPotential<T>::setCutOff(T cut) {
	Potential<T>::cutOff = sigma * cut;
	setLongRangeCorrections();
}

template<class T>
inline LJPotential<T>::~LJPotential() {

}

template<class T>
inline T LJPotential<T>::getSigma(const unsigned int& type = 0) const {
	return sigma;
}

template<class T>
inline T LJPotential<T>::getEpsilon(const unsigned int& type) const {
	return epsilon;
}

template<class T>
inline bool LJPotential<T>::validateCutoff(const T& length) const {
	return Potential<T>::cutOff <= length;
}

template<class T>
inline void MolDyn::LJPotential<T>::setLongRangeCorrections() {
	// CALCULATE LONG RANGE CORRECTIONS FOR ENERGY AND PRESSURE APRIORI
	T sl3 = pow(sigma / Potential<T>::cutOff, 3.0);
	T sl9 = pow(sl3, 3.0);
	T s3 = pow(sigma, 3.0);
	Potential<T>::ELRC = (1.0 / 3.0) * sl9 - sl3;
	Potential<T>::ELRC *= (8.0 / 3.0) * PI * epsilon * s3 * n_atoms * density;
	Potential<T>::PLRC = (2.0 / 3.0) * sl9 - sl3;
	Potential<T>::PLRC *= (16.0 / 3.0) * PI * epsilon * s3 * density * density;
}

} /* namespace MolDyn */

#endif /* LJPOTENTIAL_HPP_ */

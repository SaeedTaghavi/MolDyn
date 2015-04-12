/*
 * Integrator.hpp
 *
 *  Created on: Mar 30, 2015
 *      Author: erick
 */

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_
#include "Potential.hpp"
#include <vector>
#include <cmath>
#include "Vector3d.hpp"
#include "constants.hpp"

namespace MolDyn {

template<class T>
class Integrator {
public:
	Potential<T>& potential;
protected:
	T dt;
	T mass;
	T density;
	T temp;
public:
	Integrator(Potential<T>& p, const T& m, const T& dens, const T& tstep);
	virtual void advance(std::vector<Vector3d<T> >& r,
			std::vector<Vector3d<T> >& v, std::vector<Vector3d<T> >& f,
			T& K, T& U, T& E, T& P, const T& len, bool rescale) = 0;
	void forces(std::vector<Vector3d<T> >& r,
			std::vector<Vector3d<T> >& f, T& EPOT, T& PRESS, const T& len);
	void currentEnergy(std::vector<Vector3d<T> >& v, T& E, T& U, T& K,
			T& P);
	void GaussConstraint(std::vector<Vector3d<T> >& v,
			std::vector<Vector3d<T> >& f);
	void rescale(std::vector<Vector3d<T> >& v);
	T getMass() const {
		return mass;
	}
	void setTemp(const T& t) {
		temp = t;
	}
	T getTemp() const {
		return temp;
	}
	virtual ~Integrator();
};

template<class T>
inline MolDyn::Integrator<T>::Integrator(Potential<T>& p, const T& m,
		const T& dens, const T& tstep) :
		potential(p), mass(m), density(dens), dt(tstep) {
}

template<class T>
inline void MolDyn::Integrator<T>::forces(std::vector<Vector3d<T> >& r,
		std::vector<Vector3d<T> >& f, T& EPOT, T& PRESS, const T& length) {
	unsigned int n_atoms = r.size();
	Vector3d<T> fi;
	Vector3d<T> rij;
	Vector3d<T> fij;
	T FIJ;
	T rijlength;
	T U;

	// SET FORCES ON ALL ATOMS, POTENTIAL ENERGY, AND PRESSURE TO ZERO
	EPOT = 0;
	PRESS = 0;
	for (unsigned int i = 0; i < n_atoms; ++i) {
		f[i].zeros();
	}

	/**
	 * LOOP OVER ALL PAIRWISE INTERACTIONS TO OBTAIN POTENTIAL ENERGY,
	 * CONFIGURATIONAL PART OF PRESSURE, AND FORCES
	 */
	for (unsigned int i = 0; i < n_atoms - 1; ++i) {
		fi = f[i];
		for (unsigned int j = i + 1; j < n_atoms; ++j) {
			rij = r[i] - r[j];

			// APPLY MINIMUM IMAGE CONVENTION
			rij.modulo(length);
			rijlength = rij.length();
			/**
			 * IF DISTANCE IS LESS THAN LJ CUTOFF, CALCULATE FORCE ACCORDING TO
			 * GRADIENT OF POTENTIAL ENERGY AND THE CONFIGURATIONAL COMPONENT
			 * OF PRESSURE BY SUMMING OVER rij*fij
			 */
			if (rijlength < potential.getCutOff()) {
				potential.evaulate(rijlength, U, FIJ);
				PRESS += FIJ * rijlength;
				EPOT += U;
				fij = project(FIJ, rij);
				fi += fij;
				f[j] -= fij;
			} // end if
		} // end for
		f.at(i) = fi;
	} // end for
	PRESS *= density / (3.0 * n_atoms);
}

template<class T>
inline void MolDyn::Integrator<T>::currentEnergy(
		std::vector<Vector3d<T> >& v, T& E, T& U, T& K, T& P) {
	K = 0.0;
	for (unsigned int i = 0; i < v.size(); ++i) {
		K += v[i] * v[i];
	}
	K *= 0.5 * mass;
	E = K + U;
	P += 2.0 * K * density / (3.0 * v.size());
}

template<class T>
inline void MolDyn::Integrator<T>::GaussConstraint(
		std::vector<Vector3d<T> >& v, std::vector<Vector3d<T> >& f) {
	T denSum = 0; // Sum( pi dot F_i / m_i)
	T numSum = 0; // Sum( pi dot P_i / m_i)
	T lambda;
	// Assuming m_i = mass (same for all particles)
	for (unsigned int i = 0; i < v.size(); ++i) {
		denSum += v[i] * v[i];
		numSum += v[i] * f[i];
	}
	lambda = (denSum == 0) ? 1 : numSum / denSum;
//	std::cout << "lambda = " << lambda << std::endl;

	for (unsigned int i = 0; i < v.size(); ++i) {
		v[i] *= lambda;
//		f[i] -= v[i] * lambda * mass;
	}

}

template<class T>
inline void Integrator<T>::rescale(std::vector<Vector3d<T> >& v) {
	T velsq, factor;
	velsq = 0;
	for (unsigned int i = 0; i < v.size(); ++i) {
		velsq += Integrator<T>::mass * (v[i] * v[i]);
	}
	factor = sqrt(3.0 * v.size() * KB * temp / velsq);
	for (unsigned int i = 0; i < v.size(); ++i) {
		v[i] *= factor;
	}
}

template<class T>
inline MolDyn::Integrator<T>::~Integrator() {
}

} /* namespace MolDyn */

#endif /* INTEGRATOR_HPP_ */

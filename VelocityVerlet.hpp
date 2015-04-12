/*
 * VelocityVerlet.hpp
 *
 *  Created on: Mar 31, 2015
 *      Author: erick
 */

#ifndef VELOCITYVERLET_HPP_
#define VELOCITYVERLET_HPP_

#include "Integrator.hpp"

namespace MolDyn {

template<class T>
class VelocityVerlet: public Integrator<T> {
private:
	T dt2;
	T dt2m;
public:
	VelocityVerlet(Potential<T>& p, const T& m, const T& dens, const T& t);
	void advance(std::vector<Vector3d<T> >& r,
			std::vector<Vector3d<T> >& v, std::vector<Vector3d<T> >& f,
			T& K, T& U, T& E, T& P, const T& lbox, bool rescale);
	virtual ~VelocityVerlet();
};

template<class T>
inline MolDyn::VelocityVerlet<T>::VelocityVerlet(Potential<T>& p, const T& m,
		const T& dens, const T& tstep) :
		Integrator<T>(p, m, dens, tstep), dt2(0.5 * tstep), dt2m(
				0.5 * tstep / m) {
}

template<class T>
inline void MolDyn::VelocityVerlet<T>::advance(std::vector<Vector3d<T> >& r,
		std::vector<Vector3d<T> >& v, std::vector<Vector3d<T> >& f, T& K,
		T& U, T& E, T& P, const T& lbox, bool rescale = false) {
	unsigned int n_atoms = r.size();
	/*
	 * ADVANCE POSITIONS FROM T TO T+DT AND VELOCITIES FROM T TO T+DT/2
	 * ACCORDING TO THE VELOCITY VERLET INTEGRATOR:
	 * 	V(T+DT/2) = V(T) + 1/2*DT*A(T)
	 * 	R(T+DT)   = R(T) + DT*V(T) + 1/2*(DT**2)*A(T)
	 * 	V(T+DT)   = V(T+DT/2) + 1/2*DT*A(T+DT)
	 */
	// ADVANCE POSITIONS BY DT AND VELOCITIES BY DT/2
	for (unsigned int i = 0; i < n_atoms; ++i) {
		v[i] += f[i] * dt2m;
		r[i] += v[i] * Integrator<T>::dt;
	}

	// EVALUATE THE FORCES AT STEP T+DT
	Integrator<T>::forces(r, f, U, P, lbox);

	// ADVANCE VELOCITIES BY THE REMAINING DT/2, WITH NEW FORCES
	for (unsigned int i = 0; i < n_atoms; ++i) {
		v[i] += f[i] * dt2m;
	}

	// APPLY PERIODIC BOUNDARY CONDITIONS
	for (unsigned int i = 0; i < n_atoms; ++i) {
		r[i].modulo(lbox);
	}

	if (rescale) {
		Integrator<T>::rescale(v);
	}

	Integrator<T>::currentEnergy(v, E, U, K, P);
}

template<class T>
inline MolDyn::VelocityVerlet<T>::~VelocityVerlet() {
}

} /* namespace MolDyn */

#endif /* VELOCITYVERLET_HPP_ */

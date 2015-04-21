/*
 * VelocityVerlet.hpp
 *
 * An implementation of the Velocity-Verlet integrator to use in
 * molecular dynamics simulations
 *
 * A subclass of the abstract Integrator which implements the
 * advance method according to Velocity-Verlet algorithm
 *
 *  Created on: Mar 31, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef VELOCITYVERLET_HPP_
#define VELOCITYVERLET_HPP_

#include "Integrator.hpp"
#include "Component.hpp"

namespace MolDyn {

template<class T>
class VelocityVerlet: public Integrator<T> {
private:
	T dt2;
	std::vector<T> halfDtM;
public:
	VelocityVerlet(Potential<T>& _potential, const unsigned int& _maxSpecies,
			const T& _tStep);
	void advance(std::vector<Vector3d<T> >& _positions,
			std::vector<Vector3d<T> >& _velocities,
			std::vector<Vector3d<T> >& _forces,
			const std::vector<unsigned int>& _species, T& _eKin, T& _ePot,
			T& _eTot, T& _pressure, const T& _length, bool rescale);
	/*
	 * Add a new species of particle to the integrator.
	 * A species is define by (mass,density,number of particles)
	 *
	 * @param[in]	_mass: The mass of the species
	 * @param[in]	_density: The density of the species
	 * @param[in]	_nParticles: The number of particles of such species
	 */
	void addSpecies(const T& _mass, const T& _density,
			const unsigned int& _nParticles);
	/*
	 * Add a new species of particle to the integrator.
	 * A species is define by (mass,density,number of particles)
	 *
	 * @param[in]	_component: The component to be added
	 * @see addSpecies(const T& _mass, const T& _density,
	 * 			const unsigned int& _nParticles)
	 */
	void addSpecies(const MolDyn::Component<T>& _component);
	virtual ~VelocityVerlet();
};

template<class T>
inline MolDyn::VelocityVerlet<T>::VelocityVerlet(Potential<T>& _potential,
		const unsigned int& _maxSpecies, const T& _tStep) :
		Integrator<T>(_potential, _maxSpecies, _tStep), dt2(0.5 * _tStep), halfDtM(
				_maxSpecies) {
}

template<class T>
inline void MolDyn::VelocityVerlet<T>::advance(
		std::vector<Vector3d<T> >& _positions,
		std::vector<Vector3d<T> >& _velocities,
		std::vector<Vector3d<T> >& _forces,
		const std::vector<unsigned int>& _species, T& _eKin, T& _ePot, T& _eTot,
		T& _pressure, const T& _length, bool rescale) {
	unsigned int nAtoms = _positions.size();
	/*
	 * ADVANCE POSITIONS FROM T TO T+DT AND VELOCITIES FROM T TO T+DT/2
	 * ACCORDING TO THE VELOCITY VERLET INTEGRATOR:
	 * 	V(T+DT/2) = V(T) + 1/2*DT*A(T)
	 * 	R(T+DT)   = R(T) + DT*V(T) + 1/2*(DT**2)*A(T)
	 * 	V(T+DT)   = V(T+DT/2) + 1/2*DT*A(T+DT)
	 */
	// ADVANCE POSITIONS BY DT AND VELOCITIES BY DT/2
	for (unsigned int i = 0; i < nAtoms; ++i) {
		_velocities[i] += _forces[i] * halfDtM[_species[i]];
		_positions[i] += _velocities[i] * Integrator<T>::dt;
	}

	// EVALUATE THE FORCES AT STEP T+DT
	Integrator<T>::forces(_positions, _forces, _species, _ePot, _pressure,
			_length);

	// ADVANCE VELOCITIES BY THE REMAINING DT/2, WITH NEW FORCES
	for (unsigned int i = 0; i < nAtoms; ++i) {
		_velocities[i] += _forces[i] * halfDtM[_species[i]];
	}

	// APPLY PERIODIC BOUNDARY CONDITIONS
	for (unsigned int i = 0; i < nAtoms; ++i) {
		_positions[i].modulo(_length);
	}

	if (rescale) {
		Integrator<T>::rescale(_velocities, _species);
	}

	Integrator<T>::currentEnergy(_velocities, _species, _eTot, _ePot, _eKin,
			_pressure);
}

template<class T>
inline void VelocityVerlet<T>::addSpecies(const T& _mass, const T& _density,
		const unsigned int& _nParticles) {
	Integrator<T>::addSpecies(_mass, _density, _nParticles);
	for (unsigned int i = 0; i < Integrator<T>::maxSpecies; ++i) {
		halfDtM[i] = dt2 / Integrator<T>::mass[i];
	}
}

template<class T>
inline void VelocityVerlet<T>::addSpecies(
		const MolDyn::Component<T>& _component) {
	addSpecies(_component.mass, _component.density, _component.nParticles);
}

template<class T>
inline MolDyn::VelocityVerlet<T>::~VelocityVerlet() {

}

} /* namespace MolDyn */

#endif /* VELOCITYVERLET_HPP_ */

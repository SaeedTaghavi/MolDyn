/*
 * Integrator.hpp
 *
 * Provides an abstract class for an integrator with minimal
 * intrinsic methods.
 *
 *  Created on: Mar 30, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_
#include "Potential.hpp"
#include <cmath>
#include "Vector3d.hpp"
#include <vector>
#include "constants.hpp"

namespace MolDyn {

template<class T>
class Integrator {
public:
	/*
	 * The pairwise potential used to calculate the forces
	 * and energy in the system
	 */
	Potential<T>& potential;

protected:
	/*
	 * The number of DIFFERENT species in the system
	 */
	unsigned int maxSpecies;

	/*
	 * A counter of species.
	 * Incremented each time a new species is added and used
	 * to limit the addition of more species to the maximum
	 * defined in the constructor
	 */
	unsigned int speciesCounter;

	/*
	 * The number of particles of a defined species
	 */
	std::vector<unsigned int> nParticles;

	/*
	 * The total number of particles
	 */
	unsigned int totalParticles;

	/*
	 * A vector with the different densities of the species
	 */
	std::vector<T> densities;

	/*
	 * The mass of a determine species
	 */
	std::vector<T> mass;

	/*
	 * The partial densities defined by
	 *
	 * partialDensity = Sum(density_i / n_particles_i)
	 */
	T sumDensity;

	/*
	 * The size of the time step used for the integration
	 */
	T dt;

	/*
	 * The initial temperature of the system
	 */
	T temperature;

public:
	/*
	 * Class constructor
	 * @param 	_potential: the pairwise potential
	 * @param	_maxSpecies: The number of DIFFERENT species
	 * @param	_tStep: The time step used for the integration
	 */
	Integrator(Potential<T>& _potential, const unsigned int& _maxSpecies,
			const T& _tStep);

	/*
	 * Calculate the force on each particle due to the rest of them
	 * using the pairwise potential.
	 *
	 * @param[in,out]	_positions: An array with the positions of the particles
	 * @param[in,out]	_forces: An array with the forces felt by each particle
	 * @param[in]		_species: An array of the species of each particles
	 * @param[in,out]	_ePot: The potential energy of the system
	 * @param[in,out]	_pressure: The pressure of the system
	 * @param[in]		_length: The length of the simulation box
	 */
	void forces(std::vector<Vector3d<T> >& _positions,
			std::vector<Vector3d<T> >& _forces,
			const std::vector<unsigned int>& _species, T& _ePot, T& _pressure,
			const T& _length);

	/*
	 * An interface providing a method to advance the system a time step
	 *
	 * This is a pure virtual method that should be implemented in the
	 * actual integrator (e.g. Velocity Verlet)
	 *
	 * @param[in,out]	_positions: An array with the positions of the particles
	 * @param[in,out]	_velocities: An array with the velocities of the particles
	 * @param[in,out]	_forces: An array with the forces felt by each particle
	 * @param[in]		_species: An array of the species of each particles
	 * @param[in,out]	_eKin: The kinetic energy of the system
	 * @param[in,out]	_ePot: The potential energy of the system
	 * @param[in,out]	_eTot: The total energy of the system
	 * @param[in,out]	_pressure: The pressure of the system
	 * @param[in]		_length: The length of the simulation box
	 * @param[in]		_rescale: A flag used to re-scale the velocities (default false)
	 */
	virtual void advance(std::vector<Vector3d<T> >& _positions,
			std::vector<Vector3d<T> >& _velocities,
			std::vector<Vector3d<T> >& _forces,
			const std::vector<unsigned int>& _species, T& _eKin, T& _ePot,
			T& _eTot, T& _pressure, const T& _length, bool _rescale) = 0;

	/*
	 * Determine the current total energy of the system and pressure
	 *
	 * @param[in]		_velocities: An array with the velocities of the particles
	 * @param[in]		_species: An array of the species of each particles
	 * @param[in,out]	_eTot: The total energy of the system
	 * @param[in,out]	_ePot: The potential energy of the system
	 * @param[in,out]	_eKin: The total energy of the system
	 * @param[in,out]	_pressure: The pressure of the system
	 */
	void currentEnergy(std::vector<Vector3d<T> >& _velocities,
			const std::vector<unsigned int>& _species, T& _eTot, T& _ePot,
			T& _eKin, T& _pressure);
	/*
	 * Re-scale the velocities according to the Gauss Constraint method
	 *
	 * @param[in]		_velocities: An array with the velocities of the particles
	 * @param[in,out]	_forces: An array with the forces felt by each particle
	 */
	void GaussConstraint(std::vector<Vector3d<T> >& _velocities,
			std::vector<Vector3d<T> >& _forces,
			std::vector<unsigned int>& _species);
	/*
	 * Re-scale the velocities of the particles according to the initial temperature
	 *
	 * @param[in,out]		_velocities: An array with the velocities of the particles
	 */
	void rescale(std::vector<Vector3d<T> >& _velocities,
			const std::vector<unsigned int>& _species);
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

	void addSpecies(const MolDyn::Component<T>& component);
	/*
	 * Get the total mass of the system
	 * @return	The total mass of the system
	 */
	T getTotalMass() const {
		T sumMass = 0;
		for (unsigned int i = 0; i < maxSpecies; ++i) {
			sumMass += mass[i] * nParticles[i];
		}
		return sumMass;
	}
	/*
	 * Set the temperature of the system
	 * @param[in]	_temperature: The temperature
	 */
	void setTemperature(const T& _temp) {
		temperature = _temp;
	}
	/*
	 * @return	The system's temperature
	 */
	T getTemp() const {
		return temperature;
	}
	/*
	 * The class destructor
	 */
	virtual ~Integrator();
};

template<class T>
inline MolDyn::Integrator<T>::Integrator(Potential<T>& _potential,
		const unsigned int& _maxSpecies, const T& _tStep) :
		potential(_potential), dt(_tStep), maxSpecies(_maxSpecies), densities(
				_maxSpecies), mass(_maxSpecies), nParticles(_maxSpecies), speciesCounter(
				0), sumDensity(0.0), totalParticles(0) {
}

template<class T>
inline void MolDyn::Integrator<T>::forces(std::vector<Vector3d<T> >& _positions,
		std::vector<Vector3d<T> >& _forces,
		const std::vector<unsigned int>& _species, T& _ePot, T& _pressure,
		const T& _length) {
	unsigned int _nAtoms = _positions.size();
	Vector3d<T> fi;
	Vector3d<T> rij;
	Vector3d<T> fij;
	T fijMagnitude;
	T rijLength;
	T _potentialEnergy;

	static unsigned int calls = 0;

	// Set forces on all atoms, potential energy and pressure to zero
	_ePot = 0;
	_pressure = 0;
	fijMagnitude = 0;
	for (unsigned int i = 0; i < _nAtoms; ++i) {
		_forces[i].zeros();
	}
	/*
	 * Loop over all pairwise interactions to obtain potential energy,
	 * configurational part of pressure, and forces
	 */
	for (unsigned int i = 0; i < _nAtoms - 1; ++i) {
		fi = _forces[i];
		for (unsigned int j = i + 1; j < _nAtoms; ++j) {
			rij = _positions[i] - _positions[j];
			// Apply minimun image convention
			rij.modulo(_length);
			rijLength = rij.length();

			/*
			 * If distance is less than LJ cutoff, calculate force according to
			 * gradient of potential energy and the configurational component
			 * of pressure by summing over rij*fij
			 */
			if (rijLength < potential.getCutOff(_species[i], _species[j])) {
				potential.evaluate(rijLength, _potentialEnergy, fijMagnitude,
						_species[i], _species[j]);
				_pressure += fijMagnitude * rijLength;
				_ePot += _potentialEnergy;
				fij = project(fijMagnitude, rij);
				fi += fij;
				_forces[j] -= fij;
			} // end if
		} // end for
		_forces[i] = fi;
	} // end for
	_pressure *= sumDensity / (_nAtoms * 3.0);
}

template<class T>
inline void MolDyn::Integrator<T>::currentEnergy(
		std::vector<Vector3d<T> >& _velocities,
		const std::vector<unsigned int>& _species, T& _eTot, T& _ePot, T& _eKin,
		T& _pressure) {
	_eKin = 0.0;
	for (unsigned int i = 0; i < _velocities.size(); ++i) {
		_eKin += mass[_species[i]] * (_velocities[i] * _velocities[i]);
	}
	_eKin *= 0.5;
	_eTot = _eKin + _ePot;
	_pressure += 2.0 * _eKin * sumDensity / (_velocities.size() * 3.0);
}

template<class T>
inline void MolDyn::Integrator<T>::GaussConstraint(std::vector<Vector3d<T> >& v,
		std::vector<Vector3d<T> >& f, std::vector<unsigned int>& _species) {
	T denSum = 0; // Sum( pi dot F_i / m_i)
	T numSum = 0; // Sum( pi dot P_i / m_i)
	T lambda = 0;
// Assuming m_i = mass (same for all particles)
	for (unsigned int i = 0; i < v.size(); ++i) {
		denSum += v[i] * v[i] / mass[_species[i]];
		numSum += v[i] * f[i] / mass[_species[i]];
	}
	lambda = (denSum == 0) ? 0.0 : numSum / denSum;
//	std::cout << "lambda = " << lambda << std::endl;

	for (unsigned int i = 0; i < v.size(); ++i) {
		v[i] *= lambda;
		f[i] -= v[i] * mass[_species[i]];
	}

}

template<class T>
inline void Integrator<T>::rescale(std::vector<Vector3d<T> >& _velocities,
		const std::vector<unsigned int>& _species) {
	T velsq, factor;
	velsq = 0.0;
	for (unsigned int i = 0; i < _velocities.size(); ++i) {
		velsq += mass[_species[i]] * (_velocities[i] * _velocities[i]);
	}
	factor = sqrt(3.0 * _velocities.size() * KB * temperature / velsq);
	for (unsigned int i = 0; i < _velocities.size(); ++i) {
		_velocities[i] *= factor;
	}
}

template<class T>
inline void Integrator<T>::addSpecies(const T& _mass, const T& _dens,
		const unsigned int& _nparticles) {
	if (speciesCounter < maxSpecies) {
		mass[speciesCounter] = _mass;
		densities[speciesCounter] = _dens;
		nParticles[speciesCounter] = _nparticles;
		totalParticles += _nparticles;
		sumDensity += _dens;
		++speciesCounter;
	}
}

template<class T>
inline void Integrator<T>::addSpecies(const MolDyn::Component<T>& component) {
	addSpecies(component.mass, component.density, component.nParticles);
}

template<class T>
inline MolDyn::Integrator<T>::~Integrator() {
	// Free memory
	std::vector<T>().swap(mass);
	std::vector<T>().swap(densities);
	std::vector<unsigned int>().swap(nParticles);
}

} /* namespace MolDyn */

#endif /* INTEGRATOR_HPP_ */

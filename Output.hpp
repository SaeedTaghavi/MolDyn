/*
 * Output.hpp
 *
 *  Created on: Apr 20, 2015
 *      Author: erick
 */

#ifndef OUTPUT_HPP_
#define OUTPUT_HPP_
// ios::exceptions
#include <iostream>
#include <fstream>
#include "Configuration.hpp"
#include "SimulationBox.hpp"
#include "Vector3d.hpp"

namespace MolDyn {

template<class T>
class Output {
private:
	MolDyn::Configuration<T>& configuration;
	std::ofstream filePositions;
	std::ofstream fileVelocities;
	std::ofstream fileMovie;
	std::ofstream fileEnergies;
	std::ofstream filePressure;
	std::ofstream fileDensityDistribution;
public:
	Output(MolDyn::Configuration<T>& _config);
	void write(const MolDyn::SimulationBox<T>& simulationBox,
			const unsigned int& step, bool isEquilibration);
	virtual ~Output();
};

template<class T>
inline MolDyn::Output<T>::Output(MolDyn::Configuration<T>& _config) :
		configuration(_config) {
	if (configuration.writeEnergies) {
		fileEnergies.open(configuration.strEnergiesFilename.c_str());
		if (!fileEnergies.is_open()) {
			throw std::ofstream::failure(
					"Error opening energies output file for writing");
		}
	}

	if (configuration.writePressure) {
		filePressure.open(configuration.strPressureFilename.c_str());
		if (!filePressure.is_open()) {
			throw std::ofstream::failure(
					"Error opening pressure output file for writing");
		}
	}

	if (configuration.writePositions) {
		filePositions.open(configuration.strTapePosFilename.c_str());
		if (!filePositions.is_open()) {
			throw std::ofstream::failure(
					"Error opening positions output file for writing");
		}
	}

	if (configuration.writeVelocities) {
		fileVelocities.open(configuration.strTapeVelFilename.c_str());
		if (!fileVelocities.is_open()) {
			throw std::ofstream::failure(
					"Error opening velocities output file for writing");
		}
	}

	if (configuration.writeMovie) {
		fileMovie.open(configuration.strMovieFilename.c_str());
		if (!fileMovie.is_open()) {
			throw std::ofstream::failure(
					"Error opening movie output file for writing");
		}
	}

	if (configuration.writeDensityDistribution) {
		fileDensityDistribution.open(
				configuration.strDensityDistFilename.c_str());
		if (!fileDensityDistribution.is_open()) {
			throw std::ofstream::failure(
					"Error opening Density Distribution output file for writing");
		}
	}
}

template<class T>
inline void MolDyn::Output<T>::write(
		const MolDyn::SimulationBox<T>& simulationBox, const unsigned int& step,
		bool isEquilibration = false) {
	// A vector to temporarily store simulation data
	MolDyn::Vector3d<T> V;
	/*
	 * Check if the step matches the energy output interval
	 * and if files are available, write the corresponding output
	 */
	if ((step + 1) % configuration.energyInterval == 0) {
		if (configuration.writeEnergies) {
			fileEnergies.setf(std::ios::fixed, std::ios::floatfield);
			fileEnergies.width(5);
			fileEnergies << step + 1 << "\t ";
			fileEnergies.precision(15);
			fileEnergies.width(30);
			fileEnergies << simulationBox.getEKin() << "\t ";
			fileEnergies.width(30);
			fileEnergies << simulationBox.getEPot() << "\t ";
			fileEnergies.width(30);
			fileEnergies << simulationBox.getETot() << std::endl;
		} // end if (config.writeEnergies)

		if (configuration.writePressure) {
			filePressure.setf(std::ios::fixed, std::ios::floatfield);
			filePressure.width(10);
			filePressure << step + 1 << "\t ";
			filePressure.precision(15);
			filePressure.width(25);
			filePressure << simulationBox.getPressure() << std::endl;
		} // end if (config.writePressure)
	} // end if ((step + 1) % config.energyInterval == 0)

	/*
	 * If the call is in production mode and the configuration provides
	 * output files, write the animation, position and velocities of the
	 * simulation
	 */
	if (!isEquilibration) {
		if ((step + 1) % configuration.animationInterval == 0
				&& configuration.writeMovie) {
			fileMovie << configuration.totalParticles << std::endl << std::endl;
			for (unsigned int i = 0; i < configuration.totalParticles; ++i) {
				V = simulationBox.getPositions().at(i);
				fileMovie.setf(std::ios::fixed, std::ios::floatfield);
				fileMovie
						<< configuration.components.at(
								simulationBox.getSpecies().at(i)).identifier;
				fileMovie << "    ";
				fileMovie.precision(4);
				fileMovie.width(10);
				fileMovie << V[0] << "\t ";
				fileMovie.width(10);
				fileMovie << V[1] << "\t ";
				fileMovie.width(10);
				fileMovie << V[2] << std::endl;
			}
		} // end if ((step + 1) % config.animationInterval == 0 && config.writeMovie)

		if ((step + 1) % configuration.tapeInterval == 0) {
			for (unsigned int i = 0; i < configuration.totalParticles; ++i) {
				if (configuration.writeVelocities) {
					V = simulationBox.getVelocities().at(i);
					fileVelocities.precision(20);
					fileVelocities.setf(std::ios::scientific,
							std::ios::floatfield);
					fileVelocities.width(25);
					fileVelocities << V[0] << "\t ";
					fileVelocities.width(25);
					fileVelocities << V[1] << "\t ";
					fileVelocities.width(25);
					fileVelocities << V[2] << std::endl;
				} // end if (config.writeVelocities)
				if (configuration.writePositions) {
					V = simulationBox.getPositions().at(i);
					filePositions.precision(20);
					filePositions.setf(std::ios::scientific,
							std::ios::floatfield);
					filePositions.width(25);
					filePositions << V[0] << "\t ";
					filePositions.width(25);
					filePositions << V[1] << "\t ";
					filePositions.width(25);
					filePositions << V[2] << std::endl;
				} // end if (config.writePositions)
			} // end for
		} // end if ((step + 1) % config.tapeInterval == 0)
	} // end if (!isEquilibration)

}

template<class T>
inline MolDyn::Output<T>::~Output() {
	// Close all open streams
	if (configuration.writeEnergies) {
		fileEnergies.close();
	}
	if (configuration.writePressure) {
		filePressure.close();
	}
	if (configuration.writePositions) {
		filePositions.close();
	}
	if (configuration.writeVelocities) {
		fileVelocities.close();
	}
	if (configuration.writeMovie) {
		fileMovie.close();
	}
	if (configuration.writeDensityDistribution) {
		fileDensityDistribution.close();
	}
}

} /* namespace MolDyn */

#endif /* OUTPUT_HPP_ */

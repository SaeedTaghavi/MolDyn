/*
 * main.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author: Erick Martinez
 */
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <exception>
#include "Configuration.hpp"
#include "Output.hpp"
#include "MixedLJPotential.hpp"
#include "VelocityVerlet.hpp"
#include "SimulationBox.hpp"
#include "Component.hpp"
#include "constants.hpp"

#define MAXSHELLS 100

typedef double real;

int main(int argc, char **argv) {
	// The name of the configuration file
	std::string strConfigurationFilename;
	std::vector<MolDyn::Component<real> > components;
	real simBoxLength;

	/* Command line initialization */
	if (argc < 3) {
		std::cout << "Usage is moldyn -i <config filename>" << std::endl;
		return 0;
	}

	for (int i = 1; i < argc; ++i) {
		if (i + 1 != argc) {
			if (std::string(argv[i]).compare("-i") == 0) {
				strConfigurationFilename = argv[i + 1];
			}
		}
	}

	if (strConfigurationFilename.length() < 3) {
		std::cout << "No input files found." << std::endl;
		std::cout << "Usage is moldyn -i <config filename>" << std::endl;
		return 0;
	}

	MolDyn::Configuration<real> config(strConfigurationFilename);
	MolDyn::Output<real> output(config);

	std::cout << "----------- VARIABLES IN FILE ------------- " << std::endl;
	std::cout << "Number of particles: " << config.totalParticles << std::endl;
	std::cout << "LJCUT: " << config.LJCutoff << std::endl;
	std::cout << "Temperature (K): " << config.temperature << std::endl;
	std::cout << "dt (ps): " << config.dt << std::endl;
	std::cout << "Equilibration steps: " << config.equilibrationSteps
			<< std::endl;
	std::cout << "Production steps: " << config.productionSteps << std::endl;
	std::cout << "Energy interval (steps): " << config.energyInterval
			<< std::endl;
	std::cout << "Tape interval (steps): " << config.tapeInterval << std::endl;
	std::cout << "Animation interval (steps): " << config.animationInterval
			<< std::endl;
	std::cout << "Number of species: " << config.maxComponents << std::endl;
	for (unsigned int i = 0; i < config.maxComponents; ++i) {
		std::cout << "Component " << i + 1 << " out of "
				<< config.maxComponents;
		std::cout << std::endl;
		std::cout << "Sigma (A): " << config.components[i].sigma << std::endl;
		std::cout << "Epsilon/KB: " << config.components[i].epsilonKB
				<< std::endl;
		std::cout << "Epsilon: " << config.components[i].epsilon << std::endl;
		std::cout << "Mass (amu): " << config.components[i].mass << std::endl;
		std::cout << "Density (gm/cm3): " << config.components[i].densityGCM
				<< std::endl;
		std::cout << "Density (particles/A3): " << config.components[i].density
				<< std::endl;
	}

	simBoxLength = MolDyn::getSimBoxLength(config.components);
	MolDyn::MixedLJPotential<real> potential(config.maxComponents,
			config.LJCutoff);
	for (unsigned int i = 0; i < config.maxComponents; ++i) {
		potential.addComponent(config.components[i]);
	}
	MolDyn::VelocityVerlet<real> integrator(potential, config.maxComponents,
			config.dt);

	for (unsigned int i = 0; i < config.maxComponents; ++i) {
		integrator.addSpecies(config.components[i]);
	}

	MolDyn::SimulationBox<real> simulationBox(config.totalParticles, integrator,
			config.temperature, simBoxLength, config.sigmaCut, config.maxComponents);
	for (unsigned int i = 0; i < config.maxComponents; ++i) {
		simulationBox.addComponent(config.components[i]);
	}

	simulationBox.initialize();
	std::cout << "STARTING EQUILIBRATION RUN" << std::endl;

	for (unsigned int i = 0; i < config.equilibrationSteps; ++i) {
		simulationBox.step(true);
		if ((i + 1) % 1000 == 0) {
			std::cout << "Equilibration step " << i + 1 << std::endl;
		}
		output.write(simulationBox, i, true);
	}

	std::cout << "EQUILIBRATION OVER" << std::endl;
	std::cout << "STARTING PRODUCTION RUN" << std::endl;

	for (unsigned int i = 0; i < config.productionSteps; ++i) {
		simulationBox.step(false);
		if ((i + 1) % 1000 == 0) {
			std::cout << "Production step " << i + 1 << std::endl;
		}
		output.write(simulationBox, i + config.equilibrationSteps, false);
	} // end for

	std::cout << "PRODUCTION OVER" << std::endl;

}


/*
 * main.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author: Erick Martinez
 */
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "LJPotential.hpp"
#include "VelocityVerlet.hpp"
#include "SimulationBox.hpp"
#include "constants.hpp"

#define SCUTOFF 0.9

typedef double real;

int main(int argc, char **argv) {
	// The property tree read from the configuration file
	boost::property_tree::ptree pt;
	// The name of the configuration file
	std::string str_config_fn;
	std::string str_tape_vel_fn;
	std::string str_tape_pos_fn;
	std::string str_movie_fn;
	std::string str_energies_fn;
	std::string str_pressure_fn;
	std::ofstream file_tape_pos;
	std::ofstream file_tape_vel;
	std::ofstream file_movie;
	std::ofstream file_energies;
	std::ofstream file_pressure;
	bool write_positions = false;
	bool write_velocities = false;
	bool write_movie = false;
	bool write_energies = false;
	bool write_pressure = false;
	unsigned int output_file_cnt = 0;
	MolDyn::Vector3d<real> V;

	unsigned int NATOMS;
	real EPSILON;
	real SIGMA;
	real MASS;
	real DENS;
	real LJCUT;
	real TEMP;
	real DT;
	real LBOX;
	unsigned int EQSTEPS;
	unsigned int PRODSTEPS;
	unsigned int IE;
	unsigned int ITAPE;
	unsigned int IANIMATE;

	/* Command line initialization */
	if (argc < 3) {
		std::cout << "Usage is moldyn -i <config filename>" << std::endl;
		return 0;
	}

	for (int i = 1; i < argc; ++i) {
		if (i + 1 != argc) {
			if (std::string(argv[i]).compare("-i") == 0) {
				str_config_fn = argv[i + 1];
			}
		}
	}

	if (str_config_fn.length() < 3) {
		std::cout << "No input files found." << std::endl;
		std::cout << "Usage is moldyn -i <config filename>" << std::endl;
		return 0;
	}

	try {
		boost::property_tree::ini_parser::read_ini(str_config_fn.c_str(), pt);
		NATOMS = pt.get<unsigned int>("Params.n_atoms");
		EPSILON = pt.get<real>("Params.epsilon");
		SIGMA = pt.get<real>("Params.sigma");
		MASS = pt.get<real>("Params.mass");
		DENS = pt.get<real>("Params.density");
		LJCUT = pt.get<real>("Params.LJ_cutoff");
		TEMP = pt.get<real>("Params.temp");
		DT = pt.get<real>("Params.dt");
		EQSTEPS = pt.get<unsigned int>("Params.eq_steps");
		PRODSTEPS = pt.get<unsigned int>("Params.prod_steps");
		IE = pt.get<unsigned int>("Params.i_energy");
		ITAPE = pt.get<unsigned int>("Params.i_tape");
		IANIMATE = pt.get<unsigned int>("Params.i_animate");
	} catch (std::exception& e) {
		std::cout << "Error reading the configuration file '" << str_config_fn
				<< "'." << std::endl;
		std::cout << e.what() << std::endl;
		return 0;
	}

	try {
		str_tape_vel_fn = pt.get<std::string>("Files.tape_vel");
		write_velocities = true;
		++output_file_cnt;

	} catch (std::exception& e) {
		std::cout << "Did not find 'Files.tape_vel'" << std::endl;
	}

	try {
		str_tape_pos_fn = pt.get<std::string>("Files.tape_pos");
		write_positions = true;
		++output_file_cnt;

	} catch (std::exception& e) {
		std::cout << "Did not find 'Files.tape_pos'" << std::endl;
	}

	try {
		str_movie_fn = pt.get<std::string>("Files.movie");
		write_movie = true;
		++output_file_cnt;

	} catch (std::exception& e) {
		std::cout << "Did not find 'Files.movie'" << std::endl;
	}

	try {
		str_energies_fn = pt.get<std::string>("Files.energies");
		write_energies = true;
		++output_file_cnt;

	} catch (std::exception& e) {
		std::cout << "Did not find 'Files.energies'" << std::endl;
	}

	try {
		str_pressure_fn = pt.get<std::string>("Files.pressure");
		write_pressure = true;
		++output_file_cnt;

	} catch (std::exception& e) {
		std::cout << "Did not find 'Files.pressure'" << std::endl;
	}

	if (output_file_cnt < 1) {
		std::cout << "No output files specified." << std::endl;
		std::cout << "Program will stop now." << std::endl;
		return 0;
	}

	if (write_positions) {
		file_tape_pos.open(str_tape_pos_fn.c_str());
		if (file_tape_pos.is_open() == false) {
			std::cout << "There was an error opening the positions tape file"
					<< std::endl;
			std::cout << "File: '" << str_tape_pos_fn << "'" << std::endl;
			return 0;
		}
	}

	if (write_velocities) {
		file_tape_vel.open(str_tape_vel_fn.c_str());
		if (file_tape_vel.is_open() == false) {
			std::cout << "There was an error opening the velocities tape file"
					<< std::endl;
			std::cout << "File: '" << str_tape_vel_fn << "'" << std::endl;
			return 0;
		}
	}

	if (write_movie) {
		file_movie.open(str_movie_fn.c_str());
		if (file_movie.is_open() == false) {
			std::cout << "There was an error opening the movie file"
					<< std::endl;
			std::cout << "File: '" << str_movie_fn << "'" << std::endl;
			return 0;
		}
	}

	if (write_energies) {
		file_energies.open(str_energies_fn.c_str());
		if (file_energies.is_open() == false) {
			std::cout << "There was an error opening the movie file"
					<< std::endl;
			std::cout << "File: '" << str_energies_fn << "'" << std::endl;
			return 0;
		}
	}

	if (write_pressure) {
		file_pressure.open(str_pressure_fn.c_str());
		if (file_pressure.is_open() == false) {
			std::cout << "There was an error opening the pressure file"
					<< std::endl;
			std::cout << "File: '" << str_pressure_fn << "'" << std::endl;
			return 0;
		}
	}

	std::cout << "----------- VARIABLES IN FILE ------------- " << std::endl;
	std::cout << "Number of atoms: " << NATOMS << std::endl;
	std::cout << "Density (gm/cm3): " << DENS << std::endl;
	std::cout << "Epsilon/KB: " << EPSILON << std::endl;
	std::cout << "Sigma (A): " << SIGMA << std::endl;
	std::cout << "Mass (amu): " << MASS << std::endl;
	std::cout << "LJCUT: " << LJCUT << std::endl;
	std::cout << "Temperature (K): " << TEMP << std::endl;
	std::cout << "dt (ps): " << DT << std::endl;
	std::cout << "Equilibration steps: " << EQSTEPS << std::endl;
	std::cout << "Production steps: " << PRODSTEPS << std::endl;
	std::cout << "Energy interval (steps): " << IE << std::endl;
	std::cout << "Tape interval (steps): " << ITAPE << std::endl;
	std::cout << "Animation interval (steps): " << IANIMATE << std::endl;

	EPSILON *= MolDyn::KB;

	DENS = MolDyn::atomicDensity(DENS, MASS);
	std::cout << "Density (Atoms/ps^3): " << DENS << std::endl;
	LBOX = MolDyn::getSimBoxLength(NATOMS, DENS);
	MolDyn::LJPotential<real> potential(SIGMA, EPSILON);
	potential.setCutOff(LJCUT);

	MolDyn::VelocityVerlet<real> integrator(potential, MASS, DENS, DT);
	MolDyn::SimulationBox<real> simulationBox(NATOMS, DENS, integrator, TEMP,
			LBOX, SCUTOFF);

	simulationBox.initialize();
	std::cout << "STARTING EQUILIBRATION RUN" << std::endl;

	for (unsigned int i = 0; i < EQSTEPS; ++i) {
		simulationBox.step(true);
		if ((i + 1) % 1000 == 0) {
			std::cout << "Equilibration step " << i + 1 << std::endl;
		}
		if ((i + 1) % IE == 0) {
			if (write_energies) {
				file_energies.setf(std::ios::fixed, std::ios::floatfield);
				file_energies.width(5);
				file_energies << i + 1 << "\t ";
				file_energies.precision(15);
				file_energies.width(30);
				file_energies << simulationBox.getEKin() << "\t ";
				file_energies.width(30);
				file_energies << simulationBox.getEPot() << "\t ";
				file_energies.width(30);
				file_energies << simulationBox.getETot() << std::endl;
			}

			if (write_pressure) {
				file_pressure.setf(std::ios::fixed, std::ios::floatfield);
				file_pressure.width(10);
				file_pressure << i + 1 << "\t ";
				file_pressure.precision(15);
				file_pressure.width(25);
				file_pressure << simulationBox.getPressure() << std::endl;
			}

		} // end if ((i + 1) % IE == 0)
	}

	for (unsigned int i = 0; i < PRODSTEPS; ++i) {
		simulationBox.step(false);
		if ((i + 1) % 1000 == 0) {
			std::cout << "Production step " << i + 1 << std::endl;
		}

		if ((i + 1) % ITAPE == 0) {
			for (unsigned int j = 0; j < NATOMS; ++j) {
				if (write_velocities) {
					V = simulationBox.getVelocities().at(j);
					file_tape_vel.precision(20);
					file_tape_vel.setf(std::ios::scientific, std::ios::floatfield);
					file_tape_vel.width(25);
					file_tape_vel << V[0] << "\t ";
					file_tape_vel.width(25);
					file_tape_vel << V[1] << "\t ";
					file_tape_vel.width(25);
					file_tape_vel << V[2] << std::endl;

				}

				if (write_positions) {
					V = simulationBox.getPositions().at(j);
					file_tape_pos.precision(20);
					file_tape_pos.setf(std::ios::scientific, std::ios::floatfield);
					file_tape_pos.width(25);
					file_tape_pos << V[0] << "\t ";
					file_tape_pos.width(25);
					file_tape_pos << V[1] << "\t ";
					file_tape_pos.width(25);
					file_tape_pos << V[2] << std::endl;

				}
			} // end for (unsigned int j = 0; j < NATOMS; ++j)
		} // end if (i % ITAPE == 0)

		if ((i + 1) % IANIMATE == 0 && write_movie) {
			file_movie << NATOMS << std::endl << std::endl;
			for (unsigned int j = 0; j < NATOMS; ++j) {
				V = simulationBox.getPositions().at(j);
				file_movie.setf(std::ios::scientific, std::ios::floatfield);
				file_movie << "H    ";
				file_movie.precision(4);
				file_movie.width(10);
				file_movie << V[0] << "\t ";
				file_movie.width(10);
				file_movie << V[1] << "\t ";
				file_movie.width(10);
				file_movie << V[2] << std::endl;

			}
		} // end if (i % IANIMATE == 0 && write_movie)

		if ((i + 1) % IE == 0) {
			if (write_energies) {
				file_energies.setf(std::ios::scientific, std::ios::floatfield);
				file_energies.width(5);
				file_energies << i + 1  + EQSTEPS << "\t ";
				file_energies.precision(15);
				file_energies.width(30);
				file_energies << simulationBox.getEKin() << "\t ";
				file_energies.width(30);
				file_energies << simulationBox.getEPot() << "\t ";
				file_energies.width(30);
				file_energies << simulationBox.getETot() << std::endl;
			}

			if (write_pressure) {
				file_pressure.setf(std::ios::scientific, std::ios::floatfield);
				file_pressure.width(5);
				file_pressure << i + 1 + EQSTEPS << "\t ";
				file_pressure.precision(15);
				file_pressure.width(25);
				file_pressure << simulationBox.getPressure() << std::endl;
			}

		} // end if ((i + 1) % IE == 0)

	} // end for
	if (write_velocities) {
		file_tape_vel.close();
	}

	if (write_positions) {
		file_tape_pos.close();
	}

	if (write_movie) {
		file_movie.close();
	}

	if (write_energies) {
		file_energies.close();
	}

	if (write_pressure) {
		file_pressure.close();
	}

}


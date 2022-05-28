#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>
#include <stdio.h>
#include <mpir.h>
#include <mpirxx.h>
#include <unordered_map>
//#include <kfr/all.hpp>
//#include <kfr/dft.hpp>

#include "fftw3.h"
#include "minkowski_sum.h"
#include "bringmann_ssum.h"
#include "bellman_ssum.h"
#include "n_choose_k.h"
#include "koiliaris_xu_ssum.h"

std::pair<bool, double> test_algorithm(const std::vector<int>& set, int target, int num_runs, std::string alg_name, double delta = 0.1) {

	std::pair<bool, double> solution = {};
	double mean_time = 0;
	for (int run = 0; run < num_runs; run++) {
		if (alg_name == "bellman") {
			solution = bellman_ssum(set, target);
		}
		else if (alg_name == "koiliaris_xu") {
			solution = koiliaris_xu_ssum(set, target);
		}
		else if (alg_name == "bringmann") {
			solution = bringmann_ssum(set, target, delta);
		}
		mean_time += solution.second;
	}

	mean_time = mean_time / num_runs;

	return std::pair<bool, double>({solution.first, mean_time });
}

// SETS OF INTEGERS MUST BE SORTED!
int main(int argc, char *argv[]) {

	std::cout << "---------------------------------------------------" << std::endl;
	if (argc == 0) {
		std::cout << "Input file is not specified! " << std::endl;
		return 1;
	}
	else if (argc == 1) {
		std::cout << "Output file is not specified! " << std::endl;
		return 1;
	}
	else if (argc == 2) {
		std::cout << "Algorithm is not specified! " << std::endl;
		return 1;
	}
	else if (argc == 3) {
		std::cout << "Number of experiments is not specified! " << std::endl;
		return 1;
	}
	std::string in_filename = argv[1];
	std::string out_filename = argv[2];
	std::string algorithm = argv[3];
	int num_runs = std::stoi(argv[4]);
	std::vector<int> set;
	int target = 0;
	
	std::ifstream fin;
	std::ofstream fout;
	fin.open(in_filename);
	fout.open(out_filename);
	if (fin.is_open() && fout.is_open()) {
		std::string line;
		/* read the first line with a set */
		getline(fin, line);
		std::stringstream text_stream(line);
		std::string number;
		while (std::getline(text_stream, number, ' ')) {
			set.push_back(std::stoi(number));
		}
		/* read second line with target integer (unsafe)*/
		getline(fin, line);
		target = std::stoi(line);

		std::cout << "Input set size: " << set.size() << std::endl;
		std::cout << "Target: " << target << std::endl << std::endl;
	
		/* Check algorithm name */
		if (algorithm == "test") { // for testing
			return 0;
		}
		else if (algorithm == "colorcodinglayer") {
			std::cout << "--- Running ColorCodingLayer function ---" << std::endl;

			double mean_time = 0;
			for (int run = 0; run < num_runs; run++) {
				auto start_time = std::chrono::high_resolution_clock::now();
				ColorCodingLayer(set, target, 16, 0.1);
				auto end_time = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> running_time = end_time - start_time;
				mean_time += running_time.count();
			}
			std::cout << mean_time / num_runs << std::endl;

			fout << 0 << std::endl;
			fout << mean_time / num_runs << std::endl;
		}
		else if (algorithm == "colorcoding") {
			std::cout << "--- Running ColorCoding function ---" << std::endl;

			double mean_time = 0;
			for (int run = 0; run < num_runs; run++) {
				auto start_time = std::chrono::high_resolution_clock::now();
				ColorCoding(set, target, 10, 0.1);
				auto end_time = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> running_time = end_time - start_time;
				mean_time += running_time.count();
			}
			std::cout << mean_time / num_runs << std::endl;

			fout << 0 << std::endl;
			fout << mean_time / num_runs << std::endl;
		}

		else if (algorithm == "minkowski") {
			std::cout << "--- Running Minkowski sum ---" << std::endl;

			double mean_time = 0;
			for (int run = 0; run < num_runs; run++) {
				auto start_time = std::chrono::high_resolution_clock::now();
				//minkowski_sum_mpir(set, set, set.back());
				minkowski_sum_fft(set, set, set.back());
				
				auto end_time = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> running_time = end_time - start_time;
				mean_time += running_time.count();
				
			}
			std::cout <<  mean_time/num_runs << std::endl;

			fout << 0 << std::endl;
			fout << mean_time/ num_runs << std::endl;
		}
		else if (algorithm == "bellman") {
			std::cout << "--- Running Bellman's algorithm ---" << std::endl;
			auto bellman_solution = test_algorithm(set, target, num_runs, "bellman");
			std::cout << "Bellman solution: " << bellman_solution.first << ", time: " << bellman_solution.second << " ms" << std::endl;
			
			fout << bellman_solution.first << std::endl;
			fout << bellman_solution.second << std::endl;

			return 0;
		}
		else if (algorithm == "koiliaris_xu") {
			std::cout << "--- running koiliaris & xu's algorithm --- " << std::endl;
			auto koiliaris_xu_solution = test_algorithm(set, target, num_runs, "koiliaris_xu");
			std::cout << "koiliaris & xu solution: " << koiliaris_xu_solution.first << ", time: " << koiliaris_xu_solution.second << " ms" << std::endl;
			fout << koiliaris_xu_solution.first << std::endl;
			fout << koiliaris_xu_solution.second << std::endl;

			return 0;
		}
		else if (algorithm == "bringmann") {
			std::cout << "--- running bringmann's algorithm ---" << std::endl;
			srand(time(0));
			auto bringmann_solution = test_algorithm(set, target, num_runs, "bringmann", 0.1);
			std::cout << "bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << " ms" << std::endl;
			fout << bringmann_solution.first << std::endl;
			fout << bringmann_solution.second << std::endl;
			return 0;
		}
		else {
			std::cout << "No such algorithm! " << std::endl;
			return 1;
		}
	}
	else {
		std::cout << "File " << in_filename << " or " << out_filename << " was not found!" << std::endl;
		return 1;
	}
	return 0;
}
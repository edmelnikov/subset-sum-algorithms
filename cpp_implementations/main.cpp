#include <iostream>
#include <chrono>
#include "fftw3.h"
#include "minkowski_sum.h"
#include "bringmann_ssum.h"
#include "bellman_ssum.h"
#include "n_choose_k.h"
#include "koiliaris_xu_ssum.h"

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>

std::pair<bool, double> test_algorithm(const std::vector<int>& set, int target, int num_runs, std::string alg_name, double delta = 0.1) {

	std::pair<bool, double> solution;
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

	//std::vector<std::pair<int, int>> set1 = { {3, 1}, {4, 2}, {5, 2} };
	//std::vector<std::pair<int, int>> set2 = { {1, 1}, {2, 2}, {3, 3} };
	//	std::vector<std::pair<int, int>> set = minkowski_add_2d(set1, set2);
	//for (auto el : set) {
	//	std::cout << "(" << el.first << ", " << el.second << ") ";
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	if (argc == 0) {
		std::cout << "Input file is not specified! " << std::endl;
		return 1;
	}
	else if (argc == 1) {
		std::cout << "Algorithm is not specified! " << std::endl;
		return 1;
	}
	else if (argc == 2) {
		std::cout << "Number of experiments is not specified! " << std::endl;
		return 1;
	}
	std::string filename = argv[1];
	std::string algorithm = argv[2];
	int num_runs = std::stoi(argv[3]);
	std::vector<int> set;
	int target = 0;
	
	std::ifstream fin;
	fin.open(filename);
	if (fin.is_open()) {
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
		//for (auto el : set) {
		//	std::cout << el << " ";
		//}
		std::cout << "Target: " << target << std::endl << std::endl;

		/* Check algorithm name */
		if (algorithm == "all") {
			std::cout << "--- Running all algorithms ---" << std::endl;

			auto bellman_solution = test_algorithm(set, target, num_runs, "bellman");
			std::cout << "Bellman solution: " << bellman_solution.first << ", time: " << bellman_solution.second << " ms" << std::endl;
			bellman_solution = bellman_ssum(set, target);

			auto koiliaris_xu_solution = test_algorithm(set, target, num_runs, "koiliaris_xu");
			std::cout << "Koiliaris & Xu solution: " << koiliaris_xu_solution.first << ", time: " << koiliaris_xu_solution.second << " ms" << std::endl;
			
			srand(time(0));
			auto bringmann_solution = test_algorithm(set, target, num_runs, "bringmann", 0.2);
			std::cout << "Bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << " ms" << std::endl;
		
			return 0;
		}
		else if (algorithm == "bellman") {
			std::cout << "--- Running Bellman's algorithm ---" << std::endl;
			auto bellman_solution = test_algorithm(set, target, num_runs, "bellman");
			std::cout << "Bellman solution: " << bellman_solution.first << ", time: " << bellman_solution.second << " ms" << std::endl;

			return 0;
		}
		else if (algorithm == "koiliaris_xu") {
			std::cout << "--- Running Koiliaris & Xu's algorithm --- " << std::endl;
			auto koiliaris_xu_solution = test_algorithm(set, target, num_runs, "koiliaris_xu");
			std::cout << "Koiliaris & Xu solution: " << koiliaris_xu_solution.first << ", time: " << koiliaris_xu_solution.second << " ms" << std::endl;

			return 0;
		}
		else if (algorithm == "bringmann") {
			std::cout << "--- Running Bringmann's algorithm ---" << std::endl;
			srand(time(0));
			auto bringmann_solution = test_algorithm(set, target, num_runs, "bringmann", 0.1);
			std::cout << "Bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << " ms" << std::endl;

			return 0;
		}
		else {
			std::cout << "No such algorithm! " << std::endl;
			return 1;
		}
	}
	else {
		std::cout << "File " << filename << " was not found!" << std::endl;
		return 1;
	}


	//auto bellman_solution = bellman_ssum(set, target);
	//std::cout << "Bellman solution: " << bellman_solution.first << ", time: " << bellman_solution.second << " ms" << std::endl;

	//auto koiliaris_xu_solution = koiliaris_xu_ssum(set, target); 
	//std::cout << "Koiliaris & Xu solution: " << koiliaris_xu_solution.first << ", time: " << koiliaris_xu_solution.second << " ms" << std::endl;

	//srand(time(0));
	//auto bringmann_solution = bringmann_ssum(set, target, 0.1);
	//std::cout << "Bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << " ms" << std::endl;

	return 0;
}
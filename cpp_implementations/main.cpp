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

// SETS OF INTEGERS MUST BE SORTED!
void main() {
	// nine at the end?
	//std::vector<int> poly1 = { 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0 };
	//std::vector<int> poly2 = { 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
	////std::vector<int> poly2 = { 1, 0, 0, 0, 1, 1, 1 }; //  padding problem
	//std::vector<int> result = fft_polymul(poly1, poly2);
	//for (int i = 0; i < result.size(); i++) {
	//	std::cout << result[i] << ", ";
	//}

	//std::vector<std::vector<int>> mat1 = {
	//	{1, 1},
	//	{1, 1},
	//	{0, 3}
	//};
	//std::vector<std::vector<int>> mat2 = {
	//	{1, 0},
	//	{0, 1},
	//};

	//std::vector<std::vector<int>> res = fft_polymul_2d(mat1, mat2);
	//for (auto set : res) {
	//	for (auto el : set) {
	//		std::cout << el << " ";
	//	}
	//	std::cout << std::endl;
	//}

	//std::vector<std::pair<int, int>> mat1 = {{1, 1}, {2, 3}};
	//std::vector<std::pair<int, int>> mat2 = {{4, 2}, {2, 1}};
	//
	//for (auto pair : mat1) {
	//	std::cout << pair.first << std::endl;
	//}


	//std::vector<std::pair<int, int>> res = minkowski_add_2d(mat1, mat2);
	//
	//for (auto pair : res) {
	//	std::cout << "(" << pair.first << ", " << pair.second << ") ";
	//}
	//std::cout << std::endl;

	//std::vector<std::pair<int, int>> set1 = { {5, 1 }, {6, 1}, {11, 2} };
	//std::vector<std::pair<int, int>> set2 = { {7, 1}, {9, 1}, {16, 2} };
	//std::vector<std::pair<int, int>> res = ssum_card_set_union(set1, set2, 2, 5, 9);

	//for (auto el : res) {
	//	std::cout << "(" << el.first << ", " << el.second << ") ";
	//}
	//std::cout << std::endl;

	std::vector<int> st = {5, 6, 7, 9};

	//std::vector<std::pair<int, int>> res = compute_ssum_card_set(st, 4);

	std::vector<int> res = compute_ssums_in_range(st, 100);

	for (auto el : res) {
		std::cout << el << " ";
	}
	std::cout << std::endl;
	//for (auto el : res) {
	//	std::cout << "(" << el.first << ", " << el.second << ") ";
	//}
	//std::cout << std::endl;


	/* Minkowksy add test */
	//std::vector<int> set1 = { 2, 8, 20 };
	//std::vector<int> set2 = {3, 4, 15 };
	//
	///*std::vector<int> set1 = {1 };
	//std::vector<int> set2 = { 2};*/

	//int bound = 17;
	//std::vector<int> prod = minkowski_add(set1, set2, bound);
	//// std::cout << prod.size() << std::endl;
	//for (auto element : prod) {
	//	std::cout << element << ' ';
	//}
	//std::cout << std::endl;

	// std::vector<int> set = { 1, 5, 8, 9, 10, 12, 14, 15, 16, 18, 27, 30 };
	//std::vector<int> set = { 1, 2, 16, 17, 18, 27, 30 };
	//std::vector<int> set = { 3, 9, 15, 23 };
	
	//int target = 27;

	std::vector<int> set = { 2, 6, 7, 10, 12, 17, 18, 19, 21, 35, 36, 70, 71 };
	int target = 69; 
	//std::vector<int> set = { 3, 9, 15, 23 };
	//int target = 70;
	// in: 3, 9, 15, 23
	// out: 3 9 12 15 18 23 24 26 27 32 35 38 41 47 50

	auto bellman_solution = bellman_ssum(set, target);
	std::cout << "Bellman solution: " << bellman_solution.first << ", time: " << bellman_solution.second << " ms" << std::endl;
	

	std::cout << "Bringmann stuff " << std::endl;
	srand(time(0));
	auto bringmann_solution = bringmann_ssum(set, target, 0.1);
	std::cout << "Bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << " ms" << std::endl;

	// std::vector<int> res = ColorCoding(set, target, 4, 0.1);

	//for (auto el : res) {
	//	std::cout << el << " ";
	//}
	//std::cout << std::endl;
	
	//std::vector<int> set_ = { 4, 5, 6, 10 };
	//target = 10;
	//std::vector<int> ssums = ColorCodingLayer(set_, target, 16, 0.025);

	//for (auto el : ssums) {
	//	std::cout << el << std::endl;
	//}
	//std::cout << std::endl;

	
	
	//auto bringmann_solution = bringmann_ssum(set, target);
	//std::cout << "Bringmann solution: " << bringmann_solution.first << ", time: " << bringmann_solution.second << std::endl;

	//int k = 4;
	//std::vector<std::vector<int>> part = rand_part(set, 4, 42);

	//for (auto subset : part) {
	//	std::cout << subset.size() << std::endl;
	//	for (auto element : subset) {
	//		std::cout << element << " ";
	//	}
	//	std::cout << std::endl;
	//}


	//std::vector<int> set1 = {1, 2};
	//std::vector<int> set2 = { 4, 10 };
	//std::vector<int> set3 = { 5, 12 };
	//std::vector<int> un;

	//un.reserve(10);
	//// std::vector<int>::iterator it;

	//auto it = std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(un));
	//std::vector<int> un_old;
	//un_old = un;
	//un.clear();
	//for (auto el : un) {
	//	std::cout << ' ' << el;
	//}
	//std::cout << '\n';

	//for (auto el : un_old) {
	//	std::cout << ' ' << el;
	//}
	//std::cout << '\n';


	//std::set_union(un_old.begin(), un_old.end(), set3.begin(), set3.end(), std::back_inserter(un));


	//for (auto el : un) {
	//	std::cout << ' ' << el;
	//}
	//std::cout << '\n';

	//std::cout << n_choose_k(4, 2) << std::endl;
}
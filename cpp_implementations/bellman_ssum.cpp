#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>
#include <chrono>
typedef unsigned long long int ull;

std::pair<bool, double> bellman_ssum(const std::vector<int>& set, int target) {
	std::vector<std::vector<int>> table(set.size(), std::vector<int>(target + 1, 0));


	auto start_time = std::chrono::high_resolution_clock::now(); // time measurement

	/* Complete the first row */
	if (set[0] <= target) {
		table[0][set[0]] = 1;
	}

	/* Fill up the first column with ones */
	for (int i = 0; i < set.size(); i++) {
		table[i][0] = 1;
	}

	for (int i = 1; i < set.size(); i++) {
		for (int j = 1; j < target + 1; j++) {
			if (j < set[i]){
				table[i][j] = table[i - 1][j];
			}
			else if (set[i] == j) {
				table[i][j] = 1;
			}
			else {
				// std::cout << j - set[i] << std::endl;
				table[i][j] = std::max(table[i-1][j], table[i-1][j-set[i]]);
			}
		}
	}

	bool solution = table[set.size() - 1][target] == 1;

	auto end_time = std::chrono::high_resolution_clock::now(); // time measurement

	std::chrono::duration<double, std::milli> running_time = end_time - start_time;

	std::pair<bool, double> result(solution, running_time.count());
	
	return result;

	/*  */

}
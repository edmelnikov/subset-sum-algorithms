/*
Implementation of dynamic algorithm
*/
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>
#include <chrono>

/* Calculate a dynamic table and return it */
std::vector<std::vector<int>> calc_dyn_table(const std::vector<int>& set, int target) {
	/* Handle trivial case*/
	if (set.size() == 0) {
		return std::vector<std::vector<int>>(1, std::vector<int>({ 0 }));
	}

	/* Create table filled with zeros */
	std::vector<std::vector<int>> table(set.size(), std::vector<int>(target + 1, 0));
	
	/* Complete the first row */
	if (set[0] <= target) {
		table[0][set[0]] = 1;
	}

	/* Fill up the first column with ones */
	for (int i = 0; i < set.size(); i++) {
		table[i][0] = 1;
	}

	/* Iterate over the rest of the cells of the table*/
	for (int i = 1; i < set.size(); i++) {
		for (int j = 1; j < target + 1; j++) {
			if (j < set[i]) {
				table[i][j] = table[i - 1][j];
			}
			else if (set[i] == j) {
				table[i][j] = 1;
			}
			else {
				table[i][j] = std::max(table[i - 1][j], table[i - 1][j - set[i]]);
			}
		}
	}
	return table;
}

/* Get all subset sums from the dynamic table */
std::vector<int> get_ssums_dyn_table(const std::vector<int>& set, int target) {
	std::vector<std::vector<int>> tab = std::move(calc_dyn_table(set, target));
	std::vector<int> last_row =	std::move(tab.back());
	std::vector<int> ssums;
	for (int i = 0; i < last_row.size(); i++) {
		if (last_row[i] == 1) {
			ssums.push_back(i);
		}
	}
	return ssums;
}

/* Solve the subset sum problem with dynamic table and measure time */
std::pair<bool, double> bellman_ssum(const std::vector<int>& set, int target) {

	auto start_time = std::chrono::high_resolution_clock::now(); // time measurement
	std::vector<std::vector<int>> table = std::move(calc_dyn_table(set, target));
	bool solution = table[set.size() - 1][target] == 1;

	auto end_time = std::chrono::high_resolution_clock::now(); // time measurement

	std::chrono::duration<double, std::milli> running_time = end_time - start_time;

	std::pair<bool, double> result(solution, running_time.count());
	
	return result;
}
#pragma once
#include <vector>

std::vector<std::vector<int>> calc_dyn_table(const std::vector<int>& set, int target);
std::vector<int> get_ssums_dyn_table(const std::vector<int>& set, int target);
std::pair<bool, double> bellman_ssum(const std::vector<int>& set, int target);

#pragma once
#include <vector>

std::vector<int> fft_polymul(const std::vector<int>& poly1, const std::vector<int>& poly2);
std::vector<int> minkowski_add(const std::vector<int>& set1, const std::vector<int>& set2, int bound=-1);

std::vector<std::vector<int>> fft_polymul_2d(std::vector<std::vector<int>>& poly1, std::vector<std::vector<int>>& poly2);
std::vector<std::vector<int>> minkowski_add_2d(const std::vector<std::vector<int>>& set1, const std::vector<std::vector<int>>& set2, int bound = -1);
#pragma once
#include <vector>
#include <unordered_map>

std::unordered_map<int, std::vector<int>> rand_part(const std::vector<int>& set, int num_parts);
std::vector<int> ColorCoding(const std::vector<int>& set, int target, int k, double delta);
std::vector<int> ColorCodingLayer(const std::vector<int>& set, int target, int layer, double delta);
std::pair<bool, double> bringmann_ssum(const std::vector<int>& set, int target, double delta);


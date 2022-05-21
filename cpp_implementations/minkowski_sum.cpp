#include "fftw3.h"
#include <iostream>
#include <stdio.h>
#include <complex.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>

#include <stdio.h>
#include <mpir.h>
#include <mpirxx.h>
#include "minkowski_sum.h"
#include "convolution_fftw.h"
#include "factorize.h"


int num_mink_calls = 0;
//
//std::vector<int> fft_polymul(const std::vector<int> &poly1, const std::vector<int>& poly2) {
//	int size1 = poly1.size();
//	int size2 = poly2.size();
//	
//	//int size_mul = std::max(poly1.size(), poly2.size());
//	//int N = pow(2, ceil(log(size_mul*2 + 1) / log(2)));	
//	//int N = pow(2, ceil(log(size1 + size2 + 1) / log(2)));// size of an array with resulting polynomial coefficients rounded up to the next power of two
//	int N = size1 + size2 - 1;
//	
//	// and the number of roots of unity used for evaluation (-)
//		
//	fftw_plan plan;
//	fftw_complex* fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);  // array where input polynomial padded by zeros will be stored
//	fftw_complex* fft_output1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // output of the FFT for the first polynomial
//	fftw_complex* fft_output2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // output of the FFT for the second polynomial
//	
//	/* First polynomial evaluation  */
//	/* Fill up the array used as an fft input and pad it with zeros up to the length of N */
//	for (int i = 0; i < N; i++) {
//		if (i < size1) {
//			fft_input[i][0] = poly1[i];
//		}
//		else {
//			fft_input[i][0] = 0;
//		}
//		fft_input[i][1] = 0;
//	}
//
//	plan = fftw_plan_dft_1d(N, fft_input, fft_output1, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//
//
//	/* Second polynomial evaluation */
//	for (int i = 0; i < N; i++) {
//		if (i < size2) {
//			fft_input[i][0] = poly2[i];
//		}
//		else {
//			fft_input[i][0] = 0;
//		}
//		fft_input[i][1] = 0;
//	}
//	plan = fftw_plan_dft_1d(N, fft_input, fft_output2, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//
//	/* Polynomial evaluations multiplication */
//	fftw_complex* poly_mul_evaluations = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N));
//	for (int i = 0; i < N; i++) {
//		poly_mul_evaluations[i][0] = fft_output1[i][0] * fft_output2[i][0] - fft_output1[i][1] * fft_output2[i][1]; // real part
//		poly_mul_evaluations[i][1] = fft_output1[i][0] * fft_output2[i][1] + fft_output1[i][1] * fft_output2[i][0]; // imaginary part 
//		//std::cout << poly_mul_evaluations[i][0] << " + " << poly_mul_evaluations[i][1] << "i" << std::endl;
//	}
//
//	/* Interpolation */
//	fftw_complex* poly_mul = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N));  // multiplication product of two polynomials
//	// std::cout << poly_mul[10][0] << std::endl;
//	
//	plan = fftw_plan_dft_1d(N, poly_mul_evaluations, poly_mul, FFTW_BACKWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//
//	fftw_destroy_plan(plan);
//
//	/* Divide all the values by N since FFTW computes unnormalized DFT  */
//	std::vector<int> res(N);
//	for (int i = 0; i < N; i++) {
//		res[i] = abs(round(poly_mul[i][0] / N));
//		// std::cout << res[i] << ' ';
//	}
//	//std::cout << std::endl;
//	return res;
//}
//
//std::vector<int> minkowski_add(const std::vector<int>& set1, const std::vector<int>& set2, int bound) {
//	
//	/* initialize characteristic polynomials of sets with zeros */
//	// here, the size of a vector is "the largest integer in  a set (last item of a SORTED set)" + 1
//	// could be replaced with bound + 1
//
//	std::vector<int> set1_poly;
//	std::vector<int> set2_poly;
//
//	/* fill up both characteristic polynomials */
//	if (bound == -1) { // when there is no bound
//		set1_poly.resize(
//			(set1.size() > 0) ? (set1.back() + 1) : (1),
//			0
//		);
//		set2_poly.resize(
//			(set2.size() > 0) ? (set2.back() + 1) : (1),
//			0
//		);
//
//		for (auto element : set1) {
//			set1_poly[element] += 1;
//		}
//		set1_poly[0] = 1; // monomial of degree 0 is added by default
//		for (auto element : set2) {
//			set2_poly[element] += 1;
//		}
//		set2_poly[0] = 1; // monomial of degree 0 is added by default
//
//	}
//	else { // otherwise, if there is a specified bound 
//		set1_poly.resize(
//			(set1.size() > 0) ? (std::min(bound, set1.back()) + 1) : (1),
//			0
//		);
//		set2_poly.resize(
//			(set2.size() > 0) ? (std::min(bound, set2.back()) + 1) : (1),
//			0
//		);
//
//		for (auto element : set1) {
//			if (element <= set1_poly.size()) {
//				set1_poly[element] += 1;
//			}
//			else break;
//		}
//		set1_poly[0] = 1; // monomial of degree 0 is added by default
//		for (auto element : set2) {
//			if (element <= set2_poly.size()) {
//				set2_poly[element] += 1;
//			}
//			else break;
//		}
//		set2_poly[0] = 1; // monomial of degree 0 is added by default
//	}
//		
//	//for (auto element : set1_poly) {
//	//	std::cout << element << ' ';
//	//}
//	//std::cout << std::endl;
//	//for (auto element : set2_poly) {
//	//	std::cout << element << ' ';
//	//}
//	//std::cout << std::endl;
//
//	std::vector<int> product_poly = fft_polymul(set1_poly, set2_poly); // characteristic polynomial of minkowski addition product
//
//	//for (auto element : product_poly) {
//	//	std::cout << element << ' ';
//	//}
//	//std::cout << std::endl;
//
//	std::vector<int> product_set; // minkowski addition product (set)
//	product_set.reserve(product_poly.size()); // reserve enough space so that the vector is now resized during .push_back()'s
//
//	for (int i = 1; i < product_poly.size(); i++) {
//		if (product_poly[i] != 0 && ((bound == -1) || (bound != 1 && i <= bound))) { // if bound is not set
//			product_set.push_back(i);
//
//		}
//		//if (product_poly[i] != 0) { // if bound is not set
//		//	product_set.push_back(i);
//		//}
//	}
//	return product_set;
//
//}
//
//std::vector<std::vector<int>> fft_polymul_2d(std::vector<std::vector<int>>& poly1, std::vector<std::vector<int>>& poly2) {
//	// int max_size1 = std::max(poly1.size(), poly1[0].size());
//	// int max_size2 = std::max(poly2.size(), poly2[0].size());
//
//	//int num_rows = pow(2, ceil(log(poly1.size() + poly2.size() - 1) / log(2))); // num of rows of a resulting matrix with resulting polynomial
//	int num_rows = poly1.size() + poly2.size() - 1;
//	// coefficients rounded up to the next power of two
//	//int num_cols = pow(2, ceil(log(poly1.front().size() + poly2.front().size() - 1) / log(2))); // num of cols of a resulting matrix with resulting polynomial
//	int num_cols = poly1.front().size() + poly2.front().size() - 1;
//	// coefficients rounded up to the next power of two
//
//	fftw_plan plan;
//	fftw_complex* fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows * num_cols);  // array where input polynomial padded by zeros will be stored
//	fftw_complex* fft_output1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows * num_cols); // output of the FFT for the first polynomial
//	fftw_complex* fft_output2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows * num_cols); // output of the FFT for the second polynomial
//
//	/* first polynomial evaluation  */
//	/* fill up the array used as an fft input and pad it with zeros up to the length of n */
//	for (int i = 0; i < num_rows; i++) {
//		for (int j = 0; j < num_cols; j++) {
//			if (i < poly1.size() && j < poly1[i].size()) {
//				fft_input[i * num_cols + j][0] = poly1[i][j];
//			}
//			else {
//				fft_input[i * num_cols + j][0] = 0;
//			}
//			fft_input[i * num_cols + j][1] = 0;
//		}
//	}
//
//	//for (int i = 0; i < N; i++) {
//	//	for (int j = 0; j < N; j++) {
//	//		std::cout<< fft_input[N*i + j][0] << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//std::cout << std::endl;
//
//	plan = fftw_plan_dft_2d(num_rows, num_cols, fft_input, fft_output1, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//
//	//for (int i = 0; i < N; i++) {
//	//	for (int j = 0; j < N; j++) {
//	//		std::cout << fft_output1[N*i + j][0] << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//std::cout << std::endl;
//
//	///* Second polynomial evaluation */
//	for (int i = 0; i < num_rows; i++) {
//		for (int j = 0; j < num_cols; j++) {
//			if (i < poly2.size() && j < poly2[i].size()) {
//				fft_input[i * num_cols + j][0] = poly2[i][j];
//
//			}
//			else {
//				fft_input[i * num_cols + j][0] = 0;
//			}
//			fft_input[i * num_cols + j][1] = 0;
//		}
//	}
//
//	//for (int i = 0; i < N; i++) {
//	//	for (int j = 0; j < N; j++) {
//	//		std::cout << fft_input[N * i + j][0] << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//std::cout << std::endl;
//
//	plan = fftw_plan_dft_2d(num_rows, num_cols, fft_input, fft_output2, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//
//	//for (int i = 0; i < N; i++) {
//	//	for (int j = 0; j < N; j++) {
//	//		std::cout << fft_output2[N * i + j][0] << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//std::cout << std::endl;
//
//	/* Polynomial evaluations multiplication */
//	fftw_complex* poly_mul_evaluations = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows * num_cols);
//	 
//	for (int i = 0; i < num_rows; i++) {
//		for (int j = 0; j < num_cols; j++) {
//			poly_mul_evaluations[num_cols * i + j][0] = fft_output1[num_cols * i + j][0] * fft_output2[num_cols * i + j][0] - fft_output1[num_cols * i + j][1] * fft_output2[num_cols * i + j][1]; // real part
//			poly_mul_evaluations[num_cols * i + j][1] = fft_output1[num_cols * i + j][0] * fft_output2[num_cols * i + j][1] + fft_output1[num_cols * i + j][1] * fft_output2[num_cols * i + j][0]; // imaginary part 
//		}
//	}
//
//	/* Interpolation */
//	fftw_complex* poly_mul = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows * num_cols);  // multiplication product of two polynomials
//	plan = fftw_plan_dft_2d(num_rows, num_cols, poly_mul_evaluations, poly_mul, FFTW_BACKWARD, FFTW_ESTIMATE);
//	fftw_execute(plan);
//	fftw_destroy_plan(plan);
//
//	/* Divide all the values since FFTW computes unnormalized DFT  */
//	std::vector<std::vector<int>> res(num_rows, std::vector<int>(num_cols, 0));
//	for (int i = 0; i < num_rows; i++) {
//		for (int j = 0; j < num_cols; j++) {
//			res[i][j] = abs(round(poly_mul[num_cols * i + j][0] /(num_rows * num_cols)));
//		}
//	}
//	return res;
//
//}
//
//std::vector<std::pair<int, int>> minkowski_add_2d(const std::vector<std::pair<int, int>>& set1, const std::vector<std::pair<int, int>>& set2, int bound) {
//
//	/* initialize characteristic polynomials of sets with zeros */
//	// here, the size of a vector is "the largest integer in  a set (last item of a SORTED set)" + 1
//	// could be replaced with bound + 1
//	std::vector<std::vector<int>> set1_poly;
//	std::vector<std::vector<int>> set2_poly;
//
//	/* fill up both characteristic polynomials */
//	if (bound == -1) { // when there is no bound
//		/* set1 */
//		int max_first = 0;
//		int max_second = 0;
//		for (auto el : set1) {
//			if (el.first > max_first) {
//				max_first = el.first;
//			}
//			if (el.second > max_second) {
//				max_second = el.second;
//			}
//		}
//		set1_poly.resize(max_first + 1, std::vector<int>(max_second + 1, 0));
//
//		for (auto el : set1) {
//			set1_poly[el.first][el.second] += 1;
//		}
//		set1_poly[0][0] += 1;
//
//		//for (auto row : set1_poly) {
//		//	for (auto el : row) {
//		//		std::cout << el << " ";
//		//	}
//		//	std::cout << std::endl;
//		//}
//		//std::cout << std::endl;
//
//		/* set2 */
//		max_first = 0;
//		max_second = 0;
//		for (auto el : set2) {
//			if (el.first > max_first) {
//				max_first = el.first;
//			}
//			if (el.second > max_second) {
//				max_second = el.second;
//			}
//		}
//		set2_poly.resize(max_first + 1, std::vector<int>(max_second + 1, 0));
//
//		for (auto el : set2) {
//			set2_poly[el.first][el.second] += 1;
//		}
//		set2_poly[0][0] += 1;
//
//		//for (auto row : set2_poly) {
//		//	for (auto el : row) {
//		//		std::cout << el << " ";
//		//	}
//		//	std::cout << std::endl;
//		//}
//		//std::cout << std::endl;
//
//
//	}
//	else { // if there is a bound 
//		/* set1 */
//		int max_first = 0;
//		int max_second = 0;
//		for (auto el : set1) {
//			if (el.first > max_first) {
//				max_first = el.first;
//			}
//			if (el.second > max_second) {
//				max_second = el.second;
//			}
//		}
//		set1_poly.resize(
//			((max_first >= bound) ? (bound) : (max_first)) + 1,
//			std::vector<int>(max_second + 1, 0)
//		);
//		for (auto el : set1) {
//			if (el.first <= bound) {
//				set1_poly[el.first][el.second] += 1;
//			}
//		}
//		set1_poly[0][0] += 1;
//
//		//for (auto row : set1_poly) {
//		//	for (auto el : row) {
//		//		std::cout << el << " ";
//		//	}
//		//	std::cout << std::endl;
//		//}
//		//std::cout << std::endl;
//
//		/* set2 */
//		max_first = 0;
//		max_second = 0;
//		for (auto el : set2) {
//			if (el.first > max_first) {
//				max_first = el.first;
//			}
//			if (el.second > max_second) {
//				max_second = el.second;
//			}
//		}
//		set2_poly.resize(
//			((max_first >= bound) ? (bound) : (max_first)) + 1,
//			std::vector<int>(max_second + 1, 0)
//		);
//		for (auto el : set2) {
//			if (el.first <= bound) {
//				set2_poly[el.first][el.second] += 1;
//
//			}
//		}
//		set2_poly[0][0] += 1;
//	}
//
//	std::vector<std::vector<int>> product_poly = fft_polymul_2d(set1_poly, set2_poly); // characteristic polynomial of minkowski addition product
//
//	std::vector<std::pair<int, int>> product_set; // minkowski addition product (set)
//	// product_set.reserve(product_poly.size() * product_poly.size()); // reserve enough space so that the vector is not resized during .push_back()'s (lazy solution)
//
//	for (int i = 0; i < product_poly.size(); i++) {
//		for (int j = 0; j < product_poly[i].size(); j++) {
//			if (product_poly[i][j] != 0) {
//				// std::vector<int> pair = { i, j };
//				if (bound == -1) product_set.push_back(std::pair<int, int>({i, j})); // NOTE: (0, 0) is also appended
//				else if (bound != 1 && i <= bound) product_set.push_back(std::pair<int, int>({ i, j }));
//			}
//		}
//	}
//
//	return product_set;
//
//}

using namespace FFTW_Convolution;

std::vector<int> minkowski_sum_fft(const std::vector<int>& set1, const std::vector<int>& set2, int bound) {

	/* initialize characteristic polynomials of sets with zeros */
	// here, the size of a vector is "the largest integer in a set (last item of a SORTED set)" + 1
	// could be replaced with bound + 1

	int max;
	/* set1 */
	max = 0;
	for (auto el : set1) {
		if (el > max) {
			max = el;
		}
	}
	int set1_poly_size = ((max >= bound) ? (bound) : (max)) + 1;

	/* set2 */
	max = 0;
	for (auto el : set2) {
		if (el > max) {
			max = el;
		}
	}
	int set2_poly_size = ((max >= bound) ? (bound) : (max)) + 1;;

	double* set1_poly = new double[set1_poly_size]();
	double* set2_poly = new double[set2_poly_size]();

	for (auto element : set1) {
		if (element <= set1_poly_size - 1) {
			set1_poly[element] += 1;
		}
		else break;
	}
	set1_poly[0] = 1; // monomial of degree 0 is added by default
	for (auto element : set2) {
		if (element <= set2_poly_size - 1) {
			set2_poly[element] += 1;
		}
		else break;
	}
	set2_poly[0] = 1; // monomial of degree 0 is added by default

	/* Perform convolution */
	Workspace ws;
	init_workspace(ws, LINEAR_FULL, 1, set1_poly_size, 1, set2_poly_size);
	convolve(ws, set1_poly, set2_poly);
	
	/* Clear memory */
	delete[] set1_poly;
	delete[] set2_poly;

	std::vector<int> product_set; // minkowski addition product (set)
	//product_set.reserve(ws.w_dst - 1); // reserve enough space so that the vector is now resized during .push_back()'s

	for (int i = 1; i < ws.w_dst; i++) {
		if (abs(ws.dst[i]) > 0.001 && i <= bound) {
			product_set.push_back(i);
		}
	}
	clear_workspace(ws);
	return product_set;
}

std::vector<std::pair<int, int>> minkowski_sum_fft2d(const std::vector<std::pair<int, int>>& set1, const std::vector<std::pair<int, int>>& set2, int bound) {
	
	int max_first;
	int max_second;
	
	/* set1 */
	max_first = 0;
	max_second = 0;
	for (auto el : set1) {
		if (el.first > max_first) {
			max_first = el.first;
		}
		if (el.second > max_second) {
			max_second = el.second;
		}
	}
	int set1_poly_h = ((max_first >= bound) ? (bound) : (max_first)) + 1;
	int set1_poly_w = ((max_second >= bound) ? (bound) : (max_second)) + 1;

	/* set2 */
	max_first = 0;
	max_second = 0;
	for (auto el : set2) {
		if (el.first > max_first) {
			max_first = el.first;
		}
		if (el.second > max_second) {
			max_second = el.second;
		}
	}
	int set2_poly_h = ((max_first >= bound) ? (bound) : (max_first)) + 1;
	int set2_poly_w = ((max_second >= bound) ? (bound) : (max_second)) + 1;

	double* set1_poly = new double[set1_poly_w*set1_poly_h]();
	double* set2_poly = new double[set2_poly_w*set2_poly_h]();

	for (auto pair : set1) {
		if (pair.first < set1_poly_h && pair.second < set1_poly_w) {
			set1_poly[pair.first*set1_poly_w + pair.second] += 1;
		}
	}
	set1_poly[0] = 1;

	//for (int i = 0; i < set1_poly_h; i++) {
	//	for (int j = 0; j < set1_poly_w; j++) {
	//		std::cout << set1_poly[i* set1_poly_w + j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	for (auto pair : set2) {
		if (pair.first < set2_poly_h && pair.second < set2_poly_w) {
			set2_poly[pair.first * set2_poly_w + pair.second] += 1;
		}
	}
	set2_poly[0] = 1;

	//for (int i = 0; i < set2_poly_h; i++) {
	//	for (int j = 0; j < set2_poly_w; j++) {
	//		std::cout << set2_poly[i * set2_poly_w + j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	/* Perform convolution */
	Workspace ws;
	init_workspace(ws, LINEAR_FULL, set1_poly_h, set1_poly_w, set2_poly_h, set2_poly_w);
	convolve(ws, set1_poly, set2_poly);

	/* Clear memory */
	delete[] set1_poly;
	delete[] set2_poly;

	//for (int i = 0; i < ws.h_dst; i++) {
	//	for (int j = 0; j < ws.w_dst; j++) {
	//		if (abs(ws.dst[i * ws.w_dst + j]) > 0.001) std::cout << ws.dst[i * ws.w_dst + j] << " ";
	//		else std::cout << 0 << " ";
	//	}
	//	std::cout << std::endl;
	//}

	std::vector<std::pair<int, int>> product_set;

	for (int i = 0; i < ws.h_dst; i++) {
		for (int j = 0; j < ws.w_dst; j++) {

			if (abs(ws.dst[i * ws.w_dst + j]) > 0.001 && i <= bound) {
				product_set.push_back({i, j});
			}
		}
	}
	clear_workspace(ws);

	return product_set;
}


std::vector<int> minkowski_sum_mpir(const std::vector<int>& set1, const std::vector<int>& set2, int bound) {
	/* Increase the counter */
	num_mink_calls++;
	
	/* Handle trivial cases */
	if (set1.size() == 0 && set2.size() == 0) return set1;
	if (set1.size() == 0) return set2;
	if (set2.size() == 0) return set1;

	int padding_size = ceil(log2(bound));
	int block_size = padding_size + 1;

	// int bool_vec_size1 = ((set1.back() < bound) ? set1.back() : bound) + 1; // sizes of characteristic vectors of sets
	// int bool_vec_size2 = ((set2.back() < bound) ? set2.back() : bound) + 1;
	int bool_vec_size1 = bound + 1; // sizes of characteristic vectors of sets
	int bool_vec_size2 = bound + 1;
	int bool_vec_pad_size1 = bool_vec_size1 * block_size; // sizes of characteristic vectors with paddings
	int bool_vec_pad_size2 = bool_vec_size2 * block_size;
	
	// std::cout << bool_vec_pad_size1 << std::endl;
	
	/* Create binary vectors of zeros with paddings and populate them with ones*/
	std::string bool_vec1(bool_vec_pad_size1, '0');
	std::string bool_vec2(bool_vec_pad_size2, '0');

	bool_vec1[0] = '1';
	bool_vec2[0] = '1';

	for (auto el : set1) {
		if (el > bound) break;
		else {
			bool_vec1[el * block_size] = '1';
		}
	}
	for (auto el : set2) {
		if (el > bound) break;
		else {
			bool_vec2[el * block_size] = '1';
		}
	}
	
	// std::cout << bool_vec1 << std::endl;
	//  std::cout << bool_vec2 << std::endl;

	/* Multiply two binary numbers */
	mpz_class bin1(bool_vec1, 2);
	mpz_class bin2(bool_vec2, 2);
	std::string prod_bool_vec = mpz_class(bin1 * bin2).get_str(2);
	
	/* Get a sumset from the characteristic vector */
	int indent = prod_bool_vec.length() % block_size; // the amount of bytes for the first block representing 0 in a set
	if (indent == 0) indent = block_size;


	int num_blocks = prod_bool_vec.length() / block_size; // the amount of whole blocks 


	int block_start_ind = 0; 
	int block_end_ind = 0;

	mpz_t block_bin;
	mpz_init(block_bin);
	mpz_t zero;
	mpz_init(zero);

	std::vector<int> mink_sum_set;

	for (int block_ind = 0; block_ind < num_blocks && block_ind < bound; block_ind++) { // iterate over whole blocks 
		block_start_ind = indent + block_ind * block_size;
		
		// const char* block_char = (prod_bool_vec.substr(block_start_ind, block_size)).c_str();

		mpz_set_str(
			block_bin,
			(prod_bool_vec.substr(block_start_ind, block_size)).c_str(),
			2);
		
		if (mpz_cmp(block_bin, zero) > 0) {
			mink_sum_set.push_back(block_ind + 1);
		}

		//std::cout << mpz_cmp(block_bin, zero) << std::endl;

		//  std::cout << prod_bool_vec.substr(block_start_ind, block_size) << " ";

	}
	// std::cout << std::endl;


	// std::cout << prod_bool_vec.length() << std::endl;
	// std::cout << prod_bool_vec << std::endl;
	return mink_sum_set;
}

/* Power of two using bit shift */
int pow2_bs(int pow) { 
	return 1 << pow;
}
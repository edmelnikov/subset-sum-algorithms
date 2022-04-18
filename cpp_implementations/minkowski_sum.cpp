#include "fftw3.h"
#include <iostream>
#include <stdio.h>
#include <complex.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>

std::vector<int> fft_polymul(const std::vector<int> &poly1, const std::vector<int>& poly2) {
	int size1 = poly1.size();
	int size2 = poly2.size();
	
	//int size_mul = std::max(poly1.size(), poly2.size());
	//int N = pow(2, ceil(log(size_mul*2 + 1) / log(2)));  
	int N = pow(2, ceil(log(size1 + size2 + 1) / log(2)));// size of an array with resulting polynomial coefficients rounded up to the next power of two
	// and the number of roots of unity used for evaluation (-)
		
	fftw_plan plan;
	fftw_complex* fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);  // array where input polynomial padded by zeros will be stored
	fftw_complex* fft_output1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // output of the FFT for the first polynomial
	fftw_complex* fft_output2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); // output of the FFT for the second polynomial
	
	/* First polynomial evaluation  */
	/* Fill up the array used as an fft input and pad it with zeros up to the length of N */
	for (int i = 0; i < N; i++) {
		if (i < size1) {
			fft_input[i][0] = poly1[i];
		}
		else {
			fft_input[i][0] = 0;
		}
		fft_input[i][1] = 0;
	}	
	plan = fftw_plan_dft_1d(N, fft_input, fft_output1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	/* Second polynomial evaluation */
	for (int i = 0; i < N; i++) {
		if (i < size2) {
			fft_input[i][0] = poly2[i];
		}
		else {
			fft_input[i][0] = 0;
		}
		fft_input[i][1] = 0;
	}
	plan = fftw_plan_dft_1d(N, fft_input, fft_output2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	/* Polynomial evaluations multiplication */
	fftw_complex* poly_mul_evaluations = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N));
	for (int i = 0; i < N; i++) {
		poly_mul_evaluations[i][0] = fft_output1[i][0] * fft_output2[i][0] - fft_output1[i][1] * fft_output2[i][1]; // real part
		poly_mul_evaluations[i][1] = fft_output1[i][0] * fft_output2[i][1] + fft_output1[i][1] * fft_output2[i][0]; // imaginary part 
		//std::cout << poly_mul_evaluations[i][0] << " + " << poly_mul_evaluations[i][1] << "i" << std::endl;
	}

	/* Interpolation */
	fftw_complex* poly_mul = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N));  // multiplication product of two polynomials
	// std::cout << poly_mul[10][0] << std::endl;
	
	plan = fftw_plan_dft_1d(N, poly_mul_evaluations, poly_mul, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	fftw_destroy_plan(plan);

	/* Divide all the values by N since FFTW computes unnormalized DFT  */
	std::vector<int> res(N);
	for (int i = 0; i < N; i++) {
		res[i] = abs(round(poly_mul[i][0] / N));
		// std::cout << res[i] << ' ';
	}
	//std::cout << std::endl;
	return res;
}

std::vector<int> minkowski_add(const std::vector<int>& set1, const std::vector<int>& set2, int bound=-1) {
	
	/* initialize characteristic polynomials of sets with zeros */
	// here, the size of a vector is "the largest integer in  a set (last item of a SORTED set)" + 1
	// could be replaced with bound + 1

	std::vector<int> set1_poly;
	std::vector<int> set2_poly;

	/* fill up both characteristic polynomials */
	if (bound == -1) { // when there is no bound
		set1_poly.resize(
			(set1.size() > 0) ? (set1.back() + 1) : (1),
			0
		);
		set2_poly.resize(
			(set2.size() > 0) ? (set2.back() + 1) : (1),
			0
		);

		for (auto element : set1) {
			set1_poly[element] += 1;
		}
		set1_poly[0] = 1; // monomial of degree 0 is added by default
		for (auto element : set2) {
			set2_poly[element] += 1;
		}
		set2_poly[0] = 1; // monomial of degree 0 is added by default

	}
	else { // otherwise, if there is a specified bound 
		set1_poly.resize(
			(set1.size() > 0) ? (std::min(bound, set1.back()) + 1) : (1),
			0
		);
		set2_poly.resize(
			(set2.size() > 0) ? (std::min(bound, set2.back()) + 1) : (1),
			0
		);

		for (auto element : set1) {
			if (element <= set1_poly.size()) {
				set1_poly[element] += 1;
			}
			else break;
		}
		set1_poly[0] = 1; // monomial of degree 0 is added by default
		for (auto element : set2) {
			if (element <= set2_poly.size()) {
				set2_poly[element] += 1;
			}
			else break;
		}
		set2_poly[0] = 1; // monomial of degree 0 is added by default
	}
		
	//for (auto element : set1_poly) {
	//	std::cout << element << ' ';
	//}
	//std::cout << std::endl;
	//for (auto element : set2_poly) {
	//	std::cout << element << ' ';
	//}
	//std::cout << std::endl;

	std::vector<int> product_poly = fft_polymul(set1_poly, set2_poly); // characteristic polynomial of minkowski addition product

	//for (auto element : product_poly) {
	//	std::cout << element << ' ';
	//}
	//std::cout << std::endl;

	std::vector<int> product_set; // minkowski addition product (set)
	product_set.reserve(product_poly.size()); // reserve enough space so that the vector is now resized during .push_back()'s

	for (int i = 1; i < product_poly.size(); i++) {
		if (product_poly[i] != 0 && ((bound == -1) || (bound != 1 && i <= bound))) { // if bound is not set
			product_set.push_back(i);

		}
		//if (product_poly[i] != 0) { // if bound is not set
		//	product_set.push_back(i);
		//}
	}
	return product_set;

}

std::vector<std::vector<int>> fft_polymul_2d(std::vector<std::vector<int>>& poly1, std::vector<std::vector<int>>& poly2) {
	int max_size1 = std::max(poly1.size(), poly1[0].size());
	int max_size2 = std::max(poly2.size(), poly2[0].size());
	int N = pow(2, ceil(log(2*std::max(max_size1, max_size2) + 1) / log(2))); // dimension of a resulting matrix with resulting polynomial coefficients rounded up to the next power of two
	// std::cout << ceil(log(2 * std::max(max_size1, max_size2) + 1) / log(2)) << std::endl;

	fftw_plan plan;
	fftw_complex* fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N*N);  // array where input polynomial padded by zeros will be stored
	fftw_complex* fft_output1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N*N); // output of the FFT for the first polynomial
	fftw_complex* fft_output2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N*N); // output of the FFT for the second polynomial

	/* first polynomial evaluation  */
	/* fill up the array used as an fft input and pad it with zeros up to the length of n */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i < poly1.size() && j < poly1[i].size()) {
				fft_input[i * N + j][0] = poly1[i][j];
			
			}
			else {
				fft_input[i * N + j][0] = 0;
			}
			fft_input[i * N + j][1] = 0;
		}
	}

	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		std::cout<< fft_input[N*i + j][0] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	plan = fftw_plan_dft_2d(N, N, fft_input, fft_output1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		std::cout << fft_output1[N*i + j][0] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	///* Second polynomial evaluation */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i < poly2.size() && j < poly2[i].size()) {
				fft_input[i * N + j][0] = poly2[i][j];

			}
			else {
				fft_input[i * N + j][0] = 0;
			}
			fft_input[i * N + j][1] = 0;
		}
	}

	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		std::cout << fft_input[N * i + j][0] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	plan = fftw_plan_dft_2d(N, N, fft_input, fft_output2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		std::cout << fft_output2[N * i + j][0] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	/* Polynomial evaluations multiplication */
	fftw_complex* poly_mul_evaluations = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			poly_mul_evaluations[N * i + j][0] = fft_output1[N * i + j][0] * fft_output2[N * i + j][0] - fft_output1[N * i + j][1] * fft_output2[N * i + j][1]; // real part
			poly_mul_evaluations[N * i + j][1] = fft_output1[N * i + j][0] * fft_output2[N * i + j][1] + fft_output1[N * i + j][1] * fft_output2[N * i + j][0]; // imaginary part 
		}
	}

	/* Interpolation */
	fftw_complex* poly_mul = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);  // multiplication product of two polynomials
	plan = fftw_plan_dft_2d(N, N, poly_mul_evaluations, poly_mul, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	/* Divide all the values by N since FFTW computes unnormalized DFT  */
	std::vector<std::vector<int>> res(N, std::vector<int>(N, 0));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			res[i][j] = abs(round(poly_mul[N*i+j][0] /(N*N)));
		}
	}
	return res;

}

/* TODO: add boundaries */
std::vector<std::pair<int, int>> minkowski_add_2d(const std::vector<std::pair<int, int>>& set1, const std::vector<std::pair<int, int>>& set2, int bound = -1) {

	/* initialize characteristic polynomials of sets with zeros */
	// here, the size of a vector is "the largest integer in  a set (last item of a SORTED set)" + 1
	// could be replaced with bound + 1

	std::vector<std::vector<int>> set1_poly;
	std::vector<std::vector<int>> set2_poly;

	/* fill up both characteristic polynomials */
	if (bound == -1) { // when there is no bound
		/* set1 */
		int max_first = 0;
		int max_second = 0;
		for (auto el : set1) {
			if (el.first > max_first) {
				max_first = el.first;
			}
			if (el.second > max_second) {
				max_second = el.second;
			}
		}
		set1_poly.resize(max_first + 1, std::vector<int>(max_second + 1, 0));

		for (auto el : set1) {
			set1_poly[el.first][el.second] += 1;
		}
		set1_poly[0][0] += 1;

		//for (auto row : set1_poly) {
		//	for (auto el : row) {
		//		std::cout << el << " ";
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << std::endl;

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
		set2_poly.resize(max_first + 1, std::vector<int>(max_second + 1, 0));

		for (auto el : set2) {
			set2_poly[el.first][el.second] += 1;
		}
		set2_poly[0][0] += 1;

		//for (auto row : set2_poly) {
		//	for (auto el : row) {
		//		std::cout << el << " ";
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << std::endl;

	}

	std::vector<std::vector<int>> product_poly = fft_polymul_2d(set1_poly, set2_poly); // characteristic polynomial of minkowski addition product

	std::vector<std::pair<int, int>> product_set; // minkowski addition product (set)
	product_set.reserve(product_poly.size() * product_poly.size()); // reserve enough space so that the vector is now resized during .push_back()'s (lazy solution)

	for (int i = 0; i < product_poly.size(); i++) {
		for (int j = 0; j < product_poly[i].size(); j++) {
			if (product_poly[i][j] != 0) {
				// std::vector<int> pair = { i, j };
				product_set.push_back(std::pair<int, int>({i, j})); // NOTE: (0, 0) is also appended
			}
		}
	}

	return product_set;

}
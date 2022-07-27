# subset-sum-algorithms

This repository contains implementations of 3 algorithms solving the decision Subset sum problem (https://en.wikipedia.org/wiki/Subset_sum_problem) for my Bachelor's thesis called "FPT-algorithms for the Coin Problem and Knapsack Problem"

## Problem statement
In the Subset sum problem given a set $W$ of $n$ positive unique integers and a target value $t$, one needs to decide whether there is such subset of integers that the sum of its elements is exactly $t$.

## Implemented algorithms
The implemented algorithms are:
- Dynamic $O(nt)$ algorithm (https://en.wikipedia.org/wiki/Subset_sum_problem#Pseudo-polynomial_time_dynamic_programming_solutions)
- Koiliaris & Xu's $\tilde{O}(\sqrt{n}t)$ algorithm (https://doi.org/10.1145/3329863)
- Bringmann's $\tilde{O}(n + t)$ randomised algorithm (https://doi.org/10.48550/arXiv.1610.04712)

The implementation of FFT convolution has been retrieved from https://github.com/jeremyfix/FFTConvolution, all credit to the author. 

## Benchmarks 

All implementations were benchmarked on randomly generated data. The data is a set of different Subset sum problem instances with variable target values ($t$) and cardinalities ($n$). 
The graphs showing the comparison of the execution times of different algorithms on various data are provided below. 

### Time vs. t

![time vs. t](https://github.com/edmelnikov/subset-sum-algorithms/blob/master/img/bench_time_vs_t.png)

### Time vs. n

![time vs. n](https://github.com/edmelnikov/subset-sum-algorithms/blob/master/img/bench_time_vs_n.png)

### Time vs. nt 

![time vs. nt](https://github.com/edmelnikov/subset-sum-algorithms/blob/master/img/bench_time_vs_nt.png)
Note: the dynamic algorithm fails to solve instances with $n$ greater than 25000 due to lack of memory

# Cycle_Patterns
Python program that performs some basic operations related to cycle patterns.
This project is meant as a tool to help check statements related to cycle patterns,
and to visualize them for small graphs.

The code is not very optimized, but should be able to handle graphs with up to a couple thousand cycles.

This project uses a custom graph class made by Paul Bonsma, Pieter Bos and Tariq Bontekoe.

#Prerequisites
Standard Python packages: random, scipy, math, numpy

Other: graphviz (https://graphviz.org/download/)

#Usage

This project has the following main functionalities, use help(_function name_) for more info
- check_pattern  
  - Given a graph, a list of its cycles and a list of assigned signs, check if the cycle pattern is valid
  - If the cycle pattern is not valid, it will print a witness
  - Optionally, if the cycle pattern is valid, it can construct the partially ordered set of cycles, and find the minimal + cycles and maximal - cycles
- check_parity_valid
  - Same as check_pattern, but for parity-realizability
  - Does not have the extra options of check_pattern
  - Not very efficient (after all, problem is coNP-complete)
- make_cycle_list
  - Given a graph, make a list of all its cycles (a cycle is given as a list of vertices)
- find_smallest_witness
  - Finds the smallest witness of graph G<sub>n</sub> with a certain (unrealizable) cycle pattern. This should equal 2<sup>n</sup> for n>=2.

# Reinforcement learning driven heuristic for two-dimensional bandwidth minimization
The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper _Reinforcement learning driven heuristic for two-dimensional bandwidth minimization_ by Q. Zhou, M. Gao, and J.K. Hao.

## Cite
Please cite the paper like this:
```
@article{ZHOU2026114327,
title = {Reinforcement learning driven heuristic for two-dimensional bandwidth minimization},
journal = {Applied Soft Computing},
volume = {187},
pages = {114327},
year = {2026},
issn = {1568-4946},
doi = {https://doi.org/10.1016/j.asoc.2025.114327},
url = {https://www.sciencedirect.com/science/article/pii/S1568494625016400},
author = {Qing Zhou and Ming Gao and Jin-Kao Hao}
}
```

## How to run the programs.
** Instructions to use the source code of the proposed RLTS

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./RLTS_2DBMP ./input_file times seed ./output_file

(where input_file is the instance name, times is the number of independent runs, seed is the given random seed, such as 1, 2, ... times, output_file is the output file name)

q.zhou$

*** To clean

q.zhou$ make clean

q.zhou$

## Materials
This repository includes the following materials:

--Benchmark instances used in our paper (see the directory named instance_2DBMP).

--Source codes of the proposed PRTS (see the directory named src for the details.)   

--makefile which is used to compile the source code

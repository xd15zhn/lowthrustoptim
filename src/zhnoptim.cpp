#include <iostream>
#include <fstream>
#include "zhnoptim.hpp"
NAMESPACE_ZHNOPTIM_L

void Algorithm::Set_CostFunction(UserFunc *f) { _costFunc=f; }
double Algorithm::Get_Cost() const { return _MinCost; }
std::vector<double> Algorithm::Get_BestSolution() const { return _BestSolution; }

/***********************
Algorithm initialize
*@ solution: Initial solution vector
*@ solvelen: Length of solution vector
**********************/
Algorithm::Algorithm(const std::vector<double>& solution) {
    _costFunc = nullptr;
    for (int i=0; i<(int)solution.size(); ++i)
        _BestSolution.push_back(solution[i]);
    _MinCost = ZHNOPTIM_INFINITE; _TermCost = 1;
    std::cout.precision(8);
}
Algorithm::~Algorithm() {}

NAMESPACE_ZHNOPTIM_R

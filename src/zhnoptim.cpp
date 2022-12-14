#include <iostream>
#include "zhnoptim.hpp"
NAMESPACE_ZHNOPTIM_L

void Algorithm::Set_CostFunction(UserFunc *f) { _costFunc=f; }
double Algorithm::Get_Cost() const { return _MinCost; }
std::vector<double> Algorithm::Get_BestSolution() const { return _BestSolution; }

/***********************
Algorithm initialize
**********************/
Algorithm::Algorithm(const std::vector<double>& solution) {
    for (int i=0; i<(int)solution.size(); ++i)
        _BestSolution.push_back(solution[i]);
    _MinCost = ZHNOPTIM_INFINITE;
    std::cout.precision(8);
}
Algorithm::~Algorithm() {}
void Algorithm::Print_Solution() const {
    std::cout << _MinCost << ";  ";
    for (int i=0; i<_BestSolution.size(); ++i)
        std::cout << _BestSolution[i] << ", ";
    std::cout << std::endl;
}
NAMESPACE_ZHNOPTIM_R

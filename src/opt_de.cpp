#include <random>
#include <ctime>
#include "zhnoptim.hpp"
NAMESPACE_ZHNOPTIM_L

Differential_Evolution::~Differential_Evolution() {}
void Differential_Evolution::Set_Param(double F, double CR) { _F=F;_CR=CR; }

/***********************
Differential Evolution algorithm initialize
*@ solution: Initial solution vector
*@ solvelen: Length of solution vector
*@ popsize: number of the population
**********************/
Differential_Evolution::Differential_Evolution(const std::vector<double> &solution, int popsizeeach)
    : Algorithm(solution)
{
    double U, V;
    _solvelen = solution.size();
    _popsize = popsizeeach * _solvelen;
    _F = _CR = 0.5;
    _population.push_back(_BestSolution);
    for (int i = 1; i < _popsize; ++i){
        _population.push_back(_BestSolution);
        for (int j=0; j<_solvelen; ++j) {
            U = (rand()+1.0) / (RAND_MAX+1.0);
            V = (rand()+1.0) / (RAND_MAX+1.0);
            _population[i][j] = sqrt(-2.0 * log(U))* cos(6.283185307179586477 * V);
        }
    }
}

void Differential_Evolution::Initialize() {
    double fit;
    _MinCost = _costFunc->Function(_population[0]);
    _BestSolution = _population[0];
    for (short i = 1; i < _solvelen; i++) {
        fit = _costFunc->Function(_population[i]);
        if (fit < _MinCost) {
            _MinCost = fit;
            _BestSolution = _population[i];
        }
    }
}
void Differential_Evolution::Optimize_OneStep() {
    std::vector<double> mutant, offspring;
    double fit1, fit2;  // 适应度值
    int in1, in2;  // 2个差分向量的个体编号
    for (short i = 0; i < _popsize; i++) {
        /*变异:使每个个体变异产生变异个体*/
        in1 = rand() % _solvelen;
        in2 = rand() % _solvelen;
        for (short j = 0; j < _solvelen; j++) {
            _F = (rand()+1.0) / (RAND_MAX+1.0);
            mutant[j] = _population[i][j] + _F*(_population[in1][j] - _population[in2][j]);
        }
        /*交叉:变异个体与原个体交叉产生后代个体*/
        for (short j = 0; j < _solvelen; j++) {
            if ((rand()+1.0) / (RAND_MAX+1.0) > _CR)
                offspring[j] = _population[i][j];
            else
                offspring[j] = mutant[j];
        }
        /*选择:后代个体与原个体竞争并保留优胜个体*/
        fit1 = _costFunc->Function(_population[i]);
        fit2 = _costFunc->Function(offspring);
        if (fit2 < fit1) {
            for (short j = 0; j < _solvelen; j++)
                _population[i][j] = offspring[j];
            if (fit2 < _MinCost) {
                _MinCost = fit2;
                _BestSolution = _population[i];
            }
        }
    }
}

NAMESPACE_ZHNOPTIM_R

#ifndef __ZHNOPTIM_H
#define __ZHNOPTIM_H
#include <vector>

#define ZHNOPTIM_INFINITE             1e12
#define ZHNOPTIM_EPSILON              1e-12
#define ZHNOPTIM_ABS(x)               ((x)>=0?(x):-(x))
#ifndef PRINT_NAME_VALUE
#define PRINT_NAME_VALUE(x)           std::cout<<#x<<": "<<x<<std::endl;
#endif
#define ASSERT_ERROR_ALG(e, s)        if(!(e)){std::cout<<"Zhn algorithm Error: "<<s<<std::endl;abort();}
#define ASSERT_WARNING_ALG(e, s)      if(!(e)){std::cout<<"Zhn algorithm Warning: "<<s<<std::endl;}
#define NAMESPACE_ZHNOPTIM_L          namespace zhnoptim {
#define NAMESPACE_ZHNOPTIM_R          }
NAMESPACE_ZHNOPTIM_L

/**********************
Used to provide another method to replace the pointer to a function.
**********************/
class UserFunc
{
public:
    virtual double Function(const std::vector<double>& u) {return u[0];};
};

class Algorithm
{
public:
    Algorithm(const std::vector<double>& solution);
    virtual ~Algorithm();
    void Set_CostFunction(UserFunc *f);
    double Get_Cost() const;
    std::vector<double> Get_BestSolution() const;
    void Solution_Print() const;
    virtual void Initialize() = 0;
    virtual void Optimize_OneStep() = 0;
protected:
    UserFunc *_costFunc=nullptr;
    double _TermCost, _MinCost;  // Current and minimum cost value
    std::vector<double> _BestSolution;
};


/***********************
Differential Evolution(DE/best/2)
**********************/
class Differential_Evolution: public Algorithm
{
public:
    // *@ solution: Initial solution vector
    // *@ popsizeeach: 该数乘上待优化的参数数量为总的种群数量
    Differential_Evolution(const std::vector<double>& solution, int popsizeeach=15);
    virtual ~Differential_Evolution();
    virtual void Initialize() override;
    virtual void Optimize_OneStep() override;
    void Set_Param(double F, double CR);
private:
    int _solvelen;
    std::vector<std::vector<double>> _population;
    int _popsize;  // number of the population
    double _F, _CR;  // 缩放因子,交叉因子
};


/***********************
Pattern Search
**********************/
class Pattern_Search: public Algorithm{
public:
    Pattern_Search(const std::vector<double>& solution, double deltaStart=0.25);
    virtual ~Pattern_Search();
    void run();
    void Set_Param(double deltaStart);
private:
    int _solvelen;
    double _delta;
};

NAMESPACE_ZHNOPTIM_R
#endif  // ZHNOPTIM

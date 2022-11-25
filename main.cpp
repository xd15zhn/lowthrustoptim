#include <cmath>
#include <list>
#include <iostream>
#include "detector.hpp"
#include "parameters.hpp"
#include "zhnoptim.hpp"
using namespace zhnoptim;
// using namespace zhnoptim;

class Cubic: public simucpp::UserFunc
{
public:
    // Cubic(double a, double b, double c, double d):
    //     _a(a), _b(b),_c(c), _d(d) {}
    void Set(double a, double b, double c, double d)
        {_a=a; _b=b;_c=c; _d=d;}
    virtual double Function(double t) const {
        double ans = _a*t + _b;
        ans = ans*t + _c;
        ans = ans*t + _d;
        return ans;
    };
private:
    double _a, _b, _c, _d;
};

class OrbitalOptimFunc: public zhnoptim::UserFunc
{
public:
    OrbitalOptimFunc() {
        cubic1 = new Cubic;
        cubic2 = new Cubic;
        tgtr = Vector3d(177.411257719326, -106.002308439472, -6.57918956193115);
        tgtv = Vector3d(19.561773277048, 22.9483294968687, 4.68994956398122e-04);
    }
    ~OrbitalOptimFunc() {
        delete cubic1, cubic2;
    }
    virtual double Function(const std::vector<double>& solution) {
        int t=0;  // 总飞行时间
        solver._mssr->Set_InitialValue(Mat(vecdble{
            -113.101881186742, 101.033358324501, 0.00219405859350417}));
        solver._mssv->Set_InitialValue(Mat(vecdble{
            -19.3527970698371, -22.1182257893965, -0.000451706663468062}));
        solver._intM->Set_InitialValue(USV_M);
        /*分别设置第一段和第三段的发动机工作时间*/
        int dur1 = 25000/PI * atan(solution[16]*solution[16]);
        int dur2 = 22464 - 15000/PI * atan(solution[17]*solution[17]);
        /*第一段发动机工作，为俯仰角和方位角分别设置三次函数*/
        cubic1->Set(solution[0], solution[4], solution[8], solution[12]);
        cubic2->Set(solution[1], solution[5], solution[9], solution[13]);
        solver._inTheta->Set_Function(cubic1);
        solver._inPhi->Set_Function(cubic2);
        solver._cnstF->Set_OutValue(USV_F);
        for (; t < dur1; t++)
            solver._sim1.Simulate_OneStep();
        /*第二段发动机关闭*/
        solver._cnstF->Set_OutValue(0);
        for (; t < dur2; t++)
            solver._sim1.Simulate_OneStep();
        /*第三段发动机工作，为俯仰角和方位角分别设置新的三次函数*/
        cubic1->Set(solution[2], solution[6], solution[10], solution[14]);
        cubic2->Set(solution[3], solution[7], solution[11], solution[15]);
        solver._inTheta->Set_Function(cubic1);
        solver._inPhi->Set_Function(cubic2);
        solver._cnstF->Set_OutValue(USV_F);
        for (; t < 22464; t++)
            solver._sim1.Simulate_OneStep();
        /*评价此次飞行*/
        Vector3d errorr = Vector3d(solver._mssr->Get_OutValue()) - tgtr;
        Vector3d errorv = Vector3d(solver._mssv->Get_OutValue()) - tgtv;
        double err = errorr._x*errorr._x + errorr._y*errorr._y + errorr._z*errorr._z;
        err += 5*(errorv._x*errorv._x + errorv._y*errorv._y + errorv._z*errorv._z);
        err += (dur1 + 22464 - dur2) * 1e-2;
        cout << err << "\r";
        return err;
    };
private:
    Orbital_Solver solver;
    Cubic *cubic1=nullptr, *cubic2=nullptr;
    Vector3d tgtr, tgtv;
};

int main(void) {
    vecdble ans(18, 0);
    ans[12] = -1.57; ans[16] = 1;
    OrbitalOptimFunc *func1 = new OrbitalOptimFunc;
    vecdble bestans;
    double bestcost;
    Differential_Evolution alg(ans, 10);
    alg.Set_CostFunction(func1);
    alg.Initialize();
    cout << "First step finished." << endl;
    while (1) {
        alg.Optimize_OneStep();
        bestcost = alg.Get_Cost();
    }
    return 0;
}

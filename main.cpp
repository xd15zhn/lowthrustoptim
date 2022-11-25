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
        solver._sim1.Set_t(0);
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
        cubic2->Set(1e-5*solution[1], 1e-5*solution[5], 0.01*solution[9], 0.01*solution[13]);
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
        cubic2->Set(1e-5*solution[3], 1e-5*solution[7], 0.01*solution[11], 0.01*solution[15]);
        solver._inTheta->Set_Function(cubic1);
        solver._inPhi->Set_Function(cubic2);
        solver._cnstF->Set_OutValue(USV_F);
        for (; t < 22464; t++)
            solver._sim1.Simulate_OneStep();
        /*评价此次飞行*/
        Vector3d errorr = Vector3d(solver._mssr->Get_OutValue()) - tgtr;
        Vector3d errorv = Vector3d(solver._mssv->Get_OutValue()) - tgtv;
        double err = errorr._x*errorr._x + errorr._y*errorr._y + errorr._z*errorr._z;
        err += 1*(errorv._x*errorv._x + errorv._y*errorv._y + errorv._z*errorv._z);
        err += (dur1 + 22464 - dur2) * 1e-3;
        return err;
    };
private:
    Orbital_Solver solver;
    Cubic *cubic1=nullptr, *cubic2=nullptr;
    Vector3d tgtr, tgtv;
};

#define USE_PS
int main(void) {
    // vecdble ans(18, 0);
    // ans[12] = -1.57; ans[16] = 1;
    vecdble ans{
        // 0.327971, 2.81796, -0.523281, -0.611048, -1.17108, 1.57724,
        // 0.510948, -0.712889, 0.256603, -0.226293, 0.760507, -0.214071,
        // -1.51277, -0.300987, 0.773741, -1.33875, 1.19229, -0.752097
        0.327971, 1.06796, -0.023281, -1.361048, -1.17108, -0.17276,
        0.760948, -0.462889, 0.256603, -1.976293, 0.510507, -0.964071,
        -1.76277, -1.800987, 1.273741, -1.83875, 0.69229, -0.502097
    };  // 302.564
    cout.precision(6);
    OrbitalOptimFunc *func1 = new OrbitalOptimFunc;
#ifdef USE_PS
    Pattern_Search alg(ans);
    alg.Set_CostFunction(func1);
    alg.Initialize();
    alg.run();
#else
    int epoch=0;
    Differential_Evolution alg(ans, 10);
    alg.Set_CostFunction(func1);
    alg.Push_Solvsion(vecdble{
        0.287717, 2.77068, 1.31148, 1.01602, -0.610206, 1.1287,
        0.482372, 2.97458, -0.75428, 2.42673, 0.91558, 0.618227,
        -1.53584, -0.00251363, -1.56059, -1.15623, 0.97752, -0.484669,
    });
    alg.Push_Solvsion(vecdble{
        0.327971, 2.81796, -0.523281, -0.543126, -1.17108, 1.62028,
        0.510948, -0.712889, 0.312796, -0.226293, 0.545372, -0.420163,
        -1.51195, -0.757546, 0.773741, -1.68666, 1.19229, -0.892871,
    });
    alg.Push_Solvsion(vecdble{
        0.287717, 2.77068, 1.31148, 1.01602, -0.610206, 1.1287,
        0.482372, 2.97458, -0.75428, 2.42673, 0.91558, 0.618227,
        -1.53584, -0.00251363, -1.56059, -1.15623, 0.97752, -0.484669,
    });
    alg.Push_Solvsion(vecdble{
        0.327971, 2.81796, -0.523281, -0.543126, -1.17108, 2.10536,
        0.510948, -0.322343, 0.312796, 0.63923, 0.545372, -0.0212892,
        -1.51195, -0.783517, 0.773741, -1.68666, 1.38849, -0.892871
    });
    alg.Push_Solvsion(vecdble{
        0.327971, 2.81796, -0.523281, -0.543126, -1.17108, 1.62028,
        0.510948, -0.712889, 0.312796, -0.226293, 0.545372, -0.420163,
        -1.51195, -0.757546, 0.773741, -1.68666, 1.19229, -0.892871
    });
    alg.Push_Solvsion(vecdble{
        0.327971, 2.81796, -0.523281, -0.611048, -1.17108, 1.57724,
        0.510948, -0.712889, 0.256603, -0.226293, 0.760507, -0.214071,
        -1.51277, -0.300987, 0.773741, -1.33875, 1.19229, -0.752097,
    });
    alg.Initialize();
    cout << "\nFirst step finished." << endl;
    while (epoch < 5) {
        cout << "epoch:" << epoch << ",  best solution:        " << endl;
        alg.Print_Solution();
        alg.Optimize_OneStep();
        epoch++;
    }
#endif
    return 0;
}

#include <vector>
#include <complex>

using namespace std;

#ifndef MY_MATH
#define MY_MATH

using Cplx = complex<double>;
using Point = vector<double>;
using Triangle = vector<Point>;

//const double PI = 3.141592653589793;
const double M_C = 299792458;
const double M_MU0 = 4 * M_PI * 1e-7;
const double M_EPS0 = 8.854187817620e-12;

const vector<double> M_XI1 ({
  0.333333333333333,
  0.059715871789770,
  0.470142064105115,
  0.470142064105115,
  0.797426985353087,
  0.101286507323456,
  0.101286507323456
});
const vector<double> M_XI2 ({
  0.333333333333333,
  0.470142064105115,
  0.059715871789770,
  0.470142064105115,
  0.101286507323456,
  0.797426985353087,
  0.101286507323456
});
const vector<double> M_XI3 ({
  0.333333333333333,
  0.470142064105115,
  0.470142064105115,
  0.059715871789770,
  0.101286507323456,
  0.101286507323456,
  0.797426985353087
});
const vector<double> M_WEIGHT ({
  0.112500000000000,
  0.066197076394253,
  0.066197076394253,
  0.066197076394253,
  0.062969590272414,
  0.062969590272414,
  0.062969590272414
});

vector<double> xitor (double xi1, double xi2, double xi3, vector<double> r1, vector<double> r2, vector<double> r3);

complex<double> greenFunc (vector<double> rm, vector<double> rn, double k);
complex<double> greenFuncExtSing (vector<double> rm, vector<double> rn, double k);

double distanceTwoPoint (vector<double> p1, vector<double> p2);
double distance2 (vector<double> p1, vector<double> p2);
double triangleArea (vector<double> v1, vector<double> v2, vector<double> v3 );
vector<double> triangleCenter (vector<double> v1, vector<double> v2, vector<double> v3 );

double dGaussianQuadTriangle (vector<double> f, double area);
complex<double> cGaussianQuadTriangle (vector<complex<double> > f, double area);
vector<complex<double>> cGaussianQuadTriangleVect (vector<vector<complex<double> > > f, double area);


double ComputeIntegralType1Arcioni (const Triangle Tm);
double ComputeIntegralType1EibertHansen (const Triangle Tm);
double ComputeIntegralType1GqtAndWilton (const Triangle Tm, const Triangle Tn);
double ComputeIntegralType1GqtAndOijala (const Triangle Tm, const Triangle Tn);
double ComputeIntegralType1Gqt (const Triangle Tm, const Triangle Tn);

double ComputeIntegralType2Arcioni (const Triangle Tm, int isalphaeqbeta);
double ComputeIntegralType2EibertHansen (const Triangle Tm, vector<double> vm, vector<double> vn);
double ComputeIntegralType2GqtAndWilton(const Triangle Tm, const Triangle Tn, vector<double> vm, vector<double> vn);
double ComputeIntegralType2Gqt (const Triangle Tm, const Triangle Tn, vector<double> vm, vector<double> vn);

#endif

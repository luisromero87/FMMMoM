#include <vector>
#include <cmath>

using namespace std;

#include "my_math.h"

vector<double> xitor (double M_XI1, double M_XI2, double M_XI3, vector<double> r1, vector<double> r2, vector<double> r3) {
	vector<double> r (3);
	for (int i = 0; i < 3; i++) r[i] = M_XI1 * r1[i] + M_XI2 * r2[i] + M_XI3 * r3[i];
	return r;
}

complex<double> greenFunc (vector<double> rm, vector<double> rn, double k) {
	complex<double> G = 0;
	double R = distanceTwoPoint(rm, rn);
	G = exp( complex<double> (0, -1) * k * R) / (4. * M_PI * R);
	return G;
}

complex<double> greenFuncExtSing (vector<double> rm, vector<double> rn, double k) {
	complex<double> G = 0;
	double R = distanceTwoPoint(rm, rn);
	G = (exp( complex<double> (0, -1) * k * R) - 1.0) / (4. * M_PI* R);
	return G;
}

double distanceTwoPoint (vector<double> p1, vector<double> p2) {
	double R = 0;
	for (int i = 0; i < 3; i++) R += pow(p1[i] - p2[i], 2);
	R = sqrt(R);
	return R;
}

double distance2 (vector<double> p1, vector<double> p2) {
	double d2 = 0;
	for (int i = 0; i < 3; i++) d2 += pow(p1[i] - p2[i], 2);
	return d2;
}

double triangleArea (vector<double> v1, vector<double> v2, vector<double> v3 ) {
	double x1, x2, x3, y1, y2, y3;
	x1 = v1[0] - v2[0];
	x2 = v1[1] - v2[1];
	x3 = v1[2] - v2[2];
	y1 = v1[0] - v3[0];
	y2 = v1[1] - v3[1];
	y3 = v1[2] - v3[2];
	return 0.5 * sqrt(pow(x2 * y3 - x3 * y2, 2) + pow(x3 * y1 - x1 * y3, 2) + pow(x1 * y2 - x2 * y1, 2));
};

vector< double > triangleCenter (vector<double> v1, vector<double> v2, vector<double> v3 ) {
	vector< double > center (3, 0);
	for (int i = 0; i < 3; i++) center[i] = (v1[i] + v2[i] + v3[i]) / 3;
	return center;
};

double dGaussianQuadTriangle (vector<double> f, double area) {
	double sum = 0;
	for (int n = 0; n < 7; n++) sum += M_WEIGHT[n] * f[n];
	return 2.0 * area * sum;
}

complex<double> cGaussianQuadTriangle (vector<complex<double> > f, double area) {
	complex<double> sum = 0;
	for (int n = 0; n < 7; n++) {
		sum += M_WEIGHT[n] * f[n];
	}
	return 2.0 * area * sum;
}

vector<complex<double>> cGaussianQuadTriangleVect (vector<vector<complex<double> > > f, double area) {
	vector<complex<double>> sum (3, 0);
	for (int n = 0; n < 7; n++) {
		for (int i = 0; i < 3; i++) sum[i] += M_WEIGHT[n] * f[n][i];
	}
	for (int i = 0; i < 3; i++) sum[i] *= 2 * area;
	return sum;
}


//int(int(int(int(1/sqrt((x1-x2)**2-(y1-y2)**2),x1,1/sqrt(3)y1,1-1/sqrt(3)y1),y1,0,sqrt(3)/2),x2,-1/sqrt(3)y2,1/sqrt(3)y2-1),y2,0,-sqrt(3)/2)

//Integrate[1/Sqrt[(x1-x2)^2-(y1-y2)^2],{x1,1/Sqrt[3]y1,1-1/Sqrt[3]y1},{y1,0,Sqrt[3]/2},{x2,-1/Sqrt[3]y2,1/Sqrt[3]y2-1},{y2,0,-Sqrt[3]/2}]

// Compute the integral of the form:
// \Int_{T} \rho \cdot \Int_{T'} \rho' /R d^2r' d^2r
// R = |r'-r|
// when the two triangles overlap using Arcioni's formula.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//
// P. Arcioni, M. Bressan and L. Perregrini, “On the evaluation of double surface integrals
// arising in the application of the boundary integral method to 3-D problems”, 1997
double ComputeIntegralType2Arcioni (const Triangle Tm, int isalphaeqbeta) {
	double int_type_2 = 0;
	double a, b, c, p, area;
	a = distanceTwoPoint (Tm[0], Tm[1]);
	b = distanceTwoPoint (Tm[0], Tm[2]);
	c = distanceTwoPoint (Tm[1], Tm[2]);
	p = (a + b + c) / 2;
	area = triangleArea(Tm[0], Tm[1], Tm[2]);
	if (isalphaeqbeta) {
		int_type_2 = area * area / (120. * M_PI) * (
		                 (10 + 3 * (c * c - a * a) / (b * b) - 3 * (a * a - b * b) / (c * c)) * a -
		                 (5 - 3 * (a * a - b * b) / (c * c) - 2 * (b * b - c * c) / (a * a)) * b -
		                 (5 + 3 * (c * c - a * a) / (b * b) + 2 * (b * b - c * c) / (a * a)) * c +
		                 (a * a - 3 * b * b - 3 * c * c - 8 * area * area / (a * a)) * 2 / a * log(1 - a / p) +
		                 (a * a - 2 * b * b - 4 * c * c + 6 * area * area / (b * b)) * 4 / b * log(1 - b / p) +
		                 (a * a - 4 * b * b - 2 * c * c + 6 * area * area / (c * c)) * 4 / c * log(1 - c / p));
	}
	else {
		int_type_2 = area * area / (240. * M_PI) * (
		                 (-10 + (c * c - a * a) / (b * b) - (a * a - b * b) / (c * c)) * a +
		                 (5 + (a * a - b * b) / (c * c) - 6 * (b * b - c * c) / (a * a)) * b +
		                 (5 - (c * c - a * a) / (b * b) + 6 * (b * b - c * c) / (a * a)) * c +
		                 (2 * a * a - b * b - c * c + 4 * area * area / (a * a)) * 12 / a * log(1 - a / p) +
		                 (9 * a * a - 3 * b * b - c * c + 4 * area * area / (b * b)) * 2 / b * log(1 - b / p) +
		                 (9 * a * a - b * b - 3 * c * c + 4 * area * area / (c * c)) * 2 / c * log(1 - c / p));
	}
	return int_type_2;
}
// Compute the integral of the form:
// \Int_{T} \rho \cdot \Int_{T'} \rho' /R d^2r' d^2r
// R = |r'-r|
// when the two triangles overlap using Eibert-Hansen's formula.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
double ComputeIntegralType2EibertHansen (const Triangle Tm, vector<double> vm, vector<double> vn) {
	vector<double> v1 (Tm[1]);
	vector<double> v2 (Tm[0]);
	vector<double> v3 (Tm[2]);

	double area2 = pow(triangleArea(Tm[0], Tm[1], Tm[2]), 2);

	double l1, l2, l3;

	l1 = distanceTwoPoint(v2, v3);
	l2 = distanceTwoPoint(v3, v1);
	l3 = distanceTwoPoint(v1, v2);

	double ln1 = log(abs( (pow(l1 + l2, 2) - l3 * l3) / (l2 * l2 - pow(l3 - l1, 2)) ));
	double ln2 = log(abs( (pow(l2 + l3, 2) - l1 * l1) / (l3 * l3 - pow(l1 - l2, 2)) ));
	double ln3 = log(abs( (pow(l3 + l1, 2) - l2 * l2) / (l1 * l1 - pow(l1 - l3, 2)) ));

	double lambda2lambda2 = 1 / (20 * l1) * ln1 +
	                        (l1 * l1 + 5 * l2 * l2 - l3 * l3) / (120 * l2 * l2 * l2) * ln2 +
	                        (l1 * l1 - l2 * l2 + 5 * l3 * l3) / (120 * l3 * l3 * l3) * ln3 +
	                        (l3 - l1) / (60 * l2 * l2) +
	                        (l2 - l1) / (60 * l3 * l3);
	double lambda1lambda2 = (3 * l1 * l1 + l2 * l2 - l3 * l3) / (80 * l1 * l1 * l1) * ln1 +
	                        (l1 * l1 + 3 * l2 * l2 - l3 * l3) / (80 * l2 * l2 * l2) * ln2 +
	                        1 / (40 * l3) * ln3 +
	                        (l3 - l2) / (40 * l1 * l1) +
	                        (l3 - l1) / (40 * l2 * l2);
	double lambda2 = 1 / (8 * l1) * ln1 +
	                 (l1 * l1 + 5 * l2 * l2 - l3 * l3) / (48 * l2 * l2 * l2) * ln2 +
	                 (l1 * l1 - l2 * l2 + 5 * l3 * l3) / (48 * l3 * l3 * l3) * ln3 +
	                 (l3 - l1) / (24 * l2 * l2) +
	                 (l2 - l1) / (24 * l3 * l3);

	v1.swap(v2);
	l1 = distanceTwoPoint(v2, v3);
	l2 = distanceTwoPoint(v3, v1);
	l3 = distanceTwoPoint(v1, v2);

	ln1 = log(abs( (pow(l1 + l2, 2) - l3 * l3) / (l2 * l2 - pow(l3 - l1, 2)) ));
	ln2 = log(abs( (pow(l2 + l3, 2) - l1 * l1) / (l3 * l3 - pow(l1 - l2, 2)) ));
	ln3 = log(abs( (pow(l3 + l1, 2) - l2 * l2) / (l1 * l1 - pow(l1 - l3, 2)) ));

	double lambda1lambda1 = 1 / (20 * l1) * ln1 +
	                        (l1 * l1 + 5 * l2 * l2 - l3 * l3) / (120 * l2 * l2 * l2) * ln2 +
	                        (l1 * l1 - l2 * l2 + 5 * l3 * l3) / (120 * l3 * l3 * l3) * ln3 +
	                        (l3 - l1) / (60 * l2 * l2) +
	                        (l2 - l1) / (60 * l3 * l3);
	double lambda2lambda1 = (3 * l1 * l1 + l2 * l2 - l3 * l3) / (80 * l1 * l1 * l1) * ln1 +
	                        (l1 * l1 + 3 * l2 * l2 - l3 * l3) / (80 * l2 * l2 * l2) * ln2 +
	                        1 / (40 * l3) * ln3 +
	                        (l3 - l2) / (40 * l1 * l1) +
	                        (l3 - l1) / (40 * l2 * l2);
	double lambda1 = 1 / (8 * l1) * ln1 +
	                 (l1 * l1 + 5 * l2 * l2 - l3 * l3) / (48 * l2 * l2 * l2) * ln2 +
	                 (l1 * l1 - l2 * l2 + 5 * l3 * l3) / (48 * l3 * l3 * l3) * ln3 +
	                 (l3 - l1) / (24 * l2 * l2) +
	                 (l2 - l1) / (24 * l3 * l3);
	double intr =  1 / (3 * l1) * ln1 +
	               1 / (3 * l2) * ln2 +
	               1 / (3 * l3) * ln3;

	double v1dotv1 = 0;
	double v2dotv2 = 0;
	double v3dotv3 = 0;
	double v1dotv2 = 0;
	double v1dotv3 = 0;
	double v2dotv3 = 0;
	double v1dotvn = 0;
	double v2dotvn = 0;
	double v3dotvn = 0;
	double v1dotvm = 0;
	double v2dotvm = 0;
	double v3dotvm = 0;
	double vndotvm = 0;
	for (int i = 0; i < 3; i++) {
		v1dotv1 += v1[i] * v1[i];
		v2dotv2 += v2[i] * v2[i];
		v3dotv3 += v3[i] * v3[i];
		v1dotv2 += v1[i] * v2[i];
		v1dotv3 += v1[i] * v3[i];
		v2dotv3 += v2[i] * v3[i];
		v1dotvn += v1[i] * vn[i];
		v2dotvn += v2[i] * vn[i];
		v3dotvn += v3[i] * vn[i];
		v1dotvm += v1[i] * vm[i];
		v2dotvm += v2[i] * vm[i];
		v3dotvm += v3[i] * vm[i];
		vndotvm += vn[i] * vm[i];
	}

	double int_type_2 = lambda1lambda1 * (v1dotv1 - 2 * v1dotv3 + v3dotv3) +
	                    lambda2lambda2 * (v2dotv2 - 2 * v2dotv3 + v3dotv3) +
	                    lambda1lambda2 * (v1dotv2 - v1dotv3 - v2dotv3 + v3dotv3) +
	                    lambda2lambda1 * (v1dotv2 - v2dotv3 - v1dotv3 + v3dotv3) +
	                    lambda1 * (v1dotv3 - v3dotv3 - v1dotvn + v3dotvn) +
	                    lambda2 * (v2dotv3 - v3dotv3 - v2dotvn + v3dotvn) +
	                    lambda1 * (v1dotv3 - v3dotv3 - v1dotvm + v3dotvm) +
	                    lambda2 * (v2dotv3 - v3dotv3 - v2dotvm + v3dotvm) +
	                    intr * (v3dotv3 - v3dotvn - v3dotvm + vndotvm);


	int_type_2 *= area2 / M_PI;
	return int_type_2;
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles are close to each other.
// The inner integral with Wilton's formula and the outter integral with Gaussian quadrature
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//  * Tn : The source triangle. A tuple with the triangles vertices.
double ComputeIntegralType2GqtAndWilton(const Triangle Tm, const Triangle Tn, vector<double> vm, vector<double> vn) {
	vector<double> u_vect(3, 0);
	vector<double> n_vect(3, 0);
	vector<double> l_vect(3, 0);
	vector<double> rm(3, 0);
	vector<double> rho(3, 0);
	vector<double> rhon(3, 0);
	vector<double> rhop(3, 0);
	vector<double> rho_n(3, 0);
	vector<double> P0_vect(3, 0);
	double lp, ln, P0, Pp, Pn, d, R0, Rp, Rn, P0puntou;
	double npuntor = 0;
	vector<double> int1 (3, 0);
	double int1re = 0;
	double int2 = 0;
	vector<vector<double>> int2_vect (7, vector<double>(3, 0));
	vector<double> f2(7, 0);

	for (int i = 0; i < 3; i++) u_vect[i] = Tn[1][i] - Tn[0][i];
	for (int i = 0; i < 3; i++) n_vect[i] = Tn[2][i] - Tn[0][i];
	n_vect = {u_vect[1]*n_vect[2] - u_vect[2]*n_vect[1],
	          u_vect[2]*n_vect[0] - u_vect[0]*n_vect[2],
	          u_vect[0]*n_vect[1] - u_vect[1]*n_vect[0]
	         };
	npuntor = distanceTwoPoint(Point ({0, 0, 0}), n_vect);
	for (int i = 0; i < 3; i++) n_vect[i] /= npuntor;
	npuntor = 0;
	for (int i = 0; i < 3; i++) npuntor += n_vect[i] * vn[i];
	for (int i = 0; i < 3; i++) rho_n[i] = vn[i] - n_vect[i] * npuntor;

	for (int m = 0; m < 7; m++) {
		rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
		npuntor = 0;
		for (int i = 0; i < 3; i++) npuntor += n_vect[i] * rm[i];
		for (int i = 0; i < 3; i++) rho[i] = rm[i] - n_vect[i] * npuntor;

		int2 = 0;
		for (int i = 0; i < 3; i++) int1[i] = 0;
		for (int l = 0; l < 3; l++) {
			npuntor = 0;
			for (int i = 0; i < 3; i++) npuntor += n_vect[i] * Tn[(l + 1) % 3][i];
			for (int i = 0; i < 3; i++) rhop[i] = Tn[(l + 1) % 3][i] - n_vect[i] * npuntor;
			npuntor = 0;
			for (int i = 0; i < 3; i++) npuntor += n_vect[i] * Tn[l][i];
			for (int i = 0; i < 3; i++) rhon[i] = Tn[l][i] - n_vect[i] * npuntor;

			for (int i = 0; i < 3; i++) l_vect[i] = rhop[i] - rhon[i];
			npuntor = distanceTwoPoint(Point ({0, 0, 0}), l_vect);
			for (int i = 0; i < 3; i++) l_vect[i] /= npuntor;

			u_vect = {l_vect[1]*n_vect[2] - l_vect[2]*n_vect[1],
			          l_vect[2]*n_vect[0] - l_vect[0]*n_vect[2],
			          l_vect[0]*n_vect[1] - l_vect[1]*n_vect[0]
			         };

			lp = 0;
			for (int i = 0; i < 3; i++) lp += (rhop[i] - rho[i]) * l_vect[i];
			ln = 0;
			for (int i = 0; i < 3; i++) ln += (rhon[i] - rho[i]) * l_vect[i];

			P0 = 0;
			for (int i = 0; i < 3; i++) P0 += (rhop[i] - rho[i]) * u_vect[i];
			P0 = abs(P0);

			Pp = sqrt(P0 * P0 + lp * lp);
			Pn = sqrt(P0 * P0 + ln * ln);

			for (int i = 0; i < 3; i++) P0_vect[i] = (rhop[i] - rho[i] - lp * l_vect[i]) / P0;

			d = 0;
			for (int i = 0; i < 3; i++) d += n_vect[i] * (rm[i] - Tn[l][i]);

			R0 = sqrt(P0 * P0 + d * d);
			Rp = sqrt(Pp * Pp + d * d);
			Rn = sqrt(Pn * Pn + d * d);

			int1re = R0 * R0 * log((Rp + lp) / (Rn + ln)) + lp * Rp - ln * Rn;
			for (int i = 0; i < 3; i++) int1[i] += u_vect[i] * int1re;

			P0puntou = 0;
			for (int i = 0; i < 3; i++) P0puntou += P0_vect[i] * u_vect[i];
			int2 += P0puntou * (P0 * log((Rp + lp) / (Rn + ln)) - abs(d) * (
			                        atan((P0 * lp) / (R0 * R0 + abs(d) * Rp)) -
			                        atan((P0 * ln) / (R0 * R0 + abs(d) * Rn))));
		}

		f2[m] = 0;
		for (int i = 0; i < 3; i++)
			f2[m] += (0.5 * int1[i] + int2 * (rho[i] - rho_n[i])) * (rm[i] - vm[i]);
	}
	return dGaussianQuadTriangle(f2, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles are far enough.
// Both integrals are performed with Gaussian quadrature.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//  * Tn : The source triangle. A tuple with the triangles vertices.
double ComputeIntegralType2Gqt (const Triangle Tm, const Triangle Tn, vector<double> vm, vector<double> vn) {
	vector<double> f1(7, 0);
	vector<double> f2(7, 0);
	vector<double> rn(3, 0);
	vector<double> rm(3, 0);
	double temp_double;
	for (int m = 0; m < 7; m++) {
		rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
		for (int n = 0 ; n < 7; n++) {
			rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
			temp_double = distanceTwoPoint(rm, rn);
			f1[n] = 0;
			for (int i = 0; i < 3; i++) f1[n] += (rm[i] - vm[i]) * (rn[i] - vn[i]);
			f1[n] *= 1 / temp_double;
		}
		f2[m] = dGaussianQuadTriangle(f1, triangleArea(Tn[0], Tn[1], Tn[2]));
	}
	return dGaussianQuadTriangle(f2, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
}

// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles overlap using Arcioni's formula.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//
// P. Arcioni, M. Bressan and L. Perregrini, “On the evaluation of double surface integrals
// arising in the application of the boundary integral method to 3-D problems”, 1997
double ComputeIntegralType1Arcioni (const Triangle Tm) {
	double a = distanceTwoPoint (Tm[0], Tm[1]);
	double b = distanceTwoPoint (Tm[0], Tm[2]);
	double c = distanceTwoPoint (Tm[1], Tm[2]);
	double p = (a + b + c) / 2;
	return -pow(triangleArea(Tm[0], Tm[1], Tm[2]), 2.) / (3.*M_PI) * (
	           1. / a * log(1. - a / p) +
	           1. / b * log(1. - b / p) +
	           1. / c * log(1. - c / p));
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles overlap using Eibert-Hansen's formula.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
double ComputeIntegralType1EibertHansen (const Triangle Tm) {
	double l1 = distanceTwoPoint(Tm[1], Tm[2]);
	double l2 = distanceTwoPoint(Tm[2], Tm[0]);
	double l3 = distanceTwoPoint(Tm[0], Tm[1]);
	double int_type_1 = 1 / (3 * l1) * log(abs((pow(l1 + l2, 2) - l3 * l3) / (l2 * l2 - pow(l3 - l1, 2))))
	                    + 1 / (3 * l2) * log(abs((pow(l2 + l3, 2) - l1 * l1) / (l3 * l3 - pow(l1 - l2, 2))))
	                    + 1 / (3 * l3) * log(abs((pow(l3 + l1, 2) - l2 * l2) / (l1 * l1 - pow(l2 - l3, 2))));
	int_type_1 *= pow(triangleArea(Tm[0], Tm[1], Tm[2]), 2) / M_PI;
	return int_type_1;
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles are close to each other.
// The inner integral with Wilton's formula and the outter integral with Gaussian quadrature
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//  * Tn : The source triangle. A tuple with the triangles vertices.
double ComputeIntegralType1GqtAndWilton(const Triangle Tm, const Triangle Tn) {
	vector<double> u_vect(3, 0);
	vector<double> n_vect(3, 0);
	vector<double> l_vect(3, 0);
	vector<double> rm(3, 0);
	vector<double> rho(3, 0);
	vector<double> rhon(3, 0);
	vector<double> rhop(3, 0);
	vector<double> P0_vect(3, 0);
	double lp, ln, P0, Pp, Pn, d, R0, Rp, Rn, P0puntou;
	double npuntor = 0;
	double int2 = 0;
	vector<double> f2(7, 0);

	for (int i = 0; i < 3; i++) u_vect[i] = Tn[1][i] - Tn[0][i];
	for (int i = 0; i < 3; i++) n_vect[i] = Tn[2][i] - Tn[0][i];
	n_vect = {u_vect[1]*n_vect[2] - u_vect[2]*n_vect[1],
	          u_vect[2]*n_vect[0] - u_vect[0]*n_vect[2],
	          u_vect[0]*n_vect[1] - u_vect[1]*n_vect[0]
	         };
	npuntor = distanceTwoPoint(Point ({0, 0, 0}), n_vect);
	for (int i = 0; i < 3; i++) n_vect[i] /= npuntor;

	for (int m = 0; m < 7; m++) {
		rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
		npuntor = 0;
		for (int i = 0; i < 3; i++) npuntor += n_vect[i] * rm[i];
		for (int i = 0; i < 3; i++) rho[i] = rm[i] - n_vect[i] * npuntor;

		int2 = 0;
		for (int l = 0; l < 3; l++) {
			npuntor = 0;
			for (int i = 0; i < 3; i++) npuntor += n_vect[i] * Tn[(l + 1) % 3][i];
			for (int i = 0; i < 3; i++) rhop[i] = Tn[(l + 1) % 3][i] - n_vect[i] * npuntor;
			npuntor = 0;
			for (int i = 0; i < 3; i++) npuntor += n_vect[i] * Tn[l][i];
			for (int i = 0; i < 3; i++) rhon[i] = Tn[l][i] - n_vect[i] * npuntor;

			for (int i = 0; i < 3; i++) l_vect[i] = rhop[i] - rhon[i];
			npuntor = distanceTwoPoint(Point ({0, 0, 0}), l_vect);
			for (int i = 0; i < 3; i++) l_vect[i] /= npuntor;

			u_vect = {l_vect[1]*n_vect[2] - l_vect[2]*n_vect[1],
			          l_vect[2]*n_vect[0] - l_vect[0]*n_vect[2],
			          l_vect[0]*n_vect[1] - l_vect[1]*n_vect[0]
			         };

			lp = 0;
			for (int i = 0; i < 3; i++) lp += (rhop[i] - rho[i]) * l_vect[i];
			ln = 0;
			for (int i = 0; i < 3; i++) ln += (rhon[i] - rho[i]) * l_vect[i];

			P0 = 0;
			for (int i = 0; i < 3; i++) P0 += (rhop[i] - rho[i]) * u_vect[i];
			P0 = abs(P0);

			Pp = sqrt(P0 * P0 + lp * lp);
			Pn = sqrt(P0 * P0 + ln * ln);

			for (int i = 0; i < 3; i++) P0_vect[i] = (rhop[i] - rho[i] - lp * l_vect[i]) / P0;

			d = 0;
			for (int i = 0; i < 3; i++) d += n_vect[i] * (rm[i] - Tn[l][i]);

			R0 = sqrt(P0 * P0 + d * d);
			Rp = sqrt(Pp * Pp + d * d);
			Rn = sqrt(Pn * Pn + d * d);

			P0puntou = 0;

			for (int i = 0; i < 3; i++) P0puntou += P0_vect[i] * u_vect[i];
			int2 += P0puntou * (P0 * log((Rp + lp) / (Rn + ln)) - abs(d) * (
			                        atan((P0 * lp) / (R0 * R0 + abs(d) * Rp)) -
			                        atan((P0 * ln) / (R0 * R0 + abs(d) * Rn))));
		}
		f2[m] = int2;
	}
	return dGaussianQuadTriangle(f2, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles are close to each other.
// The inner integral with Oijala's formula and the outter integral with Gaussian quadrature
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//  * Tn : The source triangle. A tuple with the triangles vertices.
double ComputeIntegralType1GqtAndOijala (const Triangle Tm, const Triangle Tn) {
	vector<double> u_vect(3, 0);
	vector<double> n_vect(3, 0);
	vector<double> v_vect(3, 0);
	vector<double> rm(3, 0);
	double temp_double;
	double u3, v3, w0, v0, u0;
	double s1n, s2n, s3n, R1n, R2n, R3n;
	double s1p, s2p, s3p, R1p, R2p, R3p;
	double t01, t02, t03;
	double R01, R02, R03;
	vector<double> beta (3, 0);
	vector<double> Iminus1 (3, 0);
	double K1minus1 = 0, K1minus3 = 0;
	vector<double> f2 (7, 0);
	double l1 = distanceTwoPoint(Tn[1], Tn[2]);
	double l2 = distanceTwoPoint(Tn[0], Tn[2]);
	double l3 = distanceTwoPoint(Tn[0], Tn[1]);
	for (int i = 0; i < 3; i++) u_vect[i] = (Tn[1][i] - Tn[0][i]) / l3;
	for (int i = 0; i < 3; i++) n_vect[i] = (Tn[2][i] - Tn[0][i]) / l3;
	n_vect = {u_vect[1]*n_vect[2] - u_vect[2]*n_vect[1],
	          u_vect[2]*n_vect[0] - u_vect[0]*n_vect[2],
	          u_vect[0]*n_vect[1] - u_vect[1]*n_vect[0]
	         };
	temp_double = distanceTwoPoint(Point ({0, 0, 0}), n_vect);
	for (int i = 0; i < 3; i++) n_vect[i] /= temp_double;
	v_vect = {n_vect[1]*u_vect[2] - n_vect[2]*u_vect[1],
	          n_vect[2]*u_vect[0] - n_vect[0]*u_vect[2],
	          n_vect[0]*u_vect[1] - n_vect[1]*u_vect[0]
	         };
	u3 = 0;
	for (int i = 0; i < 3; i++) u3 += (Tn[2][i] - Tn[0][i]) * u_vect[i];

	v3 = 2 * triangleArea(Tn[0], Tn[1], Tn[2]) / l3;

	for (int m = 0; m < 7; m++) {
		rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);

		w0 = 0;
		for (int i = 0; i < 3; i++) w0 += (rm[i] - Tn[0][i]) * n_vect[i];
		v0 = 0;
		for (int i = 0; i < 3; i++) v0 += (rm[i] - Tn[0][i]) * v_vect[i];
		u0 = 0;
		for (int i = 0; i < 3; i++) u0 += (rm[i] - Tn[0][i]) * u_vect[i];

		s1n = -((l3 - u0) * (l3 - u3) + v0 * v3) / l1;
		s1p = s1n + l1;
		s2n = -(u3 * (u3 - u0) + v3 * (v3 - v0)) / l2;
		s2p = s2n + l2;
		s3n = -u0;
		s3p = s3n + l3;

		t01 = (v0 * (u3 - l3) + v3 * (l3 - u0)) / l1;
		t02 = (u0 * v3 - v0 * u3) / l2;
		t03 = v0;

		R01 = sqrt(t01 * t01 + w0 * w0);
		R02 = sqrt(t02 * t02 + w0 * w0);
		R03 = sqrt(t03 * t03 + w0 * w0);

		R1p = R2n =  distanceTwoPoint(rm, Tn[2]);
		R2p = R3n =  distanceTwoPoint(rm, Tn[0]);
		R3p = R1n =  distanceTwoPoint(rm, Tn[1]);

		beta[0] = atan((t01 * s1p) / (R01 * R01 + abs(w0) * R1p)) - atan((t01 * s1n) / (R01 * R01 + abs(w0) * R1n));
		beta[1] = atan((t02 * s2p) / (R02 * R02 + abs(w0) * R2p)) - atan((t02 * s2n) / (R02 * R02 + abs(w0) * R2n));
		beta[2] = atan((t03 * s3p) / (R03 * R03 + abs(w0) * R3p)) - atan((t03 * s3n) / (R03 * R03 + abs(w0) * R3n));

		Iminus1[0] = log((R1p + s1p) / (R1n + s1n));
		Iminus1[1] = log((R2p + s2p) / (R2n + s2n));
		Iminus1[2] = log((R3p + s3p) / (R3n + s3n));

		K1minus1 = t01 * Iminus1[0] + t02 * Iminus1[1] + t03 * Iminus1[3];
		if (abs(w0) > 1e-3) {
			K1minus3 = 1 / abs(w0) * (beta[0] + beta[1] + beta[2]);
			K1minus1 += -w0 * w0 * K1minus3;
		}
		f2[m] = K1minus1;
	}
	return dGaussianQuadTriangle(f2, triangleArea(Tm[0], Tm[1], Tm[2]));
}
// Compute the integral of the form:
// \Int_{T} \Int_{T'} 1/R d^2r' d^2r
// R = |r'-r|
// when the two triangles are far enough.
// Both integrals are performed with Gaussian quadrature.
//
// Params
//  * Tm : The observation triangle. A tuple with the triangles vertices.
//  * Tn : The source triangle. A tuple with the triangles vertices.
double ComputeIntegralType1Gqt (const Triangle Tm, const Triangle Tn) {
	vector<double> f1(7, 0);
	vector<double> f2(7, 0);
	vector<double> rn(3, 0);
	vector<double> rm(3, 0);
	double temp_double;
	for (int m = 0; m < 7; m++) {
		rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
		for (int n = 0 ; n < 7; n++) {
			rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
			temp_double = distanceTwoPoint(rm, rn);
			//if (temp_double < 1) {
			//	f1[n] = 0;
			//	for (int i = 0; i < 111; i++) f1[n] += pow(-1, i) * pow(-1 + temp_double, i);
			//}
			//else
			f1[n] = 1 / temp_double;
		}
		f2[m] = dGaussianQuadTriangle(f1, triangleArea(Tn[0], Tn[1], Tn[2]));
	}
	return dGaussianQuadTriangle(f2, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
}
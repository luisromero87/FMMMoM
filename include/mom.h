#include <vector>
#include <complex>
#include <armadillo>

using namespace arma;

#ifndef MOM
#define MOM

class Mesh {
public:
	unsigned int Nnodes = 0;
	unsigned int Ntriangles = 0;
	unsigned int Nrwgbf = 0;
	vector<vector<double> > node;
	vector<vector<unsigned int> > triangle;
	vector<vector<unsigned int> > rwgbf; // v1, v2, vp, vn
	vector<vector<double> > rwgbfgeo;
	vector<double> Ein;
	cx_mat zmn;
	cx_vec bm;
	//Mesh(void);
	//~Mesh(void);
	int ReadMeshFile(string mymeshfile);
	void FindRwgbf(void);
	void BuildZmnMat(void);
	void BuildBmVec(void);
	complex<double> ComputeIntA0(vector<vector<double>> Tm, vector<vector<double>> Tn);
	complex<double> ComputeIntB0(vector<vector<double>> Tm, vector<vector<double>> Tn);
	complex<double> ComputeIntA12(vector<vector<double>> Tm, vector<vector<double>> Tn);
	complex<double> ComputeIntB12(vector<vector<double>> Tm, vector<vector<double>> Tn, int vm, int vn);
	complex<double> ComputeIntA3(vector<vector<double>> Tm);
	complex<double> ComputeIntB3(vector<vector<double>> Tm, int vm, int vn);
	void SetLambda (double lambda);
	void SetFreq (double freq);
	void SetEin (double Ex, double Ey, double Ez);

private:
	double knumber = 0;
	double omega = 0;
};





#endif

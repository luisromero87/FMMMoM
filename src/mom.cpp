#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <complex>
#include <armadillo>

using namespace arma;
using namespace std;
#include "mom.h"
#include "my_math.h"

/*-------------------------------------------------------------------------*//**
 * @brief      Reads a mesh file.
 *
 * @param[in]  mymeshfile  The mymeshfile
 *
 * @return     The error code
 * - 0: OK
 */
int Mesh::ReadMeshFile(string mymeshfile) {
  int error = 0;
  ifstream meshfile (mymeshfile.c_str());
  size_t offset;
  string line;
  string search = "Vertices"; // test variable to search in file
  if (meshfile.is_open()) {
    while (!meshfile.eof()) {
      getline(meshfile, line);
      if ((offset = line.find(search, 0)) != string::npos) {
        getline(meshfile, line);
        Nnodes = stoi(line);
        node.resize(Nnodes, vector<double> (3)) ;
        for (unsigned int i = 0; i < Nnodes; i++) {
          getline(meshfile, line);
          std::istringstream in(line);      //make a stream for the line itself
          for (int j = 0; j < 3; j++)
            in >> node[i][j];
        }
      }
    }
    meshfile.clear();
    meshfile.seekg(0);
    search = "Triangles"; // test variable to search in file
    while (!meshfile.eof()) {
      getline(meshfile, line);
      if ((offset = line.find(search, 0)) != string::npos) {
        getline(meshfile, line);
        Ntriangles = stoi(line);
        triangle.resize(Ntriangles, vector<unsigned int> (3)) ;
        for (unsigned int i = 0; i < Ntriangles; i++) {
          getline(meshfile, line);
          std::istringstream in(line);      //make a stream for the line itself
          for (int j = 0; j < 3; j++) {
            in >> triangle[i][j];
            triangle[i][j]--;
          }
          sort(triangle[i].begin(), triangle[i].end());
        }
      }
    }
    meshfile.close();
  }
  else cout << "Unable to open file." << endl;
  return 0;
}

/**
 * @brief      Find the RWG basis functions in the mesh.
 */
void Mesh::FindRwgbf(void) {
  vector<unsigned int>::iterator it;
  rwgbf.resize(0, vector<unsigned int> (6)) ;
  for (unsigned int i = 0; i < Ntriangles - 1; i++) {
    for (unsigned int j = i + 1; j < Ntriangles; j++) {
      for (unsigned int m = 0; m < 3 - 1; m++)
        for (unsigned int n = m + 1; n < 3; n++) {
          it = find(triangle[j].begin(), triangle[j].end(), triangle[i][m]);
          if ( it != triangle[j].end()) {
            it = find(triangle[j].begin(), triangle[j].end(), triangle[i][n]);
            if ( it != triangle[j].end()) {
              rwgbf.resize(rwgbf.size() + 1, vector<unsigned int> {triangle[i][m], triangle[i][n], triangle[i][3 - m - n], 0, 0, 0});
              for (unsigned int o = 0; o < 3; o++) {
                if (triangle[j][o] != triangle[i][m] && triangle[j][o] != triangle[i][n]) {
                  rwgbf[rwgbf.size() - 1][3] = triangle[j][o];
                }
              }
              rwgbf[rwgbf.size() - 1][4] = m;
              rwgbf[rwgbf.size() - 1][5] = n;
            }
          }
        }
    }
  }
  Nrwgbf = rwgbf.size();
  rwgbfgeo.resize(Nrwgbf, vector<double> (3));
  for (unsigned int i = 0; i < Nrwgbf; i++) {
    double sum = 0;
    for (unsigned int j = 0; j < 3; j++) sum += pow(node[rwgbf[i][0]][j] - node[rwgbf[i][1]][j], 2);
    sum = sqrt(sum);
    rwgbfgeo[i][0] = sum;
    rwgbfgeo[i][1] = triangleArea(node[rwgbf[i][0]], node[rwgbf[i][1]], node[rwgbf[i][2]]);
    rwgbfgeo[i][2] = triangleArea(node[rwgbf[i][0]], node[rwgbf[i][1]], node[rwgbf[i][3]]);
  }
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the type A integral
 *
 * @param[in]  Tm    The observation triangle
 * @param[in]  Tn    The suorce triangle
 *
 * @return     The value of the integral
 */
complex<double> Mesh::ComputeIntA0(vector<vector<double>> Tm, vector<vector<double>> Tn) {
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
      fn[n] = greenFunc (rm, rn, knumber);
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tn[0], Tn[1], Tn[2]));
  }
  return cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the type B integral
 *
 * @param[in]  Tm    The observation triangle
 * @param[in]  Tn    The source triangle
 *
 * @return     The int b 0.
 */
complex<double> Mesh::ComputeIntB0(vector<vector<double>> Tm, vector<vector<double>> Tn) {
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
      fn[n] = 0;
      for (int i = 0; i < 3; i++) fn[n] += (rm[i] - Tm[2][i]) * (rn[i] - Tn[2][i]);
      fn[n] *= greenFunc (rm, rn, knumber);
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tn[0], Tn[1], Tn[2]));
  }
  return cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the int a 12.
 *
 * @param[in]  Tm    The observation triangle
 * @param[in]  Tn    The source triangle
 *
 * @return     The int a 12.
 */
complex<double> Mesh::ComputeIntA12(vector<vector<double>> Tm, vector<vector<double>> Tn) {
  complex<double> inta = 0;
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  double R = 0;
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
      fn[n] = greenFuncExtSing (rm, rn, knumber);
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tn[0], Tn[1], Tn[2]));
  }
  inta = cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
  inta += ComputeIntegralType1GqtAndWilton(Tm, Tn);
  return inta;
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the int b 12.
 *
 * @param[in]  Tm    The observation triangle
 * @param[in]  Tn    The source triangle
 * @param[in]  vm    The opposite vertex
 * @param[in]  vn    The opposite vertex
 *
 * @return     The int b 12.
 */
complex<double> Mesh::ComputeIntB12(vector<vector<double>> Tm, vector<vector<double>> Tn, int vm, int vn) {
  complex<double> intb = 0;
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  double R = 0;
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tn[0], Tn[1], Tn[2]);
      R = distanceTwoPoint(rm, rn);
      fn[n] = 0;
      for (int i = 0; i < 3; i++) fn[n] += (rm[i] - Tm[vm][i]) * (rn[i] - Tn[vn][i]);
      fn[n] *= greenFuncExtSing (rm, rn, knumber);
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tn[0], Tn[1], Tn[2]));
  }
  intb = cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
  intb += ComputeIntegralType2GqtAndWilton(Tm, Tn, Tm[vm], Tn[vn]) ;
  return intb;
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the int a 3.
 *
 * @param[in]  Tm    The overlapping triangle
 *
 * @return     The int a 3.
 */
complex<double> Mesh::ComputeIntA3(vector<vector<double>> Tm) {
  complex<double> inta = 0;
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  double R = 0;
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tm[0], Tm[1], Tm[2]);
      R = distanceTwoPoint(rm, rn);
      fn[n] = Cplx(0, -1) * knumber - 0.5 * knumber * knumber * R + Cplx(0, 1) * knumber * knumber * knumber * R * R / 6.0 ;
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tm[0], Tm[1], Tm[2]));
  }
  inta = cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
  //cout << inta << endl;
  //cout << ComputeIntegralType1Arcioni (Tm) << endl;
  inta += ComputeIntegralType1Arcioni (Tm);
  return inta;
}

/*-------------------------------------------------------------------------*//**
 * @brief      Calculates the int b 3.
 *
 * @param[in]  Tm    The overlapping triangle
 * @param[in]  vm    The opposite vertex
 * @param[in]  vn    The opposite vertex
 *
 * @return     The int b 3.
 */
complex<double> Mesh::ComputeIntB3(vector<vector<double>> Tm, int vm, int vn) {
  complex<double> intb = 0;
  vector<complex<double>> fn(7, 0);
  vector<complex<double>> fm(7, 0);
  vector<double> rm (3, 0);
  vector<double> rn (3, 0);
  double R = 0;
  for (int m = 0; m < 7; m++) {
    rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
    for (int n = 0; n < 7; n++) {
      rn = xitor(M_XI1[n], M_XI2[n], M_XI3[n], Tm[0], Tm[1], Tm[2]);
      R = distanceTwoPoint(rm, rn);
      fn[n] = 0;
      for (int i = 0; i < 3; i++) fn[n] += (rm[i] - Tm[vm][i]) * (rn[i] - Tm[vn][i]);
      fn[n] *= (Cplx(0, -1) * knumber - 0.5 * knumber * knumber * R + Cplx(0, 1) * knumber * knumber * knumber * R * R / 6.0 );
    }
    fm[m] = cGaussianQuadTriangle(fn, triangleArea(Tm[0], Tm[1], Tm[2]));
  }
  intb = cGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2])) / (4.*M_PI);
  //cout << intb << endl;
  //cout << ComputeIntegralType2EibertHansen (Tm, Tm[vm], Tm[vn]) << endl;
  intb += ComputeIntegralType2EibertHansen (Tm, Tm[vm], Tm[vn]);
  return intb;
}

/**
 * @brief      Builds the impedance Zmn matrix.
 */
void Mesh::BuildZmnMat(void) {
  vector<vector<double>> Tm(3, vector<double>(3, 0));
  vector<vector<double>> Tn(3, vector<double>(3, 0));
  vector<unsigned int> iTm(3, 0);
  vector<unsigned int> iTn(3, 0);
  vector<complex<double>> inta (4, 0);
  vector<complex<double>> intb (4, 0);
  complex<double> cinta = 0;
  complex<double> cintb = 0;
  int commonnodes = 0;
  //bm.set_size(Nrwgbf);
  //bm.fill(0);
  zmn.set_size(Nrwgbf, Nrwgbf);
  zmn.fill(0);
  //for (unsigned int i = 0; i < Nrwgbf; i++) zmn[i].resize(Nrwgbf, 0);

  for (unsigned int m = 0; m < Nrwgbf; m++) {
    iTm[0] = rwgbf[m][0];
    iTm[1] = rwgbf[m][1];
    for (unsigned int n = m; n < Nrwgbf; n++) {
      iTn[0] = rwgbf[n][0];
      iTn[1] = rwgbf[n][1];
      for (int i = 0; i < 4; i++) {
        inta[i] = 0.0;
        intb[i] = 0.0;
      }
      for (int tm = 0; tm < 2; tm++) {
        iTm[2] = rwgbf[m][2 + tm];
        for (int tn = 0; tn < 2; tn++) {
          iTn[2] = rwgbf[n][2 + tn];
          commonnodes = 0;
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              if (iTm[i] == iTn[j]) {
                commonnodes++;
                break;
              }
            }
          }
          Tm[0] = node[iTm[0]];
          Tm[1] = node[iTm[1]];
          Tm[2] = node[iTm[2]];
          Tn[0] = node[iTn[0]];
          Tn[1] = node[iTn[1]];
          Tn[2] = node[iTn[2]];
          //if (commonnodes != 3) commonnodes = -1;
          if (commonnodes == 0) {
            vector<double> cm (3, 0);
            vector<double> cn (3, 0);
            for (int i = 0; i < 3; i++) {
              cm[i] += (Tm[0][i] + Tm[1][i] + Tm[2][i]) / 3.0;
              cn[i] += (Tn[0][i] + Tn[1][i] + Tn[2][i]) / 3.0;
            }
            if (distanceTwoPoint (cm, cn) <= 0.1 * 2 * M_PI / knumber) {
              commonnodes = 0;
            }
            else commonnodes = -1;
          }
          switch (commonnodes) {
          case 3:
            int vn;
            for (int i = 0; i < 3; i++) {
              if (iTn[2] == iTm[i]) {
                vn = i;
                break;
              }
            }
            //cout << "Overlaping triangle: " << vn << endl;
            //cout << iTm[0] << "\t" << iTm[1] << "\t" << iTm[2] << endl << iTn[0] << "\t" << iTn[1] << "\t" << iTn[2] << endl;
            inta [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntA3(Tm);
            intb [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntB3(Tm, 2, vn);
            break;
          case 2:
          case 1:
          case 0:
            //case 0:
            //cout << "Near triangles: " << commonnodes << endl;
            //cout << iTm[0] << "\t" << iTm[1] << "\t" << iTm[2] << endl << iTn[0] << "\t" << iTn[1] << "\t" << iTn[2] << endl;
            inta [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntA12(Tm, Tn);
            intb [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntB12(Tm, Tn, 2, 2);
            break;
          case -1:
            //cout << "Far triangles: " << commonnodes << endl;
            //cout << iTm[0] << "\t" << iTm[1] << "\t" << iTm[2] << endl << iTn[0] << "\t" << iTn[1] << "\t" << iTn[2] << endl;
            inta [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntA0(Tm, Tn);
            intb [2 * tm + tn] = 1.0 / (rwgbfgeo[m][1 + tm] * rwgbfgeo[n][1 + tn]) * pow(-1, tm + tn) * ComputeIntB0(Tm, Tn);
            break;
          }
        }
      }
      //cout << "IntA: ";
      //for (int i = 0; i < 4; i++) cout << inta[i] << "\t";
      //cout << endl;
      //cout << "IntB: ";
      //for (int i = 0; i < 4; i++) cout << intb[i] << "\t";
      //cout << endl;
      cinta = cintb = 0.0;
      for (int i = 0; i < 4; i++) {
        cinta += inta[i];
        cintb += intb[i];
      }
      zmn(m, n) = complex<double>(0, 1) * omega * M_MU0 / 4.0 * cintb +
                  1.0 / (complex<double>(0, 1) * omega * M_EPS0) * cinta;
      zmn(m, n) *= rwgbfgeo[m][0] * rwgbfgeo[n][0];
      zmn(n, m) = zmn(m, n);
    }
  }
}

/**
 * @brief      Builds the bm vector.
 */
void Mesh::BuildBmVec(void) {
  bm.set_size(Nrwgbf);
  bm.fill(0);
  vector<double> fm(7, 0);
  vector<double> rm (3, 0);
  vector<vector<double>> Tm(3, vector<double>(3, 0));
  double rhodote = 0;
  for (unsigned int mm = 0; mm < Nrwgbf; mm++) {
    Tm[0] = node[rwgbf[mm][0]];
    Tm[1] = node[rwgbf[mm][1]];
    Tm[2] = node[rwgbf[mm][2]];
    for (int m = 0; m < 7; m++) {
      rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
      rhodote = 0;
      for (int i = 0; i < 3; i++) rhodote += ((rm[i] - Tm[2][i]) * Ein[i]);
      fm[m] = rhodote;
    }
    bm(mm) = rwgbfgeo[mm][0] / (2.* rwgbfgeo[mm][1]) * dGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
    Tm[0] = node[rwgbf[mm][0]];
    Tm[1] = node[rwgbf[mm][1]];
    Tm[2] = node[rwgbf[mm][3]];
    for (int m = 0; m < 7; m++) {
      rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], Tm[0], Tm[1], Tm[2]);
      rhodote = 0;
      for (int i = 0; i < 3; i++) rhodote += ((Tm[2][i] - rm[i] ) * Ein[i]);
      fm[m] = rhodote;
    }
    bm(mm) = bm(mm) + rwgbfgeo[mm][0] / (2.* rwgbfgeo[mm][2]) * dGaussianQuadTriangle(fm, triangleArea(Tm[0], Tm[1], Tm[2]));
    //bm(mm) = bm(mm) * complex<double>(0, -1) / (omega * M_MU0);
  }
}

/*-------------------------------------------------------------------------*//**
 * @brief      Sets the working lambda (wavelength).
 *
 * @param[in]  lambda  The wavelength
 */
void Mesh::SetLambda (double lambda) {
  knumber = 2 * M_PI / lambda;
  omega = knumber * M_C;
}

/*-------------------------------------------------------------------------*//**
 * @brief      Sets the working frequency.
 *
 * @param[in]  freq  The frequency
 */
void Mesh::SetFreq (double freq) {
  omega = 2 * M_PI * freq;
  knumber = omega / M_C;
}

/*-------------------------------------------------------------------------*//**
 * @brief      Sets the inciden electric field.
 *
 * @param[in]  Ex    Electric field in x direction
 * @param[in]  Ey    Electric field in y direction
 * @param[in]  Ez    Electric field in z direction
 */
void Mesh::SetEin (double Ex, double Ey, double Ez) {
  Ein.resize(3);
  Ein[0] = Ex;
  Ein[1] = Ey;
  Ein[2] = Ez;
}
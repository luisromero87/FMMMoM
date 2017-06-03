#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkMath.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vector>

#include <iostream>
#include <cmath>
#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

#include "my_math.h"
#include "mom.h"

int main(int , char *[])
{
  std::string outputFilename = "output.vtp";
  Mesh MyMesh;

  vec freq = linspace<vec>(0.01e9, 1e9, 20);
  cx_vec zin(20);
  freq.print("freq:");

  double Length = 1.0;
  double width = 0.003;

  int NN = 18;

  MyMesh.Nnodes = 3 * NN + 2;
  MyMesh.node.resize(MyMesh.Nnodes, vector<double> (3));

  for (int i = 0; i < NN; i++) {
    MyMesh.node[2 * i] = { -Length / 2 + Length / NN * i, 0, 0};
    MyMesh.node[2 * i + 1] = { -Length / 2 + Length / NN * i, width, 0};
    MyMesh.node[2 * NN + 2 + i] = { -Length / 2 + Length * (2.0 * i + 1.0) / (2.0 * NN), width / 2, 0};
  }
  MyMesh.node[2 * NN] = { Length / 2, 0, 0};
  MyMesh.node[2 * NN + 1] = { Length / 2, width, 0};


  MyMesh.Ntriangles = 4 * NN;
  MyMesh.triangle.resize(MyMesh.Ntriangles, vector<unsigned int> (3)) ;
  for (unsigned int i = 0; i < NN; i++) {
    MyMesh.triangle[4 * i + 0] = {2 * i + 0, 2 * NN + i + 2, 2 * i + 1};
    MyMesh.triangle[4 * i + 1] = {2 * i + 1, 2 * NN + i + 2, 2 * i + 3};
    MyMesh.triangle[4 * i + 2] = {2 * i + 3, 2 * NN + i + 2, 2 * i + 2};
    MyMesh.triangle[4 * i + 3] = {2 * i + 2, 2 * NN + i + 2, 2 * i + 0};
  }

  MyMesh.FindRwgbf();
  cout << "Nrwgbf :" << MyMesh.Nrwgbf << endl;

  int input;
  for (unsigned int i = 0; i < MyMesh.Nrwgbf; i++) {
    if ((MyMesh.rwgbf[i][0] == NN + 1 && MyMesh.rwgbf[i][1] == NN + 2) || (MyMesh.rwgbf[i][0] == NN + 2 && MyMesh.rwgbf[i][1] == NN + 1)) {
      //MyMesh.bm(i) = -complex<double>(0, 1) * MyMesh.rwgbfgeo[i][0] / (MyMesh.omega * M_MU0);
      input = i;
      cout << "input -> bm(" << input << ")" << endl;
    }
  }

  // Setup 3 points
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < MyMesh.Nnodes; i++) {
    points->InsertNextPoint(MyMesh.node[i][0], MyMesh.node[i][1], MyMesh.node[i][2]);
  }

  // Create a triangle
  vtkSmartPointer<vtkTriangle> triangle =
    vtkSmartPointer<vtkTriangle>::New();
  // Add the triangle to a cell array
  vtkSmartPointer<vtkCellArray> triangles =
    vtkSmartPointer<vtkCellArray>::New();
  for (unsigned int i = 0; i < MyMesh.Ntriangles; i++) {
    triangle->GetPointIds()->SetId(0, MyMesh.triangle[i][0]);
    triangle->GetPointIds()->SetId(1, MyMesh.triangle[i][1]);
    triangle->GetPointIds()->SetId(2, MyMesh.triangle[i][2]);
    triangles->InsertNextCell(triangle);
  }

  MyMesh.SetLambda(2 * Length);
  MyMesh.SetEin(1, 0, 0);

  MyMesh.BuildZmnMat();
  MyMesh.BuildBmVec();
  cx_vec a = solve(MyMesh.zmn, MyMesh.bm);
  vec aa(NN - 1);
  for (int j = 0; j < NN - 1; j++)aa(j) = abs(a(5 * j + 4));
  string filename = "am_abs";
  aa.save(filename, raw_ascii);

  vector<vector<complex<double>>> Jm(MyMesh.Ntriangles, vector<complex<double>>(3, 0));
  for (unsigned int j = 0; j < MyMesh.Nrwgbf; j++) {
    vector<double> rm (3, 0);
    vector<double> f (7, 0);
    for (int i = 0; i < 3; i++) {
      for (int m = 0; m < 7; m++) {
        rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][2]]);
        f[m] = rm[i] - MyMesh.node[MyMesh.rwgbf[j][2]][i];
      }
      Jm[MyMesh.rwgbf[j][4]][i] += MyMesh.bm[j] * dGaussianQuadTriangle(f, triangleArea(MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][2]]));
      for (int m = 0; m < 7; m++) {
        rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][3]]);
        f[m] = MyMesh.node[MyMesh.rwgbf[j][3]][i] - rm[i];
      }
      Jm[MyMesh.rwgbf[j][5]][i] += MyMesh.bm[j] * dGaussianQuadTriangle(f, triangleArea(MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][3]]));
    }
  }


  // Setup data for the triangle. Attach a value of 1.45.
  // This can be anything you wish to store with it)
  vtkSmartPointer<vtkDoubleArray> triangleData =
    vtkSmartPointer<vtkDoubleArray>::New();
  triangleData->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
  triangleData->SetName("Js_mag"); //set the name of the value
    double Jmabs=0;
  for (int i = 0; i < MyMesh.Ntriangles; i++) {
    //triangleData->InsertNextValue(1.45 / MyMesh.Ntriangles * i); //set the actual value
    Jmabs = 0;
    for (int j=0;j<3;j++){
      Jmabs+=abs(Jm[i][j])*abs(Jm[i][j]);
    }
    cout<<Jmabs<<endl;
    triangleData->InsertNextValue(sqrt(Jmabs)); //set the actual value
  }

  // Create a polydata that contains the points,
  // the triangle on those points, and the data
  // array (value) we created for the triangle
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetPolys(triangles);
  polydata->GetCellData()->AddArray(triangleData);
  //polydata->GetCellData()->SetScalars(triangleData);
  polydata->GetCellData()->SetActiveScalars("Js_mag");

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif
  writer->SetFileName(outputFilename.c_str());
  writer->Write();


  // Setup the visualization pipeline
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polydata);
#else
  mapper->SetInputData(polydata);
#endif

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);
  renderer->SetBackground(.1, .2, .3);
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

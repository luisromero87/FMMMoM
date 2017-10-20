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
#include <vtkLegendScaleActor.h>
#include <vtkScalarBarActor.h>

#include <vector>

#include <iostream>
#include <cmath>
#include <string>
#include <armadillo>

#include <cstdlib>

#include <libconfig.h++>

using namespace libconfig;
using namespace arma;
using namespace std;

#include "my_math.h"
#include "mom.h"

int main(int , char *[])
{

  Config cfg;
  try
  {
    cfg.readFile("example.cfg");
  }
  catch (const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return (EXIT_FAILURE);
  }
  catch (const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return (EXIT_FAILURE);
  }
  // Get the settings.
  try
  {
    if (cfg.exists("settings")) {
      cout << "Reading the settings..." << endl << endl;

    }
  }
  catch (const SettingNotFoundException &nfex)
  {
    cerr << "No settings found in configuration file." << endl;
  }

  const Setting& settings = cfg.getRoot();

  string inputmesh;
  try
  {
    inputmesh = cfg.lookup("settings.mesh").c_str();
    cout << "Mesh: " << inputmesh << endl << endl;
  }
  catch (const SettingNotFoundException &nfex)
  {
    cerr << "No 'mesh' setting in configuration file." << endl;
  }


  string outputFolder;
  try
  {
    outputFolder = cfg.lookup("settings.output_folder").c_str();
    cout << "output_folder: " << outputFolder << endl << endl;
  }
  catch (const SettingNotFoundException &nfex)
  {
    cerr << "No 'output' setting in configuration file." << endl;
  }

 
const int dir_err = system((string("mkdir -p ")+outputFolder).c_str());
if (-1 == dir_err)
{
    printf("Error creating directory!n");
    exit(1);
}

  vector<double> wavelength;
  if (cfg.exists("settings.wavelength")) {
    if (cfg.lookup("settings.wavelength.sweep")) {
      cout << "wavelength sweep" << endl;
    }
    else {
      int count = cfg.lookup("settings.wavelength.value").getLength();
      wavelength.resize(count, 0);
      double* pwavelength = (wavelength.data());
      for (int i = 0; i < count; i++) pwavelength[i] = cfg.lookup("settings.wavelength.value")[i];
      cout << "wavelength value: " ;
      for (int i = 0; i < count; i++) cout << wavelength[i] << endl;
    }
  }

  vector<double> angle(3, 0);
  if (cfg.exists("settings.source")) {
    angle[0] = cfg.lookup("settings.source.yaw");
    angle[1] = cfg.lookup("settings.source.pitch");
    angle[2] = cfg.lookup("settings.source.roll");
    cout << "yaw, pitch roll: " ;
    for (int i = 0; i < 3; i++) cout << angle[i] << "\t";
    string units = cfg.lookup("settings.source.units").c_str();
    if (units.compare("grades") == 0) {
      cout << "grades" << endl;
      for (int i = 0; i < 3; i++) angle[i] = angle[i] / 180.0 * M_PI;
    }
    cout << endl;
  }

  Mesh MyMesh;

  vec freq = linspace<vec>(0.01e9, 1e9, 20);
  cx_vec zin(20);
  freq.print("freq:");

  MyMesh.ReadMeshFile(inputmesh);
  MyMesh.FindRwgbf();

  cout << "Nrwgbf :" << MyMesh.Nrwgbf << endl;

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

  MyMesh.SetLambda(wavelength[0]);
  vector<double> Ein ({cos(angle[0])*sin(angle[1])*sin(angle[2]) - cos(angle[2])*sin(angle[0]),
                       cos(angle[0])*cos(angle[2]) + sin(angle[0])*sin(angle[1])*sin(angle[2]),
                       cos(angle[1])*sin(angle[2])
                      });
  cout<<"Ein: "<<Ein[0]<<"\t"<<Ein[1]<<"\t"<<Ein[2]<<endl;
  MyMesh.SetEin(Ein[0], Ein[1], Ein[2]);

  MyMesh.BuildZmnMat();
  MyMesh.BuildBmVec();
  cx_vec a = solve(MyMesh.zmn, MyMesh.bm);
  // vec aa(NN - 1);
  // for (int j = 0; j < NN - 1; j++)aa(j) = abs(a(5 * j + 4));
  // string filename = "am_abs";
  // aa.save(filename, raw_ascii);

  vector<vector<complex<double>>> Jm(MyMesh.Ntriangles, vector<complex<double>>(3, 0));
  vector<double> rm (3, 0);
  vector<double> f (7, 0);
  for (unsigned int j = 0; j < MyMesh.Nrwgbf; j++) {
    cout << "rwgbf: " << j << endl;
    cout << "T+: " << MyMesh.rwgbf[j][4] << "\tT-:" << MyMesh.rwgbf[j][5] << endl;
    for (int i = 0; i < 3; i++) {
      for (int m = 0; m < 7; m++) {
        rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][2]]);
        f[m] = rm[i] - MyMesh.node[MyMesh.rwgbf[j][2]][i];
      }
      Jm[MyMesh.rwgbf[j][4]][i] += a[j] * MyMesh.rwgbfgeo[j][0] / MyMesh.rwgbfgeo[j][1] * dGaussianQuadTriangle(f, triangleArea(MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][2]]));
      for (int m = 0; m < 7; m++) {
        rm = xitor(M_XI1[m], M_XI2[m], M_XI3[m], MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][3]]);
        f[m] = MyMesh.node[MyMesh.rwgbf[j][3]][i] - rm[i];
      }
      Jm[MyMesh.rwgbf[j][5]][i] += a[j] * MyMesh.rwgbfgeo[j][0] / MyMesh.rwgbfgeo[j][2] * dGaussianQuadTriangle(f, triangleArea(MyMesh.node[MyMesh.rwgbf[j][0]], MyMesh.node[MyMesh.rwgbf[j][1]], MyMesh.node[MyMesh.rwgbf[j][3]]));
    }
  }


  // Setup data for the triangle. Attach a value of 1.45.
  // This can be anything you wish to store with it)
  vtkSmartPointer<vtkDoubleArray> triangleData =
    vtkSmartPointer<vtkDoubleArray>::New();
  triangleData->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
  triangleData->SetName("Js_mag"); //set the name of the value
  double Jmabs = 0;
  for (int i = 0; i < MyMesh.Ntriangles; i++) {
    //triangleData->InsertNextValue(1.45 / MyMesh.Ntriangles * i); //set the actual value
    for (int j = 0; j < 3; j++) {
      Jm[i][j] = Jm[i][j] / triangleArea(MyMesh.node[MyMesh.triangle[i][0]], MyMesh.node[MyMesh.triangle[i][1]], MyMesh.node[MyMesh.triangle[i][2]]);
    }
    Jmabs = 0;
    for (int j = 0; j < 3; j++) {
      Jmabs += abs(Jm[i][j]) * abs(Jm[i][j]);
    }
    //Jmabs/=triangleArea(MyMesh.node[MyMesh.triangle[0]],MyMesh.node[MyMesh.triangle[1]],MyMesh.node[MyMesh.triangle[2]]);
    cout << i << ": " << sqrt(Jmabs) << endl;
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
  writer->SetFileName((outputFolder+"/out.vtp").c_str());
  writer->Write();


  // Setup the visualization pipeline
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polydata);
#else
  mapper->SetInputData(polydata);
#endif
  mapper->SetScalarRange(0, 6);

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
  // vtkSmartPointer<vtkLegendScaleActor> legendScaleActor =
  //  vtkSmartPointer<vtkLegendScaleActor>::New();
  //renderer->AddActor(legendScaleActor);
  vtkSmartPointer<vtkScalarBarActor> scalarBar =
    vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(mapper->GetLookupTable());
  scalarBar->SetTitle("Title");
  scalarBar->SetNumberOfLabels(4);
  scalarBar->SetMaximumWidthInPixels(25);

  renderer->AddActor(scalarBar);
  renderer->SetBackground(.1, .2, .3);
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

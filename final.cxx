/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>

#include "astro.cxx"


// Change to false for no printing
const bool P = true;


const int x = 0;
const int y = 1;
const int z = 2;


double *VectorDotProduct( const double *A, const double *B ) {
    static double ret[3] = {    
        A[x] * B[x], 
        A[y] * B[y],
        A[z] * B[z] 
    };
    return ret;
}

double *VectorCrossProduct( const double *A, const double *B ) {
    static double ret[3] = {
        (( A[y] * B[z] ) - ( A[z] * B[y] )),
        (( A[z] * B[x] ) - ( A[x] * B[z] )),
        (( A[x] * B[y] ) - ( A[y] * B[x] ))
    };
    return ret;
}

double VectorMagnitude( const double *A ) {
    return sqrt( pow(A[x], 2) + pow(A[y], 2) + pow(A[z], 2) );
}




int main()
{

    Camera camera = SetupCamera();
    
    int WIDTH = 1024;
    int HEIGHT = 1024;

    double FOV_X = camera.angle;
    double FOV_Y = camera.angle;




    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();
    cerr << "After update, file has " << rdr->GetOutput()->GetNumberOfCells() << " cells." << endl;

    // Look = Distance from focus to camera
    double look[ 3 ] = {
        camera.focus[x] - camera.position[x],
        camera.focus[y] - camera.position[y],
        camera.focus[z] - camera.position[z]
    };  if(P)printf( "look:\t\tx: %f, y: %f, z: %f\n", look[x], look[y], look[z] );
        if(P)printf( "up:  \t\tx: %f, y: %f, z: %f\n", camera.up[x], camera.up[y], camera.up[z] );

    double look_length = VectorMagnitude( look );
        if(P)printf( "look_length:\t%f\n", look_length );

    double *lookXUp = VectorCrossProduct( look, camera.up );
        if(P)printf( "lookXUp:\tx: %f, y: %f, z: %f\n", lookXUp[x], lookXUp[y], lookXUp[z] );

    double ru[ 3 ] = { 
        lookXUp[x] / VectorMagnitude( lookXUp ),
        lookXUp[y] / VectorMagnitude( lookXUp ),
        lookXUp[z] / VectorMagnitude( lookXUp )
    };  if(P)printf( "ru:\t\tx: %f, y: %f, z: %f\n", ru[x], ru[y], ru[z] );

    double rv[ 3 ] = {
    };  if(P)printf( "rv:\t\tx: %f, y: %f, z: %f\n", rv[x], rv[y], rv[z] );

    // Calculate field of view for x and y
    // Inverse Tangent // tan^-1( W / look )

    if(P)printf( "Look length:\t%f\n", look_length);

    
    //double delta_rx[ 3 ] = {
        //  (2 * tan(/2))/(







  vtkContourFilter *cf = vtkContourFilter::New();
  cf->SetNumberOfContours(1);
  cf->SetValue(0, 4.0);
  cf->SetInputConnection(rdr->GetOutputPort());
  cf->Update();

  // The mapper is responsible for pushing the geometry into the graphics
  // library. It may also do color mapping, if scalars or other attributes
  // are defined. 
  //
  vtkSmartPointer<vtkPolyDataMapper> win1Mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  win1Mapper->SetInputConnection(cf->GetOutputPort());

  vtkSmartPointer<vtkActor> win1Actor =
    vtkSmartPointer<vtkActor>::New();
  win1Actor->SetMapper(win1Mapper);

  vtkSmartPointer<vtkRenderer> ren1 =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  ren1->AddActor(win1Actor);
  ren1->SetBackground(0.0, 0.0, 0.0);
  renWin->SetSize(WIDTH, HEIGHT);

   
  ren1->GetActiveCamera()->SetFocalPoint(camera.focus[0],camera.focus[1],camera.focus[2]);
  ren1->GetActiveCamera()->SetPosition(camera.position[0],camera.position[1],camera.position[2]);
  ren1->GetActiveCamera()->SetViewUp(camera.up[0],camera.up[1],camera.up[2]);
  ren1->GetActiveCamera()->SetClippingRange(camera.near, camera.far);
  ren1->GetActiveCamera()->SetViewAngle(camera.angle);
  ren1->GetActiveCamera()->SetDistance(70);
  
  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  iren->Initialize();
  iren->Start();

  return EXIT_SUCCESS;
}





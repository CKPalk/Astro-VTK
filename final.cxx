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


const bool T = true; // Change to false for no buffer


const int x = 0;
const int y = 1;
const int z = 2;


double *VectorDotProduct( const double *A, const double *B, double ret[] ) {
    ret[x] = A[x] * B[x];
    ret[y] = A[y] * B[y];
    ret[z] = A[z] * B[z];
    return ret;
}

double *VectorCrossProduct( const double *A, const double *B, double ret[] ) {
    ret[x] = (( A[y] * B[z] ) - ( A[z] * B[y] ));
    ret[y] = (( A[z] * B[x] ) - ( A[x] * B[z] ));
    ret[z] = (( A[x] * B[y] ) - ( A[y] * B[x] ));
    return ret;
}

double VectorMagnitude( const double *A ) {
    return sqrt( pow(A[x], 2) + pow(A[y], 2) + pow(A[z], 2) );
}

double *rayForPixel( int i, int j, double *look, double *delta_rx, double *delta_ry, int width, int height, double ret[] ) {
    ret[x] = ( look[x] / VectorMagnitude( look ) ) + 
        (((( 2 * i ) + 1 - width  ) / 2 ) * delta_rx[x] ) + 
        (((( 2 * j ) + 1 - height ) / 2 ) * delta_ry[x] );
    ret[y] = ( look[y] / VectorMagnitude( look ) ) + 
        (((( 2 * i ) + 1 - width  ) / 2 ) * delta_rx[y] ) + 
        (((( 2 * j ) + 1 - height ) / 2 ) * delta_ry[y] );
    ret[z] = ( look[z] / VectorMagnitude( look ) ) + 
        (((( 2 * i ) + 1 - width  ) / 2 ) * delta_rx[z] ) + 
        (((( 2 * j ) + 1 - height ) / 2 ) * delta_ry[z] );
    return ret;
}







int main()
{

    Camera camera = SetupCamera();
    
    const int WIDTH     = 1024; // px
    const int HEIGHT    = 1024; // px;

    double FOV_X = camera.angle;
    double FOV_Y = camera.angle;


    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();
    cerr << "After update, file has " << rdr->GetOutput()->GetNumberOfCells() << " cells." << endl;

	// Gather data from file
	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
	rgrid->GetDimensions(dims);

	float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
	float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);

	float *F = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);





    // Look = Distance from focus to camera
    double look[ 3 ] = {
        camera.focus[x] - camera.position[x],
        camera.focus[y] - camera.position[y],
        camera.focus[z] - camera.position[z]
    };  if(T)printf( "look:\t\tx: %f, y: %f, z: %f\n", look[x], look[y], look[z] );
        if(T)printf( "up:  \t\tx: %f, y: %f, z: %f\n", camera.up[x], camera.up[y], camera.up[z] );

    // Magnitude of look
    double look_length = VectorMagnitude( look );
        if(T)printf( "look_length:\t%f\n", look_length );


    // vector cross product of look and up
    double lookXup_ret[3];
    double *lookXup = VectorCrossProduct( look, camera.up, lookXup_ret );
        if(T)printf( "lookXup:\tx: %f, y: %f, z: %f\n", lookXup[x], lookXup[y], lookXup[z] );

    double ru[ 3 ] = { 
        lookXup[x] / VectorMagnitude( lookXup ),
        lookXup[y] / VectorMagnitude( lookXup ),
        lookXup[z] / VectorMagnitude( lookXup )
    };  if(T)printf( "ru:\t\tx: %f, y: %f, z: %f\n", ru[x], ru[y], ru[z] );

    double lookXru_ret[3];
    double *lookXru = VectorCrossProduct( look, ru, lookXru_ret );
        if(T)printf( "lookXru:\tx: %f, y:%f, z: %f\n", lookXru[x], lookXru[y], lookXru[z] );

    double rv[ 3 ] = {
        lookXru[x] / VectorMagnitude( lookXru ),
        lookXru[y] / VectorMagnitude( lookXru ),
        lookXru[z] / VectorMagnitude( lookXru )
    };  if(T)printf( "rv:\t\tx: %f, y: %f, z: %f\n", rv[x], rv[y], rv[z] );

    double delta_rx[ 3 ] = {
        (( 2 * tan( FOV_X / 2 ) ) / WIDTH ) * ru[x],
        (( 2 * tan( FOV_X / 2 ) ) / WIDTH ) * ru[y],
        (( 2 * tan( FOV_X / 2 ) ) / WIDTH ) * ru[z]
    };  if(T)printf( "delta_rx\tx: %f, y: %f, z: %f\n", delta_rx[x], delta_rx[y], delta_rx[z] );

    double delta_ry[ 3 ] = {
        (( 2 * tan( FOV_Y / 2 ) ) / HEIGHT ) * rv[x],
        (( 2 * tan( FOV_Y / 2 ) ) / HEIGHT ) * rv[y],
        (( 2 * tan( FOV_Y / 2 ) ) / HEIGHT ) * rv[z]
    };  if(T)printf( "delta_ry\tx: %f, y: %f, z: %f\n", delta_ry[x], delta_ry[y], delta_ry[z] );


    // Iterate over pixels on screen
    int pixel_w;
    int pixel_h;
    for ( pixel_h = 0; pixel_h < HEIGHT; pixel_h++ ) {
        for ( pixel_w = 0; pixel_w < WIDTH; pixel_w++ ) {
            double ray_ret[3];
            double *ray = rayForPixel( pixel_w, pixel_h, look, delta_rx, delta_ry, WIDTH, HEIGHT, ray_ret );
                if(T)printf( "Ray for (%4d, %4d): x: %f, y: %f, z: %f\n", pixel_w, pixel_h, ray[x], ray[y], ray[z] );

			//



        }
    }






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

   
  ren1->GetActiveCamera()->SetFocalPoint(camera.focus[x], camera.focus[y], camera.focus[z]);
  ren1->GetActiveCamera()->SetPosition(camera.position[x], camera.position[y], camera.position[z]);
  ren1->GetActiveCamera()->SetViewUp(camera.up[x], camera.up[y], camera.up[z]);
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





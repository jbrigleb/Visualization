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
#include <vtkPNGWriter.h>

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
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


// ****************************************************************************
//  Function: IdentifyCase
//
//  Arguments:
//      cell: the logical index of the cell being considered
//      F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//      
//  Returns: an int value to be used as an index, case type
// ****************************************************************************

int IndentifyCase(int *cell, float iVal, const int *dims, const float *F)
{
    int v0 [2] = {};
    int v1 [2] = {};
    int v2 [2] = {};
    int v3 [2] = {};
    int * vl[4] = {v0,v1,v2,v3};
    
    
    v0[0] = cell[0];
    v0[1] = cell[1];
    
    v1[0] = cell[0] + 1;
    v1[1] = cell[1];
    
    v2[0] = cell[0];
    v2[1] = cell[1] + 1;
    
    v3[0] = cell[0] + 1;
    v3[1] = cell[1] + 1;
    
    int b0,b1,b2,b3;
    int bl[4];
    
    for(int i = 0;i < 4; i++){
        if(F[GetPointIndex(vl[i], dims)] < iVal){
            bl[i] = 0;
        }else{
            bl[i] = 1;
        }
    }
    
    return bl[3]*8 + bl[2]*4 + bl[1]*2 + bl[0]*1;
    
}
    

class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}




int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

	// ALGORITHM IMPLEMENTED BY JACK STARTS HERE
    
    //iso value
    float isoVal = 3.2;
    
    int numCells = GetNumberOfCells(dims);
    
    //lup stores the edges between which edges are drawn in the cases
    int lup [16][4] = {};
    
    
    //cases without edges, or extra spaces are set to -1
    for(int i = 0; i<16; i++){
        for(int j = 0; j < 4; j++){
            lup[i][j] = -1;
        }
    }
    
    
    // numSegments stores the number of segments in each case in lup
    int numSegments [16] = { };
    
    for(int i = 0; i<16;i++){
        numSegments[i] = 1;
    }
    
    numSegments[0] = 0;
    numSegments[15] = 0;
    numSegments[6] = 2;
    numSegments[9] = 2;
    
    //case 0 no segments, covered by initial loop
    
    //case 1 segment between 0 and 3 
    lup [1][0] = 0;
    lup [1][1] = 3;
    
    //case 2 segment between 0 and 1
    lup [2][0] = 0;
    lup [2][1] = 1;
    
    //case 3 segment between 1 and 3
    lup [3][0] = 1;
    lup [3][1] = 3;

    //case 4 segment between 2 and 3
    lup [4][0] = 2;
    lup [4][1] = 3;
    
    //case 5 segment between 0 and 2
    lup [5][0] = 0;
    lup [5][1] = 2;

    //case 6 segments between 0 and 1, 2 and 3
    lup [6][0] = 0;
    lup [6][1] = 1;
    lup [6][2] = 2;
    lup [6][3] = 3;
    
    //case 7 segment between 1 and 2
    lup [7][0] = 1;
    lup [7][1] = 2;
    
    //case 8 segment between 1 and 2
    lup [8][0] = 1;
    lup [8][1] = 2;
    
    //case 9 segments between 1 and 2, 0 and 3
    lup [9][0] = 1;
    lup [9][1] = 2;
    lup [9][2] = 0;
    lup [9][3] = 3;
    
    //case 10 segment between 0 and 2
    lup [10][0] = 0;
    lup [10][1] = 2;
    
    //case 11 segment between 2 and 3
    lup [11][0] = 2;
    lup [11][1] = 3;
    
    //case 12 segment between 1 and 3
    lup [12][0] = 1;
    lup [12][1] = 3;
    
    //case 13 segment between 0 and 1 
    lup [13][0] = 0;
    lup [13][1] = 1;
    
    //case 14 segment between 0 and 3
    lup [14][0] = 0;
    lup [14][1] = 3;
    
    //case 15 no segments, covered in initial loop
    
    
    for(int i = 0; i < numCells; i++){
        int ind[2] = {};
        GetLogicalCellIndex(ind, i, dims);
        int cellCase = IndentifyCase(ind, isoVal, dims, F);
        //cerr<<"....."<<cellCase<<"....\n";
        
        int nseg = numSegments[cellCase];
        
        for(int j = 0; j < nseg; j++){
            
            int edge1 = lup[cellCase][2*j];
            int edge2 = lup[cellCase][2*j+1];
            
            int edges[2] = {edge1, edge2};
            
            float point1[2] = {};
            float point2[2] = {};
            float *points[2] = {point1, point2};
            
            
            for(int k = 0; k < 2; k++){
            
                int v1[2] = {};
                int v2[2] = {};
                
                if(edges[k] == 0){
                    v1[0] = ind[0];
                    v1[1] = ind[1];
                    v2[0] = ind[0]+1;
                    v2[1] = ind[1];
                    
                    points[k][1] = Y[v1[1]];
                    points[k][0] = X[v1[0]] + ((isoVal - F[GetPointIndex(v1, dims)]) / (F[GetPointIndex(v2, dims)] - F[GetPointIndex(v1, dims)])) * (X[v2[0]] - X[v1[0]]);
                    
                }else if(edges[k] == 1){
                    v1[0] = ind[0]+1;
                    v1[1] = ind[1];
                    v2[0] = ind[0]+1;
                    v2[1] = ind[1]+1;
                    
                    
                    points[k][0] = X[v1[0]];
                    points[k][1] = Y[v1[1]] + ((isoVal - F[GetPointIndex(v1, dims)]) / (F[GetPointIndex(v2, dims)] - F[GetPointIndex(v1, dims)])) * (Y[v2[1]] - Y[v1[1]]);
                    
                }else if(edges[k] == 2){
                    v1[0] = ind[0];
                    v1[1] = ind[1]+1;
                    v2[0] = ind[0]+1;
                    v2[1] = ind[1]+1;
                    
                    points[k][1] = Y[v1[1]];
                    points[k][0] = X[v1[0]] + ((isoVal - F[GetPointIndex(v1, dims)]) / (F[GetPointIndex(v2, dims)] - F[GetPointIndex(v1, dims)])) * (X[v2[0]] - X[v1[0]]);
                    
                }else if(edges[k] == 3){
                    v1[0] = ind[0];
                    v1[1] = ind[1];
                    v2[0] = ind[0];
                    v2[1] = ind[1]+1;
                    
                    points[k][0] = X[v1[0]];
                    points[k][1] = Y[v1[1]] + ((isoVal - F[GetPointIndex(v1, dims)]) / (F[GetPointIndex(v2, dims)] - F[GetPointIndex(v1, dims)])) * (Y[v2[1]] - Y[v1[1]]);
                }
            }
            sl.AddSegment(points[0][0], points[0][1], points[1][0], points[1][1]); // Add segment (-10,-10) -> (+10, -10)
            
            //int edge2 = lup[cellCase][2*i+1];
            //float point2[2] = {};
            
        }
        
    }
    
    //ALGORITHM IMPLEMENTED BY JACK ENDS HERE
    
    
    
    
    
    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

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
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}

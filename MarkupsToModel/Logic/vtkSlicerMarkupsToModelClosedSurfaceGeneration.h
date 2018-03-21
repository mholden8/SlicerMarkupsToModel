#ifndef __vtkSlicerMarkupsToModelClosedSurfaceGeneration_h
#define __vtkSlicerMarkupsToModelClosedSurfaceGeneration_h

#include "vtkMRMLMarkupsToModelNode.h"

// vtk includes
#include <vtkGlyph3D.h>
#include <vtkMatrix4x4.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include "vtkSlicerMarkupsToModelModuleLogicExport.h"

class VTK_SLICER_MARKUPSTOMODEL_MODULE_LOGIC_EXPORT vtkSlicerMarkupsToModelClosedSurfaceGeneration : public vtkObject
{
  public:
    // standard vtk object methods
    vtkTypeMacro( vtkSlicerMarkupsToModelClosedSurfaceGeneration, vtkObject );
    void PrintSelf( ostream& os, vtkIndent indent ) VTK_OVERRIDE;
    static vtkSlicerMarkupsToModelClosedSurfaceGeneration *New();

    enum PointArrangement
    {
      POINT_ARRANGEMENT_SINGULAR = 0,
      POINT_ARRANGEMENT_LINEAR,
      POINT_ARRANGEMENT_PLANAR,
      POINT_ARRANGEMENT_NONPLANAR,
      POINT_ARRANGEMENT_LAST // do not set to this type, insert valid types above this line
    };

    // Generates the closed surface from the points using vtkDelaunay3D.
    static bool GenerateDelaunayClosedSurfaceModel( vtkPoints* points, vtkPolyData* outputPolyData, double delaunayAlpha, bool smoothing, bool forceConvex );

    // Generates the closed surface from the points using vtkDelaunay2D with extrusion
    static bool GenerateExtrusionClosedSurfaceModel( vtkPoints* points, vtkPolyData* outputPolyData, double extrusionDepth );

  protected:
    vtkSlicerMarkupsToModelClosedSurfaceGeneration();
    ~vtkSlicerMarkupsToModelClosedSurfaceGeneration();

  private:
    // Compute the best fit plane through the points, as well as the major and minor axes which describe variation in points.
    static void ComputeTransformMatrixFromBoundingAxes( vtkPoints* points, vtkMatrix4x4* transformFromBoundingAxes );

    // Compute the range of points along the specified axes (total lengths along which points appear)
    static void ComputeTransformedExtentRanges( vtkPoints* points, vtkMatrix4x4* transformMatrix, double outputExtentRanges[ 3 ] );

    // Compute the amount to extrude surfaces when closed surface is linear or planar.
    static double ComputeDegenerateSurfaceExtrusionAmount( const double extents[ 3 ] );

    // Compute the 2D surface poly data from a set of points
    static void GetSurfacePolyData( vtkPoints* inputPoints, vtkPolyData* surfacePolyData );

    // Compute the 2D extruded poly data from a surface and a translation vector
    static void GetExtrudedPolyData( vtkPolyData* surfacePolyData, vtkPolyData* extrudedPolyData, double surfaceNormal[ 3 ], double extrusionDepth );

    // Find the boundary points from a 2D poly data surface
    static void GetPolyDataBoundaryPoints( vtkPolyData* polyData, vtkPoints* boundaryPoints );

    // Compute the 2D surface poly data from a set of points
    static void GetSidePolyData( vtkPoints* surfaceBoundaryPoints, vtkPoints* extrudedBoundaryPoints, vtkPolyData* sidePolyData );

    // Join together and clean up 3 poly datas: the 2D top surface poly data, the 2d extruded surface poly data, and the side poly data
    static void JoinExtrusionPolyData( vtkPolyData* surfacePolyData, vtkPolyData* extrudedPolyData, vtkPolyData* sidePolyData, vtkPolyData* joinedPolyData );

    // Compute the 2D surface poly data from a set of points
    static bool IsPolyDataClosed( vtkPolyData* polyData );

    // Find out what kind of arrangment the points are in (see PointArrangementEnum above).
    // If the arrangement is planar, stores the normal of the best fit plane in planeNormal.
    // If the arrangement is linear, stores the axis of the best fit line in lineAxis.
    static PointArrangement ComputePointArrangement( const double smallestBoundingExtentRanges[ 3 ] );

    // Convert the point set into a poly data
    static void GetPolyDataFromPoints( vtkPoints* inputPoints, vtkPolyData* outputPolyData );

    // Generate basic glyphs for degenerate point sets
    static void GetPointArrangementSingularGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, int numberOfPoints );
    static void GetPointArrangementLinearGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, double smallestBoundingExtentRanges[3], vtkMatrix4x4* boundingAxesToRasTransformMatrix );
    static void GetPointArrangementPlanarGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, double smallestBoundingExtentRanges[3], vtkMatrix4x4* boundingAxesToRasTransformMatrix );
    

    // helper utility functions
    static void SetNthColumnInMatrix( vtkMatrix4x4* matrix, int n, const double axis[ 3 ] );
    static void GetNthColumnInMatrix( vtkMatrix4x4* matrix, int n, double outputAxis[ 3 ] );

    // not used
    vtkSlicerMarkupsToModelClosedSurfaceGeneration ( const vtkSlicerMarkupsToModelClosedSurfaceGeneration& ) VTK_DELETE_FUNCTION;
    void operator= ( const vtkSlicerMarkupsToModelClosedSurfaceGeneration& ) VTK_DELETE_FUNCTION;
};

#endif

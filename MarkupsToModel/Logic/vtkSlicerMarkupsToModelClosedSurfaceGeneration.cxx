#include "vtkSlicerMarkupsToModelClosedSurfaceGeneration.h"

#include "vtkMRMLModelNode.h"
#include "vtkMRMLMarkupsFiducialNode.h"

#include <vtkAppendPolyData.h>
#include <vtkButterflySubdivisionFilter.h>
#include <vtkCellArray.h>
#include <vtkCenterOfMass.h>
#include <vtkCleanPolyData.h>
#include <vtkCubeSource.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDelaunay2D.h>
#include <vtkGlyph3D.h>
#include <vtkFeatureEdges.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkLineSource.h>
#include <vtkNew.h>
#include <vtkOBBTree.h>
#include <vtkPolyDataNormals.h>
#include <vtkRegularPolygonSource.h>
#include <vtkStripper.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>

//------------------------------------------------------------------------------
// constants within this file
static const double COMPARE_TO_ZERO_TOLERANCE = 0.0001;
static const double MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT = 0.01; // if a surface is flat/linear, give it at least this much depth

//------------------------------------------------------------------------------
vtkStandardNewMacro( vtkSlicerMarkupsToModelClosedSurfaceGeneration );

//------------------------------------------------------------------------------
vtkSlicerMarkupsToModelClosedSurfaceGeneration::vtkSlicerMarkupsToModelClosedSurfaceGeneration()
{
}

//------------------------------------------------------------------------------
vtkSlicerMarkupsToModelClosedSurfaceGeneration::~vtkSlicerMarkupsToModelClosedSurfaceGeneration()
{
}

//------------------------------------------------------------------------------
bool vtkSlicerMarkupsToModelClosedSurfaceGeneration::GenerateDelaunayClosedSurfaceModel(vtkPoints* inputPoints, vtkPolyData* outputPolyData,
  double delaunayAlpha, bool smoothing, bool forceConvex)
{  
  if (inputPoints == NULL)
  {
    vtkGenericWarningMacro("Input points are null. No model generated.");
    return false;
  }
  
  if (outputPolyData == NULL)
  {
    vtkGenericWarningMacro("Output poly data is null. No model generated.");
    return false;
  }

  int numberOfPoints = inputPoints->GetNumberOfPoints();
  if (numberOfPoints == 0)
  {
    // No markup points, the output should be empty
    return true;
  }

  vtkSmartPointer< vtkCellArray > inputCellArray = vtkSmartPointer< vtkCellArray >::New();
  inputCellArray->InsertNextCell(numberOfPoints);
  for (int i = 0; i < numberOfPoints; i++)
  {
    inputCellArray->InsertCellPoint(i);
  }

  vtkSmartPointer< vtkPolyData > inputPolyData = vtkSmartPointer< vtkPolyData >::New();
  inputPolyData->SetLines(inputCellArray);
  inputPolyData->SetPoints(inputPoints);

  vtkSmartPointer< vtkDelaunay3D > delaunay = vtkSmartPointer< vtkDelaunay3D >::New();
  delaunay->SetAlpha(delaunayAlpha);
  delaunay->AlphaTrisOff();
  delaunay->AlphaLinesOff();
  delaunay->AlphaVertsOff();

  vtkSmartPointer< vtkMatrix4x4 > boundingAxesToRasTransformMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  ComputeTransformMatrixFromBoundingAxes(inputPoints, boundingAxesToRasTransformMatrix);

  vtkSmartPointer< vtkMatrix4x4 > rasToBoundingAxesTransformMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  vtkMatrix4x4::Invert(boundingAxesToRasTransformMatrix, rasToBoundingAxesTransformMatrix);

  double smallestBoundingExtentRanges[3] = { 0.0, 0.0, 0.0 }; // temporary values
  ComputeTransformedExtentRanges(inputPoints, rasToBoundingAxesTransformMatrix, smallestBoundingExtentRanges);

  PointArrangement pointArrangement = ComputePointArrangement(smallestBoundingExtentRanges);

  switch (pointArrangement)
  {
    case POINT_ARRANGEMENT_SINGULAR:
    {
      vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
      GetPointArrangementSingularGlyph(inputPolyData, glyph, numberOfPoints);
      delaunay->SetInputConnection(glyph->GetOutputPort());
      break;
    }
    case POINT_ARRANGEMENT_LINEAR:
    {
      vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
      GetPointArrangementLinearGlyph(inputPolyData, glyph, smallestBoundingExtentRanges, boundingAxesToRasTransformMatrix);
      delaunay->SetInputConnection(glyph->GetOutputPort());
      break;
    }
    case POINT_ARRANGEMENT_PLANAR:
    {
      vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
      GetPointArrangementPlanarGlyph(inputPolyData, glyph, smallestBoundingExtentRanges, boundingAxesToRasTransformMatrix);
      delaunay->SetInputConnection(glyph->GetOutputPort());
      break;
    }
    case POINT_ARRANGEMENT_NONPLANAR:
    {
      delaunay->SetInputData(inputPolyData);
      break;
    }
    default: // unsupported or invalid
    {
      vtkGenericWarningMacro("Unsupported pointArrangementType detected: " << pointArrangement << ". Aborting closed surface generation.");
      return false;
    }
  }

  vtkSmartPointer< vtkDataSetSurfaceFilter > surfaceFilter = vtkSmartPointer< vtkDataSetSurfaceFilter >::New();
  surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
  surfaceFilter->Update();

  vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetFeatureAngle(100); // TODO: This needs some justification, or set as an input parameter

  if (smoothing && pointArrangement == POINT_ARRANGEMENT_NONPLANAR)
  {
    vtkSmartPointer< vtkButterflySubdivisionFilter > subdivisionFilter = vtkSmartPointer< vtkButterflySubdivisionFilter >::New();
    subdivisionFilter->SetInputConnection(surfaceFilter->GetOutputPort());
    subdivisionFilter->SetNumberOfSubdivisions(3);
    subdivisionFilter->Update();
    if (forceConvex)
    {
      vtkSmartPointer< vtkDelaunay3D > convexHull = vtkSmartPointer< vtkDelaunay3D >::New();
      convexHull->SetInputConnection(subdivisionFilter->GetOutputPort());
      convexHull->Update();
      vtkSmartPointer< vtkDataSetSurfaceFilter > surfaceFilter = vtkSmartPointer< vtkDataSetSurfaceFilter >::New();
      surfaceFilter->SetInputData(convexHull->GetOutput());
      surfaceFilter->Update();
      normals->SetInputConnection(surfaceFilter->GetOutputPort());
    }
    else
    {
      normals->SetInputConnection(subdivisionFilter->GetOutputPort());
    }
  }
  else
  {
    vtkNew<vtkLinearSubdivisionFilter> linearSubdivision;
    linearSubdivision->SetInputConnection(surfaceFilter->GetOutputPort());
    normals->SetInputConnection(linearSubdivision->GetOutputPort());
  }
  normals->Update();

  outputPolyData->DeepCopy(normals->GetOutput());
  return true;
}

//------------------------------------------------------------------------------
bool vtkSlicerMarkupsToModelClosedSurfaceGeneration::GenerateExtrusionClosedSurfaceModel(vtkPoints* inputPoints, vtkPolyData* outputPolyData,
  double extrusionDepth, bool forcePlanar)
{  
  if (inputPoints == NULL)
  {
    vtkGenericWarningMacro("Input points are null. No model generated.");
    return false;
  }
  
  if (outputPolyData == NULL)
  {
    vtkGenericWarningMacro("Output poly data is null. No model generated.");
    return false;
  }

  int numberOfPoints = inputPoints->GetNumberOfPoints();
  if (numberOfPoints == 0)
  {
    // No markup points, the output should be empty
    return true;
  }

  vtkSmartPointer< vtkPolyData > inputPolyData = vtkSmartPointer< vtkPolyData >::New();
  GetPolyDataFromPoints( inputPoints, inputPolyData );

  vtkSmartPointer< vtkMatrix4x4 > boundingAxesToRasTransformMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  ComputeTransformMatrixFromBoundingAxes(inputPoints, boundingAxesToRasTransformMatrix);

  vtkSmartPointer< vtkMatrix4x4 > rasToBoundingAxesTransformMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  vtkMatrix4x4::Invert(boundingAxesToRasTransformMatrix, rasToBoundingAxesTransformMatrix);

  double smallestBoundingExtentRanges[3] = { 0.0, 0.0, 0.0 }; // temporary values
  ComputeTransformedExtentRanges(inputPoints, rasToBoundingAxesTransformMatrix, smallestBoundingExtentRanges);

  PointArrangement pointArrangement = ComputePointArrangement(smallestBoundingExtentRanges);

  // Deal with cases of degenerate point arrangements
  if (pointArrangement == POINT_ARRANGEMENT_SINGULAR || pointArrangement == POINT_ARRANGEMENT_LINEAR)
  {
    return GenerateDelaunayClosedSurfaceModel(inputPoints, outputPolyData, 0.0, false, true); // TODO: Should this be handled otherwise?
  }

  // Find the normal vector to the surface (this is the direction in which we will extrude)
  double surfaceNormal[3] = { 0.0, 0.0, 0.0 }; // temporary values
  const int PLANE_NORMAL_INDEX = 2; // The plane normal has the smallest variation, and is stored in the last column
  GetNthColumnInMatrix(boundingAxesToRasTransformMatrix, PLANE_NORMAL_INDEX, surfaceNormal);

  // Compute the poly data on the collected surface of the model
  vtkSmartPointer<vtkPolyData> surfacePolyData = vtkSmartPointer<vtkPolyData>::New();
  if (forcePlanar)
  {
    GetPlanarSurfacePolyData(inputPoints, surfacePolyData, surfaceNormal);
  }
  else
  {
    GetSurfacePolyData(inputPoints, surfacePolyData);
  }

  // Compute the poly data that is on the extruded surface
  vtkSmartPointer<vtkPolyData> extrudedPolyData = vtkSmartPointer<vtkPolyData>::New();
  GetExtrudedPolyData(surfacePolyData, extrudedPolyData, surfaceNormal, extrusionDepth);

  // Compute the boundary points on each surface - these will be used to join the surfaces
  vtkSmartPointer<vtkPoints> surfaceBoundaryPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> extrudedBoundaryPoints = vtkSmartPointer<vtkPoints>::New();
  GetPolyDataBoundaryPoints(surfacePolyData, surfaceBoundaryPoints);
  GetPolyDataBoundaryPoints(extrudedPolyData, extrudedBoundaryPoints);

  // If the boundary points could not be created, this is because the pointset is degenerate
  if( surfaceBoundaryPoints->GetNumberOfPoints() == 0 || extrudedBoundaryPoints->GetNumberOfPoints() == 0 )
  {
    vtkGenericWarningMacro("Could not create surface boundary. Possibly due to point arrangement on coordinate plane. Aborting closed surface generation.");
    return false;
  }

  // Create the poly data for the "sides" of the surface (joining the top surface and extruded surface)
  vtkSmartPointer<vtkPolyData> sidePolyData = vtkSmartPointer<vtkPolyData>::New();
  GetSidePolyData(surfaceBoundaryPoints, extrudedBoundaryPoints, sidePolyData);

  // Finally, put it all together into a closed surface model
  vtkSmartPointer<vtkPolyData> joinedPolyData = vtkSmartPointer<vtkPolyData>::New();
  JoinExtrusionPolyData(surfacePolyData, extrudedPolyData, sidePolyData, joinedPolyData);

  if(!IsPolyDataClosed(joinedPolyData))
  {
    vtkGenericWarningMacro("Could not create closed surface. Aborting closed surface generation.");
    return false;
  }

  outputPolyData->DeepCopy(joinedPolyData);
  return true;
}

//------------------------------------------------------------------------------
// Compute the principal axes of the point cloud. The x axis represents the axis
// with maximum variation, and the z axis has minimum variation.
// This function is currently implemented using the vtkOBBTree object.
// There are two limitations with this approach:
// 1. vtkOBBTree may have a performance impact
// 2. The axes returned are based on variation of coordinates, not the range
//    (so the return result is not necessarily intuitive, variation != length).
// Neither of these limitations will prevent the overall logic from functioning
// correctly, but it is worth keeping in mind, and worth changing should a need 
// arise
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::ComputeTransformMatrixFromBoundingAxes(vtkPoints* points, vtkMatrix4x4* boundingAxesToRasTransformMatrix)
{
  if (points == NULL)
  {
    vtkGenericWarningMacro("Points object is null. Cannot compute best fit planes.");
    return;
  }

  if (boundingAxesToRasTransformMatrix == NULL)
  {
    vtkGenericWarningMacro("Output matrix object is null. Cannot compute best fit planes.");
    return;
  }

  // the output matrix should start as identity, so no translation etc.
  boundingAxesToRasTransformMatrix->Identity();

  // Compute the plane using the smallest bounding box that can have arbitrary axes
  vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();
  double cornerOBBOrigin[3] = { 0.0, 0.0, 0.0 }; // unused
  double variationMaximumOBBAxis[3] = { 0.0, 0.0, 0.0 };
  double variationMediumOBBAxis[3] = { 0.0, 0.0, 0.0 };
  double variationMinimumOBBAxis[3] = { 0.0, 0.0, 0.0 };
  double relativeAxisSizes[3] = { 0.0, 0.0, 0.0 }; // unused, the values represented herein are unclear
  obbTree->ComputeOBB(points, cornerOBBOrigin, variationMaximumOBBAxis, variationMediumOBBAxis, variationMinimumOBBAxis, relativeAxisSizes);

  // now to store the desired results in the appropriate axis of the output matrix.
  // must check each axis to make sure it was actually computed (non-zero)
  // do the maxmimum variation axis
  if (vtkMath::Norm(variationMaximumOBBAxis) < COMPARE_TO_ZERO_TOLERANCE)
  {
    // there is no variation in the points whatsoever.
    // i.e. all points are in a single position.
    // return arbitrary orthonormal axes (the standard axes will do).
    boundingAxesToRasTransformMatrix->Identity();
    return;
  }
  vtkMath::Normalize(variationMaximumOBBAxis);
  SetNthColumnInMatrix(boundingAxesToRasTransformMatrix, 0, variationMaximumOBBAxis);

  // do the medium variation axis
  if (vtkMath::Norm(variationMediumOBBAxis) < COMPARE_TO_ZERO_TOLERANCE)
  {
    // the points are colinear along only the maximum axis
    // any two perpendicular orthonormal vectors will do for the remaining axes.
    double thetaAngle = 0.0; // this can be arbitrary
    vtkMath::Perpendiculars(variationMaximumOBBAxis, variationMediumOBBAxis, variationMinimumOBBAxis, thetaAngle);
  }
  vtkMath::Normalize(variationMediumOBBAxis);
  SetNthColumnInMatrix(boundingAxesToRasTransformMatrix, 1, variationMediumOBBAxis);

  // do the minimum variation axis
  if (vtkMath::Norm(variationMinimumOBBAxis) < COMPARE_TO_ZERO_TOLERANCE)
  {
    // all points lie exactly on a plane.
    // the remaining perpendicular vector found using cross product.
    vtkMath::Cross(variationMaximumOBBAxis, variationMediumOBBAxis, variationMinimumOBBAxis);
  }
  vtkMath::Normalize(variationMinimumOBBAxis);
  SetNthColumnInMatrix(boundingAxesToRasTransformMatrix, 2, variationMinimumOBBAxis);
}

//------------------------------------------------------------------------------
// It is assumed that sortedExtentRanges is pre-sorted in descending order (largest to smallest)
vtkSlicerMarkupsToModelClosedSurfaceGeneration::PointArrangement vtkSlicerMarkupsToModelClosedSurfaceGeneration::ComputePointArrangement(const double sortedExtentRanges[3])
{
  if (sortedExtentRanges == NULL)
  {
    vtkGenericWarningMacro("Input sortedExtentRanges is null. Returning singularity result.");
    return POINT_ARRANGEMENT_SINGULAR;
  }

  double longestExtentRange = sortedExtentRanges[0];
  double mediumExtentRange = sortedExtentRanges[1];
  double shortestExtentRange = sortedExtentRanges[2];

  // sanity checking
  bool longestExtentSmallerThanMedium = longestExtentRange >= COMPARE_TO_ZERO_TOLERANCE && longestExtentRange < mediumExtentRange;
  bool longestExtentSmallerThanShortest = longestExtentRange >= COMPARE_TO_ZERO_TOLERANCE && longestExtentRange < shortestExtentRange;
  bool mediumExtentSmallerThanShortest = mediumExtentRange >= COMPARE_TO_ZERO_TOLERANCE && mediumExtentRange < shortestExtentRange;
  if (longestExtentSmallerThanMedium || longestExtentSmallerThanShortest || mediumExtentSmallerThanShortest)
  {
    // Don't correct the problem here. Code external to this function should pass
    // extent ranges already sorted, so it indicates a problem elsewhere.
    vtkGenericWarningMacro("Extent ranges not provided in order largest to smallest. Unexpected results may occur.");
  }

  if (longestExtentRange < COMPARE_TO_ZERO_TOLERANCE)
  {
    return POINT_ARRANGEMENT_SINGULAR;
  }

  // We need to compare relative lengths of the short and medium axes against
  // the longest axis.
  double mediumToLongestRatio = mediumExtentRange / longestExtentRange;

  // The Delaunay3D class tends to fail with thin planes/lines, so it is important
  // to capture these cases, even liberally. It was experimentally determined that
  // extents less than 1/10th of the maximum extent tend to produce errors.
  const double RATIO_THRESHOLD = 0.1;

  if (mediumToLongestRatio < RATIO_THRESHOLD)
  {
    return POINT_ARRANGEMENT_LINEAR;
  }

  double shortestToLongestRatio = shortestExtentRange / longestExtentRange;
  if (shortestToLongestRatio < RATIO_THRESHOLD)
  {
    return POINT_ARRANGEMENT_PLANAR;
  }

  return POINT_ARRANGEMENT_NONPLANAR;
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::ComputeTransformedExtentRanges(vtkPoints* points, vtkMatrix4x4* transformMatrix, double outputExtentRanges[3])
{
  if (points == NULL)
  {
    vtkGenericWarningMacro("points is null. Aborting output extent computation.");
    return;
  }

  if (transformMatrix == NULL)
  {
    vtkGenericWarningMacro("transformMatrix is null. Aborting output extent computation.");
    return;
  }

  if (outputExtentRanges == NULL)
  {
    vtkGenericWarningMacro("outputExtentRanges is null. Aborting output extent computation.");
    return;
  }

  vtkSmartPointer< vtkTransform > transform = vtkSmartPointer< vtkTransform >::New();
  transform->SetMatrix(transformMatrix);
  transform->Update();

  // can't transform points directly, so need to store in a container
  vtkSmartPointer< vtkPolyData > polyDataWithPoints = vtkSmartPointer< vtkPolyData >::New();
  polyDataWithPoints->SetPoints(points);

  vtkSmartPointer< vtkTransformFilter > transformFilter = vtkSmartPointer< vtkTransformFilter >::New();
  transformFilter->SetTransform(transform);
  transformFilter->SetInputData(polyDataWithPoints);
  transformFilter->Update();

  // the extent can be extracted from the output points object (poly data bounds does not work)
  vtkPoints* transformedPoints = transformFilter->GetPolyDataOutput()->GetPoints();
  transformedPoints->ComputeBounds();
  double* extents = transformedPoints->GetBounds(); // { xmin, xmax, ymin, ymax, zmin, zmax }

  for (int i = 0; i < 3; i++)
  {
    double axisIMin = extents[2 * i];
    double axisIMax = extents[2 * i + 1];
    double axisIRange = axisIMax - axisIMin;
    outputExtentRanges[i] = axisIRange;
  }
}

//------------------------------------------------------------------------------
double vtkSlicerMarkupsToModelClosedSurfaceGeneration::ComputeDegenerateSurfaceExtrusionAmount(const double extents[3])
{
  // MINIMUM_SURFACE_EXTRUSION_AMOUNT is the value returned by default, and the final result cannot be less than this.
  if (extents == NULL)
  {
    vtkGenericWarningMacro("extents is null. Returning MINIMUM_SURFACE_EXTRUSION_AMOUNT: " << MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT << ".");
    return MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT;
  }

  double normOfExtents = vtkMath::Norm(extents);
  const double SURFACE_EXTRUSION_NORM_MULTIPLIER = 0.01; // this value is observed to produce generally acceptable results
  double surfaceExtrusionAmount = normOfExtents * SURFACE_EXTRUSION_NORM_MULTIPLIER;

  if (surfaceExtrusionAmount < MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT)
  {
    vtkGenericWarningMacro("Surface extrusion amount smaller than " << MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT << " : " << surfaceExtrusionAmount << ". "
      << "Consider checking the points for singularity. Setting surface extrusion amount to default "
      << MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT << ".");
    surfaceExtrusionAmount = MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT;
  }
  return surfaceExtrusionAmount;
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetSurfacePolyData(vtkPoints* inputPoints, vtkPolyData* surfacePolyData)
{
  vtkNew<vtkPolyData> surfacePointsPolyData;
  surfacePointsPolyData->SetPoints(inputPoints);

  // Use 2D delaunay triangulation to find the top surface
  vtkNew<vtkDelaunay2D> delaunayFilter;
  delaunayFilter->SetProjectionPlaneMode(VTK_BEST_FITTING_PLANE);
  delaunayFilter->SetInputData(surfacePointsPolyData.GetPointer());

  // Clean the surface
  vtkNew<vtkCleanPolyData> cleanSurfaceFilter;
  cleanSurfaceFilter->SetTolerance( COMPARE_TO_ZERO_TOLERANCE );
  cleanSurfaceFilter->SetInputConnection(delaunayFilter->GetOutputPort());

  cleanSurfaceFilter->Update();  
  surfacePolyData->DeepCopy(cleanSurfaceFilter->GetOutput());
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPlanarSurfacePolyData(vtkPoints* inputPoints, vtkPolyData* surfacePolyData, double surfaceNormal[ 3 ])
{
  vtkNew<vtkPolyData> surfacePointsPolyData;
  surfacePointsPolyData->SetPoints(inputPoints);

  // Compute the centre of mass
  vtkNew<vtkCenterOfMass> centerOfMassFilter;
  centerOfMassFilter->SetInputData(surfacePointsPolyData.GetPointer());
  centerOfMassFilter->SetUseScalarsAsWeights(false);
  centerOfMassFilter->Update();
  double centerOfMass[3] = {0.0, 0.0, 0.0};
  centerOfMassFilter->GetCenter(centerOfMass);

  // Find the projection of each point onto the plane of best fit (defined by center of mass and normal vector)
  vtkNew<vtkPoints> planePoints;
  for (int i = 0; i < inputPoints->GetNumberOfPoints(); i++)
  {
    double currentPoint[3] = {0.0, 0.0, 0.0};
    inputPoints->GetPoint(i, currentPoint);
    double currentRelativePoint[3] = {0.0, 0.0, 0.0};
    vtkMath::Subtract(currentPoint, centerOfMass, currentRelativePoint);
    double normalLength = vtkMath::Dot(currentRelativePoint, surfaceNormal);
    double normalComponent[3] = {surfaceNormal[0], surfaceNormal[1], surfaceNormal[2]};
    vtkMath::MultiplyScalar(normalComponent, normalLength);
    double currentRelativePlanePoint[3] = {0.0, 0.0, 0.0};
    vtkMath::Subtract(currentRelativePoint, normalComponent, currentRelativePlanePoint);
    double currentPlanePoint[3] = {0.0, 0.0, 0.0};
    vtkMath::Add(centerOfMass, currentRelativePlanePoint, currentPlanePoint);

    planePoints->InsertNextPoint(currentPlanePoint);
  }

  // Use the regular surface polydata method on the projected points
  GetSurfacePolyData(planePoints.GetPointer(), surfacePolyData);
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetExtrudedPolyData(vtkPolyData* surfacePolyData, vtkPolyData* extrudedPolyData, double surfaceNormal[ 3 ], double extrusionDepth)
{
  // Get the translation vector for the extruded surface poly data
  vtkMath::MultiplyScalar( surfaceNormal, extrusionDepth );

  vtkNew<vtkTransform> extrusionTransform;
  extrusionTransform->Translate(surfaceNormal);

  vtkNew<vtkTransformPolyDataFilter> extrusionTransformFilter;
  extrusionTransformFilter->SetInputData(surfacePolyData);
  extrusionTransformFilter->SetTransform(extrusionTransform.GetPointer());

  extrusionTransformFilter->Update();
  extrudedPolyData->DeepCopy(extrusionTransformFilter->GetOutput());
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPolyDataBoundaryPoints(vtkPolyData* polyData, vtkPoints* boundaryPoints)
{
  // Find the boundary edges
  vtkNew<vtkFeatureEdges> featureEdgesFilter;
  featureEdgesFilter->FeatureEdgesOff();
  featureEdgesFilter->NonManifoldEdgesOff();
  featureEdgesFilter->ManifoldEdgesOff();
  featureEdgesFilter->BoundaryEdgesOn();
  featureEdgesFilter->SetInputData(polyData);

  // Cleanup the boundary edges so that we just have the relevant points
  vtkNew<vtkStripper> stripperFilter;
  stripperFilter->SetInputConnection(featureEdgesFilter->GetOutputPort());
  
  vtkNew<vtkCleanPolyData> cleanPointsFilter;
  cleanPointsFilter->SetInputConnection(stripperFilter->GetOutputPort());

  cleanPointsFilter->Update();
  boundaryPoints->DeepCopy(cleanPointsFilter->GetOutput()->GetPoints());
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetSidePolyData(vtkPoints* surfaceBoundaryPoints, vtkPoints* extrudedBoundaryPoints, vtkPolyData* sidePolyData)
{
  // If there are different numbers of boundary points for the two poly datas, then we cannot join them. Abort.
  if (surfaceBoundaryPoints->GetNumberOfPoints() != extrudedBoundaryPoints->GetNumberOfPoints())
  {
    vtkGenericWarningMacro("Top and deep surfaces have different number of boundary points.");
    return;
  }
  int numberOfBoundaryPoints = surfaceBoundaryPoints->GetNumberOfPoints();

  vtkNew<vtkAppendPolyData> sidePolyDataAppend;

  for ( int i = 0; i < numberOfBoundaryPoints; i++ )
  {
    vtkNew<vtkPoints> currentStripPoints;

    double point1[ 3 ] = { 0, 0, 0 }; // From top surface
    surfaceBoundaryPoints->GetPoint( i % numberOfBoundaryPoints, point1 );
    currentStripPoints->InsertNextPoint( point1 );

    double point2[ 3 ] = { 0, 0, 0 }; // From extruded surface
    extrudedBoundaryPoints->GetPoint( i % numberOfBoundaryPoints, point1 );
    currentStripPoints->InsertNextPoint( point1 );
    
    double point3[ 3 ] = { 0, 0, 0 }; // From top surface
    surfaceBoundaryPoints->GetPoint( ( i + 1 ) % numberOfBoundaryPoints, point1 );
    currentStripPoints->InsertNextPoint( point1 );

    double point4[ 3 ] = { 0, 0, 0 }; // From extruded surface
    extrudedBoundaryPoints->GetPoint( ( i + 1 ) % numberOfBoundaryPoints, point1 );
    currentStripPoints->InsertNextPoint( point1 );
      
    // We know what the triangles should look like, so create them manually
    vtkNew<vtkTriangle> triangle1;
    triangle1->GetPointIds()->SetId( 0, 0 );
    triangle1->GetPointIds()->SetId( 1, 1 );
    triangle1->GetPointIds()->SetId( 2, 2 );
      
    vtkNew<vtkTriangle> triangle2;
    triangle2->GetPointIds()->SetId( 0, 3 );
    triangle2->GetPointIds()->SetId( 1, 2 );
    triangle2->GetPointIds()->SetId( 2, 1 );
      
    vtkNew<vtkCellArray> triangles;
    triangles->InsertNextCell(triangle1.GetPointer());
    triangles->InsertNextCell(triangle2.GetPointer());
      
    vtkNew<vtkPolyData> currentStripPolyData;
    currentStripPolyData->SetPoints(currentStripPoints.GetPointer());
    currentStripPolyData->SetPolys(triangles.GetPointer());
         
    sidePolyDataAppend->AddInputData(currentStripPolyData.GetPointer());
  }
    
  sidePolyDataAppend->Update();
  sidePolyData->DeepCopy(sidePolyDataAppend->GetOutput());
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::JoinExtrusionPolyData(vtkPolyData* surfacePolyData, vtkPolyData* extrudedPolyData, vtkPolyData* sidePolyData, vtkPolyData* joinedPolyData)
{
  // Append the poly data into one single poly data Put the poly data together
  vtkNew<vtkAppendPolyData> appendDirtyPolyData;
  appendDirtyPolyData->AddInputData(surfacePolyData);
  appendDirtyPolyData->AddInputData(extrudedPolyData);
  appendDirtyPolyData->AddInputData(sidePolyData);
    
  // Clean up so the new poly data so that it is closed and the normal face inward
  vtkNew<vtkCleanPolyData> preNormalCleanFilter;
  preNormalCleanFilter->SetInputConnection(appendDirtyPolyData->GetOutputPort());

  vtkNew<vtkPolyDataNormals> normalsFilter;
  normalsFilter->AutoOrientNormalsOn();
  normalsFilter->SetInputConnection(preNormalCleanFilter->GetOutputPort());

  vtkNew<vtkCleanPolyData> postNormalCleanFilter;
  postNormalCleanFilter->SetInputConnection(normalsFilter->GetOutputPort());

  postNormalCleanFilter->Update();
  joinedPolyData->DeepCopy(postNormalCleanFilter->GetOutput());
}

//------------------------------------------------------------------------------
bool vtkSlicerMarkupsToModelClosedSurfaceGeneration::IsPolyDataClosed(vtkPolyData* polyData)
{
 vtkNew<vtkFeatureEdges> edgesFilter;
  edgesFilter->FeatureEdgesOff();
  edgesFilter->BoundaryEdgesOn();
  edgesFilter->NonManifoldEdgesOn();
  edgesFilter->SetInputData(polyData);
  edgesFilter->Update();
    
  return (edgesFilter->GetOutput()->GetNumberOfCells() == 0);
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::SetNthColumnInMatrix(vtkMatrix4x4* matrix, int n, const double axis[3])
{
  if (matrix == NULL)
  {
    vtkGenericWarningMacro("No matrix provided as input. No operation performed.");
    return;
  }

  if (n < 0 || n >= 3)
  {
    vtkGenericWarningMacro("Axis n " << n << " is out of bounds. Valid values are 0, 1, and 2. No operation performed.");
    return;
  }

  if (axis == NULL)
  {
    vtkGenericWarningMacro("Axis is null. No operation performed.");
    return;
  }

  matrix->SetElement(0, n, axis[0]);
  matrix->SetElement(1, n, axis[1]);
  matrix->SetElement(2, n, axis[2]);
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetNthColumnInMatrix(vtkMatrix4x4* matrix, int n, double outputAxis[3])
{
  if (matrix == NULL)
  {
    vtkGenericWarningMacro("No matrix provided as input. No operation performed.");
    return;
  }

  if (n < 0 || n >= 3)
  {
    vtkGenericWarningMacro("Axis n " << n << " is out of bounds. Valid values are 0, 1, and 2. No operation performed.");
    return;
  }

  if (outputAxis == NULL)
  {
    vtkGenericWarningMacro("Axis is null. No operation performed.");
    return;
  }

  outputAxis[0] = matrix->GetElement(0, n);
  outputAxis[1] = matrix->GetElement(1, n);
  outputAxis[2] = matrix->GetElement(2, n);
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPolyDataFromPoints( vtkPoints* inputPoints, vtkPolyData* outputPolyData )
{
  int numberOfPoints = inputPoints->GetNumberOfPoints();

  vtkSmartPointer< vtkCellArray > inputCellArray = vtkSmartPointer< vtkCellArray >::New();
  inputCellArray->InsertNextCell(numberOfPoints);
  for (int i = 0; i < numberOfPoints; i++)
  {
    inputCellArray->InsertCellPoint(i);
  }

  outputPolyData->SetLines(inputCellArray);
  outputPolyData->SetPoints(inputPoints);
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPointArrangementSingularGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, int numberOfPoints )
{
  vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
  // there is only one point, we cannot compute extent or extrusion from this.
  double extrusionMagnitude = MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT;
  if ( numberOfPoints > 1 )
  {
    vtkGenericWarningMacro( "There is more than one input point, but they form a singularity. " <<
                            "Giving depth of " << MINIMUM_DEGENERATE_SURFACE_EXTRUSION_AMOUNT << "." );
  }
  cubeSource->SetBounds(-extrusionMagnitude, extrusionMagnitude,
    -extrusionMagnitude, extrusionMagnitude,
    -extrusionMagnitude, extrusionMagnitude);

  outputGlyph->SetSourceConnection(cubeSource->GetOutputPort());
  outputGlyph->SetInputData(inputPolyData);
  outputGlyph->Update();
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPointArrangementLinearGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, double smallestBoundingExtentRanges[3], vtkMatrix4x4* boundingAxesToRasTransformMatrix )
{
  // draw a "square" around the line (make it a rectangular prism)
  vtkSmartPointer<vtkRegularPolygonSource> squareSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
  squareSource->SetCenter(0.0, 0.0, 0.0);
  double extrusionMagnitude = ComputeDegenerateSurfaceExtrusionAmount(smallestBoundingExtentRanges); // need to give some depth
  squareSource->SetRadius(extrusionMagnitude);
  squareSource->SetNumberOfSides(4);
  double lineAxis[3] = { 0.0, 0.0, 0.0 }; // temporary values
  const int LINE_AXIS_INDEX = 0; // The largest (and only meaningful) axis is in the 0th column
  // the bounding axes are stored in the columns of transformFromBoundingAxes
  GetNthColumnInMatrix(boundingAxesToRasTransformMatrix, LINE_AXIS_INDEX, lineAxis);
  squareSource->SetNormal(lineAxis);

  outputGlyph->SetSourceConnection(squareSource->GetOutputPort());
  outputGlyph->SetInputData(inputPolyData);
  outputGlyph->Update();
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::GetPointArrangementPlanarGlyph( vtkPolyData* inputPolyData, vtkGlyph3D* outputGlyph, double smallestBoundingExtentRanges[3], vtkMatrix4x4* boundingAxesToRasTransformMatrix )
{
  // extrude additional points on either side of the plane
  vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
  double planeNormal[3] = { 0.0, 0.0, 0.0 }; // temporary values
  const int PLANE_NORMAL_INDEX = 2; // The plane normal has the smallest variation, and is stored in the last column
  // the bounding axes are stored in the columns of transformFromBoundingAxes
  GetNthColumnInMatrix(boundingAxesToRasTransformMatrix, PLANE_NORMAL_INDEX, planeNormal);
  double extrusionMagnitude = ComputeDegenerateSurfaceExtrusionAmount(smallestBoundingExtentRanges); // need to give some depth
  double point1[3] = { planeNormal[0], planeNormal[1], planeNormal[2] };
  vtkMath::MultiplyScalar(point1, extrusionMagnitude);
  lineSource->SetPoint1(point1);
  double point2[3] = { planeNormal[0], planeNormal[1], planeNormal[2] };
  vtkMath::MultiplyScalar(point2, -extrusionMagnitude);
  lineSource->SetPoint2(point2);

  outputGlyph->SetSourceConnection(lineSource->GetOutputPort());
  outputGlyph->SetInputData(inputPolyData);
  outputGlyph->Update();
}

//------------------------------------------------------------------------------
void vtkSlicerMarkupsToModelClosedSurfaceGeneration::PrintSelf( ostream &os, vtkIndent indent )
{
  Superclass::PrintSelf( os, indent );
}

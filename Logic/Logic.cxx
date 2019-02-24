#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "LogicCLP.h"

// bender includes
#include "benderIOUtils.h"

// ITK includes
#include <itkDirectory.h>
#include <itkImageTransformer.h>
#include <itkResampleImageFilter.h>

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkVolume.h>
#include <vtkNrrdReader.h>
#include <vtkImageWriter.h>

typedef itk::Matrix<double, 2, 4> Mat24;

typedef unsigned char CharType;
typedef unsigned short LabelType;

typedef itk::Image<unsigned short, 3>  LabelImage;
typedef itk::Image<bool, 3>  BoolImage;
typedef itk::Image<float, 3>  WeightImage;

typedef itk::Index<3> Voxel;
typedef itk::Offset<3> VoxelOffset;
typedef itk::ImageRegion<3> Region;

typedef itk::Matrix<double, 3, 3> Mat33;
typedef itk::Matrix<double, 4, 4> Mat44;

typedef itk::Vector<double, 3> Vec3;
typedef itk::Vector<double, 4> Vec4;

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
	template<class T>
	int DoIt( int argc, char * argv[])
	{
		  PARSE_ARGS;

		  std::cout << "Read input image";

		  vtkSmartPointer<vtkNrrdReader> reader = vtkSmartPointer<vtkNrrdReader>::New();
		  reader->SetFileName(inputVolume.c_str());
		  reader->Update();

		  std::cout << "Input Image read: fname: " << inputVolume << " result: " << reader->GetOutput() << std::endl;

		  //
		  // READ ARMATURE
		  //

		  vtkSmartPointer<vtkPolyData> armature;
		  armature.TakeReference(bender::IOUtils::ReadPolyData(ArmaturePoly.c_str(), false));
		  //double restArmatureBounds[6] = { 0., -1., 0., -1., 0., -1. };
		  //armature->GetBounds(restArmatureBounds);
		  //std::cout << "Rest armature bounds: "
		//	  << restArmatureBounds[0] << ", " << restArmatureBounds[1] << ", "
			//  << restArmatureBounds[2] << ", " << restArmatureBounds[3] << ", "
			//  << restArmatureBounds[4] << ", " << restArmatureBounds[5] << std::endl;

 
		  Vec3 fixedA;
		  Vec3 fixedB;
		  Vec3 rotateA;
		  Vec3 rotateB;

		  vtkPoints* inPoints = armature->GetPoints();
		  vtkCellArray* armatureSegments = armature->GetLines();
		  vtkCellData* armatureCellData = armature->GetCellData();
		  vtkNew<vtkIdList> cell;
		  armatureSegments->InitTraversal();
		  int edgeId(0);
		  int i = 2;
		  while (armatureSegments->GetNextCell(cell.GetPointer()))
		  {
			  vtkIdType a = cell->GetId(0);
			  vtkIdType b = cell->GetId(1);

			  Vec3 ax(inPoints->GetPoint(a));
			  Vec3 bx(inPoints->GetPoint(b));

			  if (i == ComponentFixed) {
				  fixedA = ax;
				  fixedB = bx;
			  }
			  if (i == ComponentToRotate) {
				  rotateA = ax;
				  rotateB = bx;
			  }

			  //std::cout << "Segment " << i << "A : " << ax << " B : " << bx<< std::endl;
			  i++;
		  }

		  std::cout << "FixedPart " << ComponentFixed << "A : " << fixedA << " FixedPart B : " << fixedB << std::endl;
		  std::cout << "RotatePart " << ComponentToRotate << " A : " << rotateA << " RotatePart B : " << rotateB << std::endl;

		  Vec3 rotPart;
		  Vec3 fixedPart;
		  Vec3 fixedPoint;
		  vtkSmartPointer<vtkTransform> transformPointer = vtkSmartPointer<vtkTransform>::New();
  

		  if (fixedA == rotateB) {
			  rotPart = rotateA - fixedB;
			  fixedPart = fixedA - fixedB;
			  fixedPoint = fixedA;
		  }
		  else if (fixedB == rotateA) {
			  rotPart = rotateB - fixedA;
			  fixedPart = fixedB - fixedA;
			  fixedPoint = fixedB;
		  }

		  std::cout << "Two vecs: rot: " << rotPart << " fix: " << fixedPart << std::endl;

		  // Translate to origin
		  transformPointer->Translate(-fixedPoint.GetElement(0), -fixedPoint.GetElement(1), -fixedPoint.GetElement(2));
		 // transformPointer->PrintSelf(std::cout, *vtkIndent::New());

		  // Rotation around XY:
		  double va = atan2(rotPart.GetElement(0),rotPart.GetElement(1));
		  double ua = atan2(fixedPart.GetElement(0), fixedPart.GetElement(1));
		  double A = va - ua;
		  A = A - 360 * (A > 180) + 360 * (A < -180);
		  double thetaXY = -A;
		  transformPointer->RotateZ(thetaXY);
		//  transformPointer->PrintSelf(std::cout, *vtkIndent::New());

		  // Rotation around YZ:
		  va = atan2(rotPart.GetElement(1), rotPart.GetElement(2));
		  ua = atan2(fixedPart.GetElement(1), fixedPart.GetElement(2));
		  A = va - ua;
		  A = A - 360 * (A > 180) + 360 * (A < -180);
		  double thetaYZ = A;
		  transformPointer->RotateX(thetaYZ);
		  //transformPointer->PrintSelf(std::cout, *vtkIndent::New());

		  // Rotation around XZ:
		  va = atan2(rotPart.GetElement(0), rotPart.GetElement(2));
		  ua = atan2(fixedPart.GetElement(0), fixedPart.GetElement(2));
		  A = va - ua;
		  A = A - 360 * (A > 180) + 360 * (A < -180);
		  double thetaXZ = A;
		  transformPointer->RotateY(thetaXZ);
		 // transformPointer->PrintSelf(std::cout, *vtkIndent::New());

		  // Translate Back
		  transformPointer->Translate(fixedPoint.GetElement(0), fixedPoint.GetElement(1), fixedPoint.GetElement(2));
		  transformPointer->PrintSelf(std::cout, *vtkIndent::New());

		  // Apply Transformation to the Input Image
		  vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
		  transformFilter->SetTransform(transformPointer);
		  transformFilter->SetInputData((vtkDataObject*)reader->GetOutput());

		  // Write the applied transformation to the Output Image
		  vtkSmartPointer<vtkImageWriter> writer = vtkSmartPointer<vtkImageWriter>::New();
		  writer->SetFileName(outputVolume.c_str());
		  writer->SetInputData(transformFilter->GetOutput());
		  writer->Write();

		  std::cout << "Wrote output image, fname: " << outputVolume << " filter output:" << transformFilter->GetOutput() << std::endl;
		 
		  return EXIT_SUCCESS;
	}

} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  try
  {
	  itk::ImageIOBase::IOPixelType     pixelType;
	  itk::ImageIOBase::IOComponentType componentType;

	  itk::GetImageType(inputVolume, pixelType, componentType);

	  // This filter handles all types on input, but only produces
	  // signed types
	  switch (componentType)
	  {
	  case itk::ImageIOBase::UCHAR:
		  return DoIt<unsigned char>(argc, argv);
		  break;
	  case itk::ImageIOBase::CHAR:
		  return DoIt<char>(argc, argv);
		  break;
	  case itk::ImageIOBase::USHORT:
		  return DoIt<unsigned short>(argc, argv);
		  break;
	  case itk::ImageIOBase::SHORT:
		  return DoIt<short>(argc, argv);
		  break;
	  case itk::ImageIOBase::UINT:
		  return DoIt<unsigned int>(argc, argv);
		  break;
	  case itk::ImageIOBase::INT:
		  return DoIt<int>(argc, argv);
		  break;
	  case itk::ImageIOBase::ULONG:
		  return DoIt<unsigned long>(argc, argv);
		  break;
	  case itk::ImageIOBase::LONG:
		  return DoIt<long>(argc, argv);
		  break;
	  case itk::ImageIOBase::FLOAT:
		  return DoIt<float>(argc, argv);
		  break;
	  case itk::ImageIOBase::DOUBLE:
		  return DoIt<double>(argc, argv);
		  break;
	  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
	  default:
		  std::cerr << "Unknown component type: " << componentType << std::endl;
		  break;
	  }
  }

  catch (itk::ExceptionObject & excep)
  {
	  std::cerr << argv[0] << ": exception caught !" << std::endl;
	  std::cerr << excep << std::endl;
	  return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

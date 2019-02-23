#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "LogicCLP.h"

// bender includes
#include "benderIOUtils.h"

// ITK includes
#include <itkDirectory.h>

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>

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

  typedef itk::ImageFileWriter<OutputImageType> WriterType;

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
  vtkSmartPointer<vtkMatrix4x4> beforeTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  beforeTransformMatrix->Identity();
  vtkSmartPointer<vtkMatrix4x4> afterTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  afterTransformMatrix->Identity();

  if (fixedA == rotateB) {
	  rotPart = rotateA - fixedB;
	  fixedPart = fixedA - fixedB;

	  beforeTransformMatrix->SetElement(0, 3, -fixedA.GetElement(0));
	  beforeTransformMatrix->SetElement(1, 3, -fixedA.GetElement(1));
	  beforeTransformMatrix->SetElement(2, 3, -fixedA.GetElement(2));
	  afterTransformMatrix->SetElement(0, 3, fixedA.GetElement(0));
	  afterTransformMatrix->SetElement(1, 3, fixedA.GetElement(1));
	  afterTransformMatrix->SetElement(2, 3, fixedA.GetElement(2));
  }
  else if (fixedB == rotateA) {
	  rotPart = rotateB - fixedA;
	  fixedPart = fixedB - fixedA;

	  beforeTransformMatrix->SetElement(0,3,-fixedB.GetElement(0));
	  beforeTransformMatrix->SetElement(1,3,-fixedB.GetElement(1));
	  beforeTransformMatrix->SetElement(2,3,-fixedB.GetElement(2));
	  afterTransformMatrix->SetElement(0, 3, fixedB.GetElement(0));
	  afterTransformMatrix->SetElement(1, 3, fixedB.GetElement(1));
	  afterTransformMatrix->SetElement(2, 3, fixedB.GetElement(2));
  }

  std::cout << "Two vecs: rot: " << rotPart << " fix: " << fixedPart << std::endl;
  beforeTransformMatrix->PrintSelf(std::cout, *vtkIndent::New());
  afterTransformMatrix->PrintSelf(std::cout, *vtkIndent::New());


  vtkSmartPointer<vtkTransform> transfromBefore = vtkSmartPointer<vtkTransform>::New();
 
  transfrom->SetMatrix(beforeTransformMatrix);
  vtkSmartPointer<vtkTransform> transfromBefore = vtkSmartPointer<vtkTransform>::New();
  transfrom->SetMatrix(beforeTransformMatrix);

  typename WriterType::Pointer writer = WriterType::New();
  itk::PluginFilterWatcher watchWriter(writer,
	  "Write Transformation",
	  CLPProcessInformation);
  writer->SetFileName(beforeTransform.c_str());
  writer->SetInput(beforeTransformMatrix);
  writer->Update();

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

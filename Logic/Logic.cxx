#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "LogicCLP.h"

// bender includes
#include "benderIOUtils.h"
#include "benderWeightMap.h"
#include "benderWeightMapIO.h"
#include "benderWeightMapMath.h"

// ITK includes
#include <itkDirectory.h>

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


  vtkSmartPointer<vtkPolyData> armature;
  armature.TakeReference(bender::IOUtils::ReadPolyData(ArmaturePoly.c_str(), false));
  double restArmatureBounds[6] = { 0., -1., 0., -1., 0., -1. };
  armature->GetBounds(restArmatureBounds);
  std::cout << "Rest armature bounds: "
	  << restArmatureBounds[0] << ", " << restArmatureBounds[1] << ", "
	  << restArmatureBounds[2] << ", " << restArmatureBounds[3] << ", "
	  << restArmatureBounds[4] << ", " << restArmatureBounds[5] << std::endl;

 
  vtkPoints* inPoints = armature->GetPoints();
  vtkCellArray* armatureSegments = armature->GetLines();
  vtkCellData* armatureCellData = armature->GetCellData();
  vtkNew<vtkIdList> cell;
  armatureSegments->InitTraversal();
  int edgeId(0);
  int i = 2;
  
  Vec3 fixedA = new Vec3();
  Vec3 fixedB = new Vec3();
  Vec3 rotateA = new Vec3();
  Vec3 rotateB = new Vec3();

  while (armatureSegments->GetNextCell(cell.GetPointer()))
  {
	  vtkIdType a = cell->GetId(0);
	  vtkIdType b = cell->GetId(1);

	  Vec3 ax(inPoints->GetPoint(a));
	  Vec3 bx(inPoints->GetPoint(b));

	  if (i  == ComponentFixed) {

	  }
	  if()
  }

  ComponentFixed->

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

void GetRotationalPartAndFixedPart(vtkSmartPointer<vtkPolyData> armature, int fixed, int rotation, Vec3* fixedA, Vec3* fixedB, Vec3* rotateA, Vec3* rotateB) {
}

//----------------------------------------------------------------------------
void GetWeightFileNames(const std::string& dirName, std::vector<std::string>& fnames)
{
	fnames.clear();
	itk::Directory::Pointer dir = itk::Directory::New();
	dir->Load(dirName.c_str());
	for (unsigned int i = 0; i < dir->GetNumberOfFiles(); ++i)
	{
		std::string name = dir->GetFile(i);
		if (strstr(name.c_str(), ".mha"))
		{
			std::string fname = dirName;
			fname.append("/");
			fname.append(name);
			fnames.push_back(fname);
		}
	}

	std::sort(fnames.begin(), fnames.end());
}

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

#define _USE_MATH_DEFINES

#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "LogicCLP.h"

#include <math.h>

// bender includes
#include "benderIOUtils.h"

// ITK includes
#include <itkDirectory.h>
#include <itkImageTransformer.h>
#include <itkAffineTransform.h>
#include <itkCenteredEuler3DTransform.h>
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
	int DoIt(int argc, char * argv[])
	{
		PARSE_ARGS;


		typedef    T InputPixelType;
		typedef    T OutputPixelType;

		typedef itk::Image<InputPixelType, 3> InputImageType;
		typedef itk::Image<OutputPixelType, 3> OutputImageType;

		typedef itk::ImageFileReader<InputImageType>  ReaderType;
		typedef itk::ImageFileWriter<OutputImageType> WriterType;

		typedef itk::ResampleImageFilter<InputImageType, OutputImageType>                                    ResampleType;

		typename ReaderType::Pointer reader = ReaderType::New();
		itk::PluginFilterWatcher watchReader(reader, "Read Volume",
			CLPProcessInformation);

		using ScalarType = double;
		constexpr unsigned int Dimension = 3;

		reader->SetFileName(inputVolume.c_str());

		std::cout << "Read input image" << std::endl;
		reader->Update();
		std::cout << "Input image read" << std::endl;
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

		Vec3 bellyA;
		Vec3 bellyB;

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

			if (i == 6) {
				// Get the belly for the legs
				bellyA = ax;
				bellyB = bx;
			}

			//std::cout << "Segment " << i << " A : " << ax << " B : " << bx<< std::endl;
			i++;
		}

		Vec3 rotPart;
		Vec3 fixedPart;
		Vec3 fixedPoint;
		vtkSmartPointer<vtkTransform> transformPointer = vtkSmartPointer<vtkTransform>::New();

		if (fixedA == rotateB) {
			rotPart = rotateA - rotateB;
			fixedPart = fixedA - fixedB;
			fixedPoint = fixedA;
		}
		else if (fixedB == rotateA) {
			rotPart = rotateB - rotateA;
			fixedPart = fixedB - fixedA;
			fixedPoint = fixedB;
		}

		//	std::cout << "FixedPart " << ComponentFixed << "A : " << fixedA << " FixedPart B : " << fixedB << std::endl;
		//	std::cout << "RotatePart " << ComponentToRotate << " A : " << rotateA << " RotatePart B : " << rotateB << std::endl;


		if (ComponentToRotate == 10 || ComponentToRotate == 13) {
			// The legs have to be rotated agains the belly NOT the pelvis
			fixedPart = bellyA - bellyB;
		}

		// Normalize the vectors
		rotPart = rotPart / rotPart.GetNorm();
		fixedPart = fixedPart / fixedPart.GetNorm();

		//std::cout << "Two vecs: rot: " << rotPart << " fix: " << fixedPart << std::endl;

		/*
		  Translation and rotaion of the Rotational part according to the given fixed Part
		  done by using an affine transformation.
		*/
		using TransformType = itk::CenteredEuler3DTransform<ScalarType>;
		TransformType::Pointer eulerTransform = TransformType::New();

		const InputImageType::SizeType& size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
		using ResampleImageFilterType = itk::ResampleImageFilter< InputImageType, OutputImageType >;

		/*
		  ========
		  ROTATION BEGIN
		  =========
		*/
		// Rotation around XY:
		double thetaXY = getTheta(rotPart.GetElement(0), rotPart.GetElement(1),
			fixedPart.GetElement(0), fixedPart.GetElement(1));

		//std::cout << "thetaXY: " << thetaXY << std::endl;
		double axisZ[3] = { 0.0, 0.0, 1.0 };
		//  transformRotate->Rotate3D(Vec3(axisZ), thetaXY, false);

		  // calculate the new coordinates of the rotational part in order to carry on with rotational calculation
		double arrayRz[3][3] = { {cos(thetaXY), sin(thetaXY), 0},
							  {-sin(thetaXY),cos(thetaXY),0},
							  { 0,0,1} };
		Mat33 Rz = ToItkMatrix(arrayRz);
		rotPart = Rz * rotPart;
		//	std::cout << "rotPart: " << rotPart << std::endl;

			// Rotation around YZ:
		double thetaYZ = getTheta(rotPart.GetElement(1), rotPart.GetElement(2),
			fixedPart.GetElement(1), fixedPart.GetElement(2));
		//	std::cout << "thetaYZ: " << thetaYZ << std::endl;
		double axisX[3] = { 1.0, 0.0, 0.0 };

		// calculate the new coordinates of the rotational part in order to carry on with rotational calculation
		double arrayRx[3][3] = { {1, 0, 0},
							  {0,cos(thetaYZ),sin(thetaYZ)},
							  {0,-sin(thetaYZ),cos(thetaYZ)} };
		Mat33 Rx = ToItkMatrix(arrayRx);
		rotPart = Rx * rotPart;
		//std::cout << "rotPart: " << rotPart << std::endl;

		// Rotation around XZ:
		double thetaXZ = -getTheta(rotPart.GetElement(0), rotPart.GetElement(2),
			fixedPart.GetElement(0), fixedPart.GetElement(2));
		//	std::cout << "thetaXZ: " << thetaXZ << std::endl;
		double axisY[3] = { 0.0, 1.0, 0.0 };

		/*
		  ========
		  ROTATION END
		  =========
		*/
		double parametersArray[9] = { -thetaYZ, thetaXZ, -thetaXY, abs(fixedPoint[0]),abs(fixedPoint[1]), abs(fixedPoint[2]), 0,0,0 };
		
		std::cout << "BEFORE ADAPTION:" <<endl<< "FixedPart " << ComponentFixed << " vec: " << fixedPart <<
			" RotatePart " << ComponentToRotate << " vec: " << rotPart << endl <<
			"Rx: " << parametersArray[0] << " Ry: " << parametersArray[1] << " Rz: " << parametersArray[2] << endl <<
			"Cross Rx: " << getCross(rotPart.GetElement(1), rotPart.GetElement(2),
				fixedPart.GetElement(1), fixedPart.GetElement(2)) <<
			" Cross Ry: " << getCross(rotPart.GetElement(0), rotPart.GetElement(2),
				fixedPart.GetElement(0), fixedPart.GetElement(2)) <<
			" Cross Rz: " << getCross(rotPart.GetElement(0), rotPart.GetElement(1),
				fixedPart.GetElement(0), fixedPart.GetElement(1)) << endl <<
			"Cx: " << parametersArray[3] <<
			" Cy: " << parametersArray[4] <<
			" Cz: " << parametersArray[5] << std::endl;

		/*
		====================================================
		Some adaptions that have to be done for the skeleton
		====================================================
		*/
		specialAdaptions(parametersArray, ComponentToRotate);

		TransformType::ParametersType params(9);
		for (int i = 0; i < 9; i++) {
			params[i] = parametersArray[i];
		}
		eulerTransform->SetParameters(params);

		std::cout << "AFTER ADAPTION:" << endl << "FixedPart " << ComponentFixed << " vec: " << fixedPart <<
			" RotatePart " << ComponentToRotate << " vec: " << rotPart << endl <<
			"Rx: " << parametersArray[0] << " Ry: " << parametersArray[1] << " Rz: " << parametersArray[2] << endl <<
			"Cross Rx: " << getCross(rotPart.GetElement(1), rotPart.GetElement(2),
				fixedPart.GetElement(1), fixedPart.GetElement(2)) <<
			" Cross Ry: " << getCross(rotPart.GetElement(0), rotPart.GetElement(2),
				fixedPart.GetElement(0), fixedPart.GetElement(2)) <<
			" Cross Rz: " << getCross(rotPart.GetElement(0), rotPart.GetElement(1),
				fixedPart.GetElement(0), fixedPart.GetElement(1)) << endl <<
			"Cx: " << parametersArray[3] <<
			" Cy: " << parametersArray[4] <<
			" Cz: " << parametersArray[5] << std::endl;

		//eulerTransform->SetCenter(fixedPoint);
	//	std::cout << "Euler Transform with Center Set to fixed Point";
	//	eulerTransform->Print(std::cout);

		/*
		  Translation and rotation finished
		*/

		/*
		  Result Resampling and writing
		*/
		typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
		resample->SetInput(reader->GetOutput());
		resample->SetTransform(eulerTransform);
		resample->SetReferenceImage(reader->GetOutput());
		resample->UseReferenceImageOn();
		resample->SetSize(size);

		typename WriterType::Pointer writer = WriterType::New();
		itk::PluginFilterWatcher watchWriter(writer,
			"Write Volume",
			CLPProcessInformation);
		writer->SetFileName(outputVolume.c_str());
		writer->SetInput(resample->GetOutput());
		writer->SetUseCompression(1);
		writer->Update();


		return EXIT_SUCCESS;
	}

} // end of anonymous namespace

double get180adaption(double angle) {
	if (angle < 0) {
		return M_PI + angle;
	}
	return M_PI - angle;
}

void specialAdaptions(double* parametersArray, int ComponentToRotate) {

	if (ComponentToRotate == 3) { // Nack
		if (parametersArray[1] > 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 21) { // LEFT HAND
		if (parametersArray[1] > 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 20) { // LEFT UNDER ARM
		if (parametersArray[1] > 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 19) { // LEFT UPPER ARM
		if (parametersArray[0] < 0) {
			parametersArray[0] = -parametersArray[0];
		}
		if (parametersArray[1] > 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 14) { // LEFT FEMURE
		parametersArray[2] = get180adaption(parametersArray[2]);
		if (parametersArray[0] < 0) {
			parametersArray[0] = -parametersArray[0];
		}
	}
	else if (ComponentToRotate == 13) { // LEFT Tibea
		parametersArray[0] = get180adaption(parametersArray[0]);
		if (parametersArray[0] > 0) {
			parametersArray[0] = -parametersArray[0];
		}
		if (parametersArray[1] < 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 17) { // RIGHT UNDER ARM
		if (parametersArray[1] > 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 16) { // RIGHT UPPER ARM
		if (parametersArray[0] > 0) {
		}
		if (parametersArray[1] < 0) {
			parametersArray[1] = -parametersArray[1];
		}
		if (parametersArray[2] < 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 10) { // RIGHT TIBEA
		parametersArray[0] = get180adaption(parametersArray[0]);
		if (parametersArray[0] > 0) {
			parametersArray[0] = -parametersArray[0];
		}
		if (parametersArray[2] > 0) {
			parametersArray[2] = -parametersArray[2];
		}
	}
	else if (ComponentToRotate == 11) { // RIGHT FEMURE
		parametersArray[2] = get180adaption(parametersArray[2]);
		if (parametersArray[0] < 0) {
			parametersArray[0] = -parametersArray[0];
		}
	}
}

double getCross(double vx, double vy, double ux, double uy) {
	return vx * uy - vy * ux;
}

double getTheta(double vx, double vy, double ux, double uy) {
	double va = -atan2(vx, vy)*180.0 / M_PI;
	double ua = -atan2(ux, uy)*180.0 / M_PI;
	double A = (va - ua);
	A = A - 360 * (A > 180) + 360 * (A < -180);
	//	std::cout << "vx: " << vx << " vy: " << vy << " va : " << va
	//		<< " ux: " << ux << " uy: " << uy << " ua: " << ua << std::endl;
	return A * M_PI / 180;
}

Mat33 ToItkMatrix(double M[3][3])
{
	Mat33 itkM;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			itkM(i, j) = M[i][j];
		}
	}

	return itkM;
}

int main(int argc, char * argv[])
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

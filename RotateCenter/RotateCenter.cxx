#define _USE_MATH_DEFINES

#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include <math.h>

#include "RotateCenterCLP.h"

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

		double parametersArray[9] = { Rx, Ry, Rz, Cx,Cy,Cz, 0,0,0 };
		if (NegateRx) {
			parametersArray[0] = -parametersArray[0];
		}
		if (NegateRy) {
			parametersArray[1] = -parametersArray[1];
		}
		if (NegateRz) {
			parametersArray[2] = -parametersArray[2];
		}
		if (NegateCx) {
			parametersArray[3] = -parametersArray[3];
		}
		if (NegateCy) {
			parametersArray[4] = -parametersArray[4];
		}
		if (NegateCz) {
			parametersArray[5] = -parametersArray[5];
		}

		TransformType::ParametersType params(9);
		for (int i = 0; i < 9; i++) {
			params[i] = parametersArray[i];
		}
		eulerTransform->SetParameters(params);
		//eulerTransform->SetCenter(fixedPoint);
		std::cout << "Euler Transform with Center Set to fixed Point";
		eulerTransform->Print(std::cout);

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

#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "ExpandCLP.h"

// ITK
#include <itkPasteImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkCenteredEuler3DTransform.h>

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

	template <typename TPixel>
	int DoIt(int argc, char * argv[], TPixel)
	{
		PARSE_ARGS;

		typedef TPixel InputPixelType;
		typedef TPixel OutputPixelType;

		const unsigned int Dimension = 3;

		typedef itk::Image<InputPixelType, Dimension> InputImageType;
		typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

		typedef itk::ImageFileReader<InputImageType>  ReaderType;

		typename ReaderType::Pointer reader = ReaderType::New();

		reader->SetFileName(inputVolume.c_str());
		reader->Update();

		itk::PasteImageFilter< InputImageType, InputImageType, InputImageType >::Pointer paster =
			itk::PasteImageFilter< InputImageType, InputImageType, InputImageType >::New();

		itk::ResampleImageFilter< InputImageType, InputImageType >::Pointer transformer =
			itk::ResampleImageFilter< InputImageType, InputImageType >::New();

		itk::CenteredEuler3DTransform< double >::Pointer transform =
			itk::CenteredEuler3DTransform< double >::New();

		typename InputImageType::Pointer inputImage = reader->GetOutput();

		InputImageType::PointType origin = inputImage->GetOrigin();
		InputImageType::SpacingType spacing = inputImage->GetSpacing();
		InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();

		// Make a bounding box large enough for the potentally largest rotated volume
		InputImageType::SizeType outputsize;
		outputsize[0] = size[0] * 1.6;
		outputsize[1] = size[1] * 1.6;
		outputsize[2] = size[2] * 1.6;
		//outputsize[1] = sqrt(pow((double)size[1] * spacing[1], 2.0) +
	//		pow((double)size[2] * spacing[2], 2.0)) / spacing[1];
	//	outputsize[2] = (double)outputsize[1] * spacing[1] / spacing[2];


		// Calculate the centre of the volume
		InputImageType::PointType centre;
		centre[0] = origin[0] + spacing[0] * (double)(outputsize[0] - 1) / 2.0;
		centre[1] = origin[1] + spacing[1] * (double)(outputsize[1] - 1) / 2.0;
		centre[2] = origin[2] + spacing[2] * (double)(outputsize[2] - 1) / 2.0;

		// Shift the origin over to get the whole image in
		origin[0] -= spacing[0] * static_cast<double>(outputsize[0] - size[0]) / 2.0;
		origin[1] -= spacing[1] * static_cast<double>(outputsize[1] - size[1]) / 2.0;
		origin[2] -= spacing[2] * static_cast<double>(outputsize[2] - size[2]) / 2.0;

		// Create a blank image volume to extend the bounds
		InputImageType::Pointer boundingimg = InputImageType::New();
		boundingimg->SetSpacing(spacing);
		//boundingimg->SetOrigin(origin);
		InputImageType::IndexType index;
		index[0] = 0;
		index[1] = 0;
		index[2] = 0;
		InputImageType::RegionType region;
		region.SetSize(outputsize);
		region.SetIndex(index);
		boundingimg->SetRegions(region);
		boundingimg->Allocate();
		boundingimg->FillBuffer(0);
		boundingimg->Update();

		// Paste the image in the middle of the new blank image
		index[0] = (outputsize[0] - size[0]) / 2 + origin[0];
		index[1] = (outputsize[1] - size[1]) / 2 + origin[1];
		index[2] = (outputsize[2] - size[2]) / 2 + origin[2];
		paster->SetDestinationIndex(index);
		paster->SetDestinationImage(boundingimg);
		paster->SetSourceImage(inputImage);
		paster->SetSourceRegion(inputImage->GetLargestPossibleRegion());
		paster->Update();

		typedef itk::ImageFileWriter<OutputImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(outputVolume.c_str());
		writer->SetInput(paster->GetOutput());
		writer->SetUseCompression(1);
		writer->Update();

		return EXIT_SUCCESS;
	}

} // end of anonymous namespace

int main(int argc, char * argv[])
{
	PARSE_ARGS;

	itk::ImageIOBase::IOPixelType     pixelType;
	itk::ImageIOBase::IOComponentType componentType;

	try
	{
		itk::GetImageType(inputVolume, pixelType, componentType);

		// This filter handles all types on input, but only produces
		// signed types
		switch (componentType)
		{
		case itk::ImageIOBase::UCHAR:
			return DoIt(argc, argv, static_cast<unsigned char>(0));
			break;
		case itk::ImageIOBase::CHAR:
			return DoIt(argc, argv, static_cast<signed char>(0));
			break;
		case itk::ImageIOBase::USHORT:
			return DoIt(argc, argv, static_cast<unsigned short>(0));
			break;
		case itk::ImageIOBase::SHORT:
			return DoIt(argc, argv, static_cast<short>(0));
			break;
		case itk::ImageIOBase::UINT:
			return DoIt(argc, argv, static_cast<unsigned int>(0));
			break;
		case itk::ImageIOBase::INT:
			return DoIt(argc, argv, static_cast<int>(0));
			break;
		case itk::ImageIOBase::ULONG:
			return DoIt(argc, argv, static_cast<unsigned long>(0));
			break;
		case itk::ImageIOBase::LONG:
			return DoIt(argc, argv, static_cast<long>(0));
			break;
		case itk::ImageIOBase::FLOAT:
			return DoIt(argc, argv, static_cast<float>(0));
			break;
		case itk::ImageIOBase::DOUBLE:
			return DoIt(argc, argv, static_cast<double>(0));
			break;
		case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
		default:
			std::cerr << "Unknown input image pixel component type: ";
			std::cerr << itk::ImageIOBase::GetComponentTypeAsString(componentType);
			std::cerr << std::endl;
			return EXIT_FAILURE;
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

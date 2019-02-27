#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "AddVolumeCLP.h"

// ITK includes
#include <itkDirectory.h>
#include <itkConstrainedValueAdditionImageFilter.h>
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

  typedef itk::ResampleImageFilter<InputImageType, OutputImageType> ResampleType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  typename ReaderType::Pointer reader2 = ReaderType::New();

  typedef itk::ConstrainedValueAdditionImageFilter<InputImageType, OutputImageType, OutputImageType> FilterType;

  using ScalarType = double;
  constexpr unsigned int Dimension = 3;

  reader1->SetFileName(inputVolume1.c_str());
  reader1->Update();
  reader2->SetFileName(inputVolume2.c_str());
  reader2->Update();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(reader1->GetOutput());
  filter->SetInput2(reader2->GetOutput());

  /*
	Result Resampling and writing
  */
  const InputImageType::SizeType& size = reader1->GetOutput()->GetLargestPossibleRegion().GetSize();

  typename ResampleType::Pointer resample = ResampleType::New();
  resample->SetInput(filter->GetOutput());
  resample->SetReferenceImage(filter->GetOutput());
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

		itk::GetImageType(inputVolume1, pixelType, componentType);

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

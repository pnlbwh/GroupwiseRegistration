#include "WarpVolumeCLP.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include "itkWarpImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

std::string WarpedImageName(std::string outputDir, std::string filename)
{
  std::stringstream result;
  result << outputDir << "/" << itksys::SystemTools::GetFilenameWithoutExtension(filename) << "_warped.nrrd";
  return result.str();
}


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::cout << "output directory:" << resultsDirectory << std::endl;
  std::cout << "warp:" << warp << std::endl;
  std::cout << "input volume:" << input << std::endl;


  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Vector<float, Dimension>  VectorPixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::Image<VectorPixelType, Dimension>  DeformationFieldType;
  typedef itk::WarpImageFilter <ImageType, ImageType, DeformationFieldType>  WarperType;
  typedef itk::ImageFileReader< DeformationFieldType >    DeformationReaderType;
  typedef itk::ImageFileReader< ImageType >                     ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >                     WriterType;

  DeformationReaderType::Pointer   fieldReader = DeformationReaderType::New();
  fieldReader->SetFileName( warp );

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( input );
  imageReader->Update();

  WarperType::Pointer   warper = WarperType::New();
  warper->SetDeformationField( fieldReader->GetOutput() );
  warper->SetInput( imageReader->GetOutput() );
  warper->SetOutputSpacing( imageReader->GetOutput()->GetSpacing() );
  warper->SetOutputOrigin( imageReader->GetOutput()->GetOrigin() );
  warper->SetOutputDirection( imageReader->GetOutput()->GetDirection() );

  WriterType::Pointer  writer =  WriterType::New();
  writer->SetFileName( WarpedImageName(resultsDirectory, input) );
  writer->SetInput( warper->GetOutput() );
  //writer->SetInput( imageReader->GetOutput() );
  writer->SetUseCompression( true );
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Could not write warped image" << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }

  return EXIT_SUCCESS;
}


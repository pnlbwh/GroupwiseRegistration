#include "WarpVolumeCLP.h"
#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkImageIOBase.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <vector>
#include "itkWarpImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


struct parameters
{
  std::string inputVolume;
  std::string warp;
  std::string resultsDirectory;
};



/* 
 * This function was taken from ResampleVolume2.cxx
 *
 * Separate the vector image into a vector of images
 */
template< class PixelType >
int SeparateImages( const typename itk::VectorImage< PixelType , 3 >
                    ::Pointer &imagePile ,
                    std::vector< typename itk::OrientedImage< PixelType , 3 >::Pointer > &vectorImage
                  )
{
   typedef itk::OrientedImage< PixelType , 3 > ImageType;
   typedef itk::VectorImage< PixelType , 3 > VectorImageType;
   typename itk::VectorImage< PixelType , 3 >::SizeType size;
   typename itk::VectorImage< PixelType , 3 >::DirectionType direction;
   typename itk::VectorImage< PixelType , 3 >::PointType origin;
   typename itk::VectorImage< PixelType , 3 >::SpacingType spacing;
   size = imagePile->GetLargestPossibleRegion().GetSize();
   direction=imagePile->GetDirection();
   origin=imagePile->GetOrigin();
   spacing=imagePile->GetSpacing();
   typename itk::ImageRegionIterator< VectorImageType > in( imagePile ,
                                                            imagePile->GetLargestPossibleRegion() );
   typedef typename itk::ImageRegionIterator< ImageType > IteratorImageType;
   std::vector< IteratorImageType > out;
   for( unsigned int i = 0; i < imagePile->GetVectorLength() ; i++ )
   {
      typename ImageType::Pointer imageTemp = ImageType::New();
      imageTemp->SetRegions( size );
      imageTemp->SetOrigin( origin );
      imageTemp->SetDirection( direction );
      imageTemp->SetSpacing( spacing );
      imageTemp->Allocate();
      vectorImage.push_back( imageTemp );
      IteratorImageType outtemp( imageTemp , imageTemp->GetLargestPossibleRegion() );
      outtemp.GoToBegin();
      out.push_back( outtemp );
   }
   for( in.GoToBegin(); !in.IsAtEnd() ; ++in )
   {
      itk::VariableLengthVector< PixelType > value = in.Get();
      for( unsigned int i = 0; i < imagePile->GetVectorLength() ; i++ )
      {
         out[ i ].Set( value[ i ] );
         ++out[ i ];
      }
   }
   return EXIT_SUCCESS;
}

/*
 * This function was taken from ResampleVolume2.cxx
 *
 * Write back the vector of images into a image vector
 */
template<class PixelType>
int AddImage( typename itk::VectorImage< PixelType, 3 >
              ::Pointer &imagePile,
              const std::vector< typename itk::OrientedImage< PixelType , 3 > ::Pointer > &vectorImage
            )
{
   typedef itk::OrientedImage< PixelType , 3 > ImageType;
   imagePile->SetRegions( vectorImage.at( 0 )->GetLargestPossibleRegion().GetSize() );
   imagePile->SetOrigin( vectorImage.at( 0 )->GetOrigin() );
   imagePile->SetDirection( vectorImage.at( 0 )->GetDirection() );
   imagePile->SetSpacing( vectorImage.at( 0 )->GetSpacing() );
   imagePile->SetVectorLength( vectorImage.size() );
   imagePile->Allocate();
   typename itk::ImageRegionIterator< itk::VectorImage< PixelType , 3 > > out( imagePile ,
                                                                               imagePile->GetLargestPossibleRegion()
                                                                             );
   typedef typename itk::ImageRegionIterator< ImageType > IteratorImageType;
   std::vector< IteratorImageType > in;
   for( unsigned int i = 0; i < imagePile->GetVectorLength() ; i++ )
   {
      IteratorImageType intemp( vectorImage.at( i ) , vectorImage.at( i )->GetLargestPossibleRegion() );
      intemp.GoToBegin();
      in.push_back( intemp );
   }
   itk::VariableLengthVector< PixelType > value;
   value.SetSize( vectorImage.size() );
   for( out.GoToBegin(); !out.IsAtEnd() ; ++out )
   {
      for( unsigned int i = 0; i < imagePile->GetVectorLength() ; i++ )
      {
         value.SetElement( i , in.at( i ).Get() );
         ++in[ i ];
      }
    out.Set( value );
    }
  return EXIT_SUCCESS;
}

std::string WarpedImageName(std::string outputDir, std::string filename)
{
  std::stringstream result;
  result << outputDir << "/" << itksys::SystemTools::GetFilenameWithoutExtension(filename) << "_warped.nrrd";
  return result.str();
}

void GetImageType( std::string fileName ,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType
                 )
{
   typedef itk::Image< unsigned char , 3 > ImageType;
   itk::ImageFileReader< ImageType >::Pointer imageReader;
   imageReader = itk::ImageFileReader< ImageType >::New();
   imageReader->SetFileName( fileName.c_str() );
   imageReader->UpdateOutputInformation();
   pixelType = imageReader->GetImageIO()->GetPixelType();
   componentType = imageReader->GetImageIO()->GetComponentType();
}



template< class PixelType > 
int Warp( parameters &list )
{
  const unsigned int Dimension = 3;
  typedef itk::Vector<float, Dimension>  VectorPixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::Image<VectorPixelType, Dimension>  DeformationFieldType;
  typedef itk::WarpImageFilter <ImageType, ImageType, DeformationFieldType>  WarperType;
  typedef itk::ImageFileReader< DeformationFieldType >    DeformationReaderType;
  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >   WriterType;

  DeformationReaderType::Pointer   fieldReader = DeformationReaderType::New();
  fieldReader->SetFileName( list.warp.c_str() );

  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( list.inputVolume.c_str() );
  imageReader->Update();

  typename WarperType::Pointer   warper = WarperType::New();
  warper->SetDeformationField( fieldReader->GetOutput() );
  warper->SetInput( imageReader->GetOutput() );
  warper->SetOutputSpacing( imageReader->GetOutput()->GetSpacing() );
  warper->SetOutputOrigin( imageReader->GetOutput()->GetOrigin() );
  warper->SetOutputDirection( imageReader->GetOutput()->GetDirection() );

  typename WriterType::Pointer  writer =  WriterType::New();
  writer->SetFileName( WarpedImageName(list.resultsDirectory, list.inputVolume) );
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
}


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::cout << "output directory:" << resultsDirectory << std::endl;
  std::cout << "warp:" << warp << std::endl;
  std::cout << "input volume:" << inputVolume << std::endl;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  GetImageType( inputVolume , pixelType , componentType );

  parameters list;
  list.resultsDirectory = resultsDirectory;
  list.warp = warp;
  list.inputVolume = inputVolume;

  switch( componentType )
   {
      case itk::ImageIOBase::UCHAR:
         return Warp< unsigned char >( list );
         break;
      case itk::ImageIOBase::CHAR:
         return Warp< char >( list );
         break;
      case itk::ImageIOBase::USHORT:
         return Warp< unsigned short >( list );
         break;
      case itk::ImageIOBase::SHORT:
         return Warp< short >( list );
         break;
      case itk::ImageIOBase::UINT:
         return Warp< unsigned int >( list );
         break;
      case itk::ImageIOBase::INT:
         return Warp< int >( list );
         break;
      case itk::ImageIOBase::ULONG:
         return Warp< unsigned long >( list );
         break;
      case itk::ImageIOBase::LONG:
         return Warp< long >( list );
         break;
      case itk::ImageIOBase::FLOAT:
         return Warp< float >( list );
         break;
      case itk::ImageIOBase::DOUBLE:
         return Warp< double >( list );
         break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
         std::cerr << "Unknown component type" << std::endl;
         break;
   }
   return EXIT_FAILURE;
}
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <itkMetaDataObject.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkVectorImage.h>

#include "itkPluginFilterWatcher.h"
#include "itkPluginUtilities.h"
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include "itkAddImageFilter.h"

#include "itkOrientedImage.h"
#include "itkOrientImageFilter.h"
#include "itkMultiResolutionPDEDeformableRegistration2.h"
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include "itkCommand.h"
#include "itkWarpJacobianDeterminantFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkWarpSmoothnessCalculator.h"
#include "itkWarpJacobianDeterminantFilter.h"
#include "itkGridForwardWarpImageFilter.h"
#include "itkVectorCentralDifferenceImageFunction.h"

#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include "itkMultiplyByConstantImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "GroupWiseRegistrationCLP.h"

#include "itkExponentialDeformationFieldImageFilter.h"
#include "itkESMInvConDemonsRegistrationFunction.h"
#include <itkNeighborhoodAlgorithm.h>
#include "itkSmoothingRecursiveGaussianImageFilter.h"

namespace 
{
  const unsigned int Dimension = 3;
  typedef float PixelType;

  typedef itk::Vector<float, Dimension>  VectorPixelType;

  typedef itk::Image< PixelType, Dimension >                    ImageType;
  typedef itk::ImageFileReader< ImageType >                     ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >                     WriterType;
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType>  AdderType;
  typedef itk::MultiplyByConstantImageFilter< ImageType, float, ImageType >   MultiplyByConstantType;
  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType>             DivideByImageType;

  typedef itk::Image<VectorPixelType, Dimension>          DeformationFieldType;
  typedef itk::ImageFileReader< DeformationFieldType >    DeformationReaderType;
  typedef itk::ImageFileWriter< DeformationFieldType >    DeformationWriterType;
  typedef itk::AddImageFilter<DeformationFieldType, DeformationFieldType, DeformationFieldType>     DeformationAdderType;
  typedef itk::MultiplyByConstantImageFilter< DeformationFieldType, float, DeformationFieldType >   DeformationMultiplyByConstantType;

  typedef itk::WarpImageFilter <ImageType, ImageType, DeformationFieldType>     WarperType;

  typedef itk::ESMInvConDemonsRegistrationFunction <ImageType, ImageType, DeformationFieldType>  DemonsRegistrationFunctionType;
  typedef itk::DiffeomorphicDemonsRegistrationFilter < ImageType, ImageType, DeformationFieldType>   ActualRegistrationFilterType;
  typedef itk::MultiResolutionPDEDeformableRegistration2< ImageType, ImageType, DeformationFieldType, PixelType >   MultiResRegistrationFilterType;


struct arguments
{
   std::string outputImageFile;
   unsigned int numLevels;
   unsigned int numIterations;
   std::vector<std::string> volumeFileNames;   
   std::string resultsDirectory;
   unsigned int numOuterIterations;
   float initialSigmaDiff;
   float finalSigmaDiff;
   float regWeight;
   bool useJac;
   float sigmaDef;        
   float sigmaUp;          
   float maxStepLength;     
   unsigned int verbosity;
   // unsigned int gradientType;

   // std::string  fixedImageFile;  /* -f option */
   // std::string  movingImageFile; /* -m option */
   // std::string  outputFieldFile; /* -O option */
   // std::string  trueFieldFile;   /* -r option */
   // std::vector<unsigned int> numIterations;   /* -i option */
   // bool useVanillaDem;           /* -a option */
   // bool useHistogramMatching;    /* -e option */

   arguments () :
     numLevels(6u),
     numIterations(50u),
     resultsDirectory("./"),
     numOuterIterations(50),
     initialSigmaDiff(10),
     finalSigmaDiff(2),
     regWeight(1e1),
     useJac(true),
     outputImageFile("output.nii.gz"),
     sigmaDef(3.0f),
     sigmaUp(0.0f),
     maxStepLength(2.0f),
     verbosity(true)
   {
     volumeFileNames = std::vector<std::string>(1, "");
   }

   // friend std::ostream& operator<< (std::ostream& o, const arguments& args)
   // {
    // std::ostringstream osstr;
    // for (unsigned int i=0; i<args.numIterations.size(); ++i)
    //    osstr<<args.numIterations[i]<<" ";
    // std::string iterstr = "[ " + osstr.str() + "]";

    // std::string gtypeStr;
    // switch (args.gradientType)
    // {
    // case 0:
    //    gtypeStr = "symmetrized";
    //    break;
    // case 1:
    //    gtypeStr = "fixed image";
    //    break;
    // case 2:
    //    gtypeStr = "moving image";
    //    break;
    // default:
    //    gtypeStr = "unsuported";
    // }
       
    // return o
    //    <<"Arguments structure:"<<std::endl
    //    <<"  Fixed image file: "<<args.fixedImageFile<<std::endl
    //    <<"  Moving image file: "<<args.movingImageFile<<std::endl
    //    <<"  Output image file: "<<args.outputImageFile<<std::endl
    //    <<"  Output field file: "<<args.outputFieldFile<<std::endl
    //    <<"  True field file: "<<args.trueFieldFile<<std::endl
    //    <<"  Number of multiresolution levels: "<<args.numLevels<<std::endl
    //    <<"  Number of demons iterations: "<<args.numIterations<<std::endl //<<"  Number of demons iterations: "<<iterstr<<std::endl
    //    <<"  Deformation field sigma: "<<args.sigmaDef<<std::endl
    //    <<"  Update field sigma: "<<args.sigmaUp<<std::endl
    //    <<"  Maximum update step length: "<<args.maxStepLength<<std::endl
    //    <<"  Use vanilla demons: "<<args.useVanillaDem<<std::endl
    //    <<"  Type of gradient: "<<gtypeStr<<std::endl
    //    <<"  Use histogram matching: "<<args.useHistogramMatching<<std::endl
    //    <<"  Verbosity: "<<args.verbosity;
   // }
};

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
template <class TPixel=float, unsigned int VImageDimension=3>
class CommandIterationUpdate : public itk::Command 
{
public:
   typedef  CommandIterationUpdate   Self;
   typedef  itk::Command             Superclass;
   typedef  itk::SmartPointer<Self>  Pointer;

   typedef itk::Image< TPixel, VImageDimension > InternalImageType;
   typedef itk::Vector< TPixel, VImageDimension >    VectorPixelType;
   typedef itk::Image<  VectorPixelType, VImageDimension > DeformationFieldType;

   typedef itk::DiffeomorphicDemonsRegistrationFilter<
      InternalImageType,
      InternalImageType,
      DeformationFieldType>   DiffeomorphicDemonsRegistrationFilterType;

   // typedef itk::FastSymmetricForcesDemonsRegistrationFilter<
   //    InternalImageType,
   //    InternalImageType,
   //    DeformationFieldType>   FastSymmetricForcesDemonsRegistrationFilterType;

   typedef itk::MultiResolutionPDEDeformableRegistration2<
      InternalImageType, InternalImageType,
      DeformationFieldType, TPixel >   MultiResRegistrationFilterType;

   typedef itk::WarpJacobianDeterminantFilter<DeformationFieldType, InternalImageType> JacobianFilterType;
   
   typedef itk::MinimumMaximumImageCalculator<InternalImageType> MinMaxFilterType;

   typedef itk::WarpSmoothnessCalculator<DeformationFieldType>
      SmoothnessCalculatorType;

   typedef itk::VectorCentralDifferenceImageFunction<DeformationFieldType>
      WarpGradientCalculatorType;

   typedef typename WarpGradientCalculatorType::OutputType WarpGradientType;
   
   itkNewMacro( Self );

private:
   std::ofstream m_Fid;
   bool m_headerwritten;
   typename JacobianFilterType::Pointer m_JacobianFilter;
   typename MinMaxFilterType::Pointer m_Minmaxfilter;
   typename SmoothnessCalculatorType::Pointer m_SmothnessCalculator;
   typename DeformationFieldType::ConstPointer m_TrueField;
   typename WarpGradientCalculatorType::Pointer m_TrueWarpGradientCalculator;
   typename WarpGradientCalculatorType::Pointer m_CompWarpGradientCalculator;

public:
   void SetTrueField(const DeformationFieldType * truefield)
   {
      m_TrueField = truefield;

      m_TrueWarpGradientCalculator = WarpGradientCalculatorType::New();
      m_TrueWarpGradientCalculator->SetInputImage( m_TrueField );

      m_CompWarpGradientCalculator =  WarpGradientCalculatorType::New();
   }
   
   void Execute(itk::Object *caller, const itk::EventObject & event)
   {
      Execute( (const itk::Object *)caller, event);
   }

   void Execute(const itk::Object * object, const itk::EventObject & event)
   {
      if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
         return;
      }

      typename DeformationFieldType::ConstPointer deffield = 0;
      unsigned int iter = -1;
      double metricbefore = -1.0;
      
      if ( const DiffeomorphicDemonsRegistrationFilterType * filter1 = 
           dynamic_cast< const DiffeomorphicDemonsRegistrationFilterType * >( object ) )
      {
         iter = filter1->GetElapsedIterations() - 1;
         metricbefore = filter1->GetMetric();
         deffield = const_cast<DiffeomorphicDemonsRegistrationFilterType *>
            (filter1)->GetDeformationField();
      }
      // else if ( const FastSymmetricForcesDemonsRegistrationFilterType * filter2 = 
      //      dynamic_cast< const FastSymmetricForcesDemonsRegistrationFilterType * >( object ) )
      // {
      //    iter = filter2->GetElapsedIterations() - 1;
      //    metricbefore = filter2->GetMetric();
      //    deffield = const_cast<FastSymmetricForcesDemonsRegistrationFilterType *>
      //       (filter2)->GetDeformationField();
      // }
      else if ( const MultiResRegistrationFilterType * multiresfilter = 
           dynamic_cast< const MultiResRegistrationFilterType * >( object ) )
      {
         std::cout<<"Finished Multi-resolution iteration :"<<multiresfilter->GetCurrentLevel()-1<<std::endl;
         std::cout<<"=============================="<<std::endl<<std::endl;
      }
      else
      {
         return;
      }

      if (deffield)
      {
         std::cout<<iter<<": MSE "<<metricbefore<<" - ";

         double fieldDist = -1.0;
         double fieldGradDist = -1.0;
         double tmp;
         if (m_TrueField)
         {
            typedef itk::ImageRegionConstIteratorWithIndex<DeformationFieldType> FieldIteratorType;
            FieldIteratorType currIter( deffield, deffield->GetLargestPossibleRegion() );
            FieldIteratorType trueIter( m_TrueField, deffield->GetLargestPossibleRegion() );

            m_CompWarpGradientCalculator->SetInputImage( deffield );

            fieldDist = 0.0;
            fieldGradDist = 0.0;
            for ( currIter.GoToBegin(), trueIter.GoToBegin();
                  !currIter.IsAtEnd(); ++currIter, ++trueIter )
            {
               fieldDist += (currIter.Value() - trueIter.Value()).GetSquaredNorm();

               // No need to add Id matrix here as we do a substraction
               tmp = ( ( m_CompWarpGradientCalculator->EvaluateAtIndex(currIter.GetIndex()) -m_TrueWarpGradientCalculator->EvaluateAtIndex(trueIter.GetIndex())
                     ).GetVnlMatrix() ).frobenius_norm();
               fieldGradDist += tmp*tmp;
            }
            fieldDist = sqrt( fieldDist/ (double)( deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );
            fieldGradDist = sqrt( fieldGradDist/ (double)( deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );
            
            std::cout<<"d(.,true) "<<fieldDist<<" - ";
            std::cout<<"d(.,Jac(true)) "<<fieldGradDist<<" - ";
         }
         
         m_SmothnessCalculator->SetImage( deffield );
         m_SmothnessCalculator->Compute();
         const double harmonicEnergy = m_SmothnessCalculator->GetSmoothness();
         std::cout<<"harmo. "<<harmonicEnergy<<" - ";

         
         m_JacobianFilter->SetInput( deffield );
         m_JacobianFilter->UpdateLargestPossibleRegion();

        
         const unsigned int numPix = m_JacobianFilter-> GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
         
         TPixel* pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
         TPixel* pix_end = pix_start + numPix;

         TPixel* jac_ptr;

         // Get percentage of det(Jac) below 0
         unsigned int jacBelowZero(0u);
         for (jac_ptr=pix_start; jac_ptr!=pix_end; ++jac_ptr)
         {
            if ( *jac_ptr<=0.0 ) ++jacBelowZero;
         }
         const double jacBelowZeroPrc = static_cast<double>(jacBelowZero) / static_cast<double>(numPix);
         

         // Get min an max jac
         /*
         std::pair<TPixel*, TPixel*> minmax_res =
            boost::minmax_element(pix_start, pix_end);
         */

         //const double minJac = *(minmax_res.first);
         //const double maxJac = *(minmax_res.second);

         const double minJac = *(std::min_element (pix_start, pix_end));
         const double maxJac = *(std::max_element (pix_start, pix_end));

         // Get some quantiles
         // We don't need the jacobian image
         // we can modify/sort it in place
         jac_ptr = pix_start + static_cast<unsigned int>(0.002*numPix);
         std::nth_element(pix_start, jac_ptr, pix_end);
         const double Q002 = *jac_ptr;

         jac_ptr = pix_start + static_cast<unsigned int>(0.01*numPix);
         std::nth_element(pix_start, jac_ptr, pix_end);
         const double Q01 = *jac_ptr;

         jac_ptr = pix_start + static_cast<unsigned int>(0.99*numPix);
         std::nth_element(pix_start, jac_ptr, pix_end);
         const double Q99 = *jac_ptr;

         jac_ptr = pix_start + static_cast<unsigned int>(0.998*numPix);
         std::nth_element(pix_start, jac_ptr, pix_end);
         const double Q998 = *jac_ptr;
         

         std::cout<<"max|Jac| "<<maxJac<<" - "
                  <<"min|Jac| "<<minJac<<" - "
                  <<"ratio(|Jac|<=0) "<<jacBelowZeroPrc<<std::endl;
         
         

         if (this->m_Fid.is_open())
         {
            if (!m_headerwritten)
            {
               this->m_Fid<<"Iteration"
                          <<", MSE before"
                          <<", Harmonic energy"
                          <<", min|Jac|"
                          <<", 0.2% |Jac|"
                          <<", 01% |Jac|"
                          <<", 99% |Jac|"
                          <<", 99.8% |Jac|"
                          <<", max|Jac|"
                          <<", ratio(|Jac|<=0)";
               
               if (m_TrueField)
               {
                  this->m_Fid<<", dist(warp,true warp)"
                             <<", dist(Jac,true Jac)";
               }
               
               this->m_Fid<<std::endl;
               
               m_headerwritten = true;
            }
            
            this->m_Fid<<iter
                       <<", "<<metricbefore
                       <<", "<<harmonicEnergy
                       <<", "<<minJac
                       <<", "<<Q002
                       <<", "<<Q01
                       <<", "<<Q99
                       <<", "<<Q998
                       <<", "<<maxJac
                       <<", "<<jacBelowZeroPrc;

            if (m_TrueField)
            {
               this->m_Fid<<", "<<fieldDist
                          <<", "<<fieldGradDist;
            }
            
            this->m_Fid<<std::endl;
         }
      }
   }
   
protected:   
     
   // Kilian : Do not currently write results back to file  
   //  m_Fid( "metricvalues.csv" ),

   CommandIterationUpdate() :
      m_headerwritten(false)
   {
      m_JacobianFilter = JacobianFilterType::New();
      m_JacobianFilter->SetUseImageSpacing( true );
      m_JacobianFilter->ReleaseDataFlagOn();
      
      m_Minmaxfilter = MinMaxFilterType::New();

      m_SmothnessCalculator = SmoothnessCalculatorType::New();

      m_TrueField = 0;
      m_TrueWarpGradientCalculator = 0;
      m_CompWarpGradientCalculator = 0;
   };

   ~CommandIterationUpdate()
   {
      this->m_Fid.close();
   }
};


void MakeZeroVolume(ImageType::Pointer &zero_image, ImageType::Pointer image)
{
  ImageType::RegionType            region;
  ImageType::IndexType             start;
  region.SetSize( image->GetLargestPossibleRegion().GetSize() );
  start.Fill(0);
  region.SetIndex( start );
  zero_image->SetDirection( image->GetDirection() );
  zero_image->SetOrigin( image->GetOrigin() );
  zero_image->SetSpacing( image->GetSpacing());
  zero_image->SetRegions( region );
  zero_image->Allocate();
  zero_image->FillBuffer( 0.0 );
}


void GetImage(std::string filename, ImageType::Pointer& image)
{
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(filename);
  try
  {
    imageReader->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Could not load volume from disk" << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }
  image = imageReader->GetOutput();
  image->DisconnectPipeline();
}


std::string TemplateName(int i)
{
  std::stringstream result;
  result << "template" << i << ".nii.gz";
  return result.str();
}


void UpdateWriter(WriterType::Pointer &writer, std::string error_msg)
{
  writer->SetUseCompression( true );
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << error_msg << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }
}


void ComputeMean( std::vector<std::string> filenames, std::string filename)
{
  AdderType::Pointer          adder = AdderType::New();
  ImageType::Pointer          template_vol =  ImageType::New();
  ImageType::Pointer          image;

  std::cout << "Computing mean:" << std::endl;
  for (unsigned int i=0; i < filenames.size(); i++)
  {
    std::cout << filenames[i] << std::endl;

    GetImage(filenames[i], image);
    std::cout << image->GetOrigin() << std::endl;

    if (i==0)
    {
      MakeZeroVolume(template_vol, image);
    }
    adder->SetInput1( template_vol );
    adder->SetInput2( image );
    adder->Update();
    template_vol = adder->GetOutput();
  }

  MultiplyByConstantType::Pointer  multiplier = MultiplyByConstantType::New();
  multiplier->SetConstant( 1.0/filenames.size() );
  multiplier->SetInput( template_vol );
  WriterType::Pointer         writer =  WriterType::New();
  writer->SetFileName( filename );
  writer->SetUseCompression( true );
  writer->SetInput( multiplier->GetOutput() );
  UpdateWriter(writer, "Could not write template to disk");
}


std::string FieldName(std::string imageName, int i)
{
  std::stringstream result;
  result << imageName << "_" << i << "_deformation.nii.gz";
  return result.str();
}


std::string WarpedImageName(std::string imageName, int i)
{
  std::stringstream result;
  result << imageName << "_" << i << "_warped.nii.gz";
  return result.str();
}

void ConfigureWarper(WarperType::Pointer& warper, ImageType::Pointer image)
{
  warper->SetInput( image );
  warper->SetOutputSpacing( image->GetSpacing() );
  warper->SetOutputOrigin( image->GetOrigin() );
  warper->SetOutputDirection( image->GetDirection() );
}


void ComputeTemplateFromWarps( std::vector<std::string> filenames, int iter )
{
  typedef itk::DisplacementFieldJacobianDeterminantFilter <DeformationFieldType, float>    JacobianDeterminantFilterType;
  JacobianDeterminantFilterType::Pointer jacdetfilter = JacobianDeterminantFilterType::New();
  jacdetfilter->SetUseImageSpacing( false );
  DeformationReaderType::Pointer   fieldReader = DeformationReaderType::New();
  WriterType::Pointer              writer =  WriterType::New();
  writer->SetUseCompression( true );
  AdderType::Pointer               adder = AdderType::New();
  AdderType::Pointer               jacDetAdder = AdderType::New();
  WarperType::Pointer              warper = WarperType::New();
  ImageType::Pointer               image;
  ImageType::Pointer               imgSum =  ImageType::New();
  ImageType::Pointer               jacDetSum = ImageType::New();

  for (unsigned int i=0; i < filenames.size(); i++)
  {
    std::string field_name = FieldName(filenames[i], iter-1);
    std::string warped_image_name = WarpedImageName(filenames[i], iter-1);
    GetImage(filenames[i], image);
    if (i==0)
    {
      MakeZeroVolume(imgSum, image);
      MakeZeroVolume(jacDetSum, image);
    }
    fieldReader->SetFileName( field_name );
    warper->SetDeformationField( fieldReader->GetOutput() );
    ConfigureWarper(warper, image);
    
    writer->SetFileName( warped_image_name );
    writer->SetInput( warper->GetOutput() );
    UpdateWriter(writer, "Could not write warped image to disk");

    //jacdetfilter->UpdateLargestPossibleRegion();
    jacdetfilter->SetInput( fieldReader->GetOutput() );
    jacDetAdder->SetInput1( jacDetSum );
    jacDetAdder->SetInput2( jacdetfilter->GetOutput() );
    jacDetAdder->Update();
    jacDetSum = jacDetAdder->GetOutput();

    adder->SetInput1( imgSum );
    adder->SetInput2( warper->GetOutput() );
    adder->Update();
    imgSum = adder->GetOutput();
  }
  DivideByImageType::Pointer divider = DivideByImageType::New();
  divider->SetInput1( imgSum );
  divider->SetInput2( jacDetSum );
  writer->SetFileName( TemplateName(iter) );
  writer->SetInput( divider->GetOutput() );
  UpdateWriter(writer, "Could not save the template.");
}


void ComputeTemplate( std::vector<std::string> filenames, int iter )
{
  if (iter == 0)
  {
    ComputeMean( filenames, TemplateName(0) );
    return;
  }
  ComputeTemplateFromWarps( filenames, iter);
}


void NormalizeWarps( std::vector<std::string> filenames, int iter )
{
  DeformationWriterType::Pointer    writer =  DeformationWriterType::New();
  DeformationReaderType::Pointer    defFieldReader = DeformationReaderType::New();
  DeformationAdderType::Pointer     adder = DeformationAdderType::New();
  DeformationFieldType::Pointer     avg_def = DeformationFieldType::New();
  DeformationFieldType::Pointer     image;
  DeformationFieldType::RegionType  region;
  DeformationFieldType::IndexType   start;

  for (unsigned int i=0; i < filenames.size(); i++)
  {
    /* Read in the current image's warp */
    std::stringstream deformation_name;
    deformation_name << filenames[i] << "_" << iter << "_deformation.nii.gz";
    defFieldReader->SetFileName( deformation_name.str() );
    try
    {
      defFieldReader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the deformation fields." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }
    image = defFieldReader->GetOutput();

    if (i==0)
    {
      region.SetSize( image->GetLargestPossibleRegion().GetSize() );
      start.Fill(0);
      region.SetIndex( start );
      avg_def->SetDirection( image->GetDirection() );
      avg_def->SetOrigin( image->GetOrigin() );
      avg_def->SetSpacing( image->GetSpacing() );
      avg_def->SetRegions( region );
      avg_def->Allocate();
      avg_def->FillBuffer( 0.0 );
    }

    adder->SetInput1( avg_def );
    adder->SetInput2( defFieldReader->GetOutput() );
    adder->Update();
    avg_def = adder->GetOutput();
    avg_def->DisconnectPipeline();
  }

  /* Divide the total by the number of volumes, and make it negative */
  DeformationMultiplyByConstantType::Pointer  multiplier = DeformationMultiplyByConstantType::New();
  std::cout << "1 / number of files = " << 1.0/filenames.size() << std::endl;
  multiplier->SetConstant( -1 * 1.0/filenames.size() );
  multiplier->SetInput( avg_def );
  multiplier->Update();
  avg_def = multiplier->GetOutput();


  /* DEBUG Write the avagerage volume to file */
  std::stringstream deformation_name;
  deformation_name << "neg_warp_avg_" << iter << ".nii.gz";
  writer->SetFileName( deformation_name.str() );
  writer->SetUseCompression( true );
  writer->SetInput( adder->GetOutput() );
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }


  for (unsigned int i=0; i < filenames.size(); i++)
  {
    /* Read in the current image's warp */
    std::stringstream deformation_name;
    deformation_name << filenames[i] << "_" << iter << "_deformation.nii.gz";
    defFieldReader->SetFileName( deformation_name.str() );
    try
    {
      defFieldReader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the deformation fields." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }

    adder->SetInput1( avg_def );
    adder->SetInput2( defFieldReader->GetOutput() );
    adder->Update();

    /* Write the new warp to file */
    writer->SetFileName( deformation_name.str() );
    writer->SetUseCompression( true );
    writer->SetInput( adder->GetOutput() );
    try
    {
      writer->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }

  }

}


void WriteDeformationField(DeformationFieldType::Pointer image, std::string filename)
{
  DeformationWriterType::Pointer  writer =  DeformationWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput( image  );
  writer->SetUseCompression( true );
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }
}


std::string DeformationName(std::string filename, int i)
{
  std::stringstream result;
  result << filename << "_" << i << "_" "deformation.nii.gz";
  return result.str();
}


void SetFilterSmoothing(ActualRegistrationFilterType::Pointer& filter, float sigmaDiff, float sigmaUp)
{
  if ( sigmaDiff > 0.1 ) // if ( args.sigmaDef > 0.1 )
  {
    filter->SmoothDeformationFieldOn();
    filter->SetStandardDeviations( sigmaDiff );  //TODO Update with a range
  }
  else
  {
    filter->SmoothDeformationFieldOff();
  }

  if ( sigmaUp > 0.1 )
  {
    filter->SmoothUpdateFieldOn();
    filter->SetUpdateFieldStandardDeviations( sigmaUp );
  }
  else
  {
    filter->SmoothUpdateFieldOff();
  }
}


void ConfigureFilter(ActualRegistrationFilterType::Pointer& filter, float sigmaDiff, float sigmaUp)
{
  filter->SetMaximumUpdateStepLength( 2.0 );
  filter->SetUseGradientType( DemonsRegistrationFunctionType::Symmetric ); //Symmetric is used in ESMInvConDemonsRegistrationFunction 
  SetFilterSmoothing(filter, sigmaDiff, sigmaUp);    
}


void ConfigureMultiResFilter(MultiResRegistrationFilterType::Pointer& multires, arguments args)
{
  multires->SetNumberOfLevels( args.numLevels );
  std::vector<unsigned int> numIterationsVector= std::vector<unsigned int>(args.numLevels, args.numIterations);
  multires->SetNumberOfIterations( &numIterationsVector[0]); // multires->SetNumberOfIterations( &args.numIterations[0] );
}


std::vector<float> Range(float initial, float final, int num)
{  
  float delta = (initial - final)/(num - 1);
  std::vector<float> result(num, 0.0);
  if (delta < 0.0)
  {
    delta = 0.0;
  }
  for (unsigned int i = 0; i < num; ++i)
  {
    result[i] = initial - i * delta;
  }
  return result;
}


void DoGroupWiseRegistration( arguments args, std::string p_arrVolumeNames[], std::string resultsDirectory, std::string outputVolume, bool useJacFlag )
{
  WriterType::Pointer                    writer =  WriterType::New();
  ImageReaderType::Pointer               imageReader = ImageReaderType::New();
  ActualRegistrationFilterType::Pointer  filter;
  std::stringstream deformation_name;
  ImageType::Pointer template_vol = 0; 
  DeformationFieldType::Pointer field = 0;
  std::vector<float> sigma_diff = Range(args.initialSigmaDiff, args.finalSigmaDiff, args.numOuterIterations);

  //{//for mem allocations


  ComputeTemplate(args.volumeFileNames, 0);

  for (unsigned int j = 0; j < args.numOuterIterations; ++j)
  {
    std::cout << "=========================================================================================" << std::endl;
    std::cout << " Starting Outer Iteration  " << j << std::endl;
    std::cout << "=========================================================================================" << std::endl;

    GetImage(TemplateName(j), template_vol);

    for (unsigned int i = 0; i < args.volumeFileNames.size(); ++i)
    {
      std::cout << "=========================================================================================" << std::endl;
      std::cout << " Registering " << args.volumeFileNames[i] << " to template" << j << std::endl;
      std::cout << "=========================================================================================" << std::endl;

      imageReader->SetFileName( args.volumeFileNames[i] );
      
      filter = ActualRegistrationFilterType::New();
      ConfigureFilter(filter, sigma_diff[j], args.sigmaUp);
      std::cout << "sigma_diff[" << j << "] = " << sigma_diff[j] << std::endl;

      MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();
      multires->SetRegistrationFilter( filter );
      ConfigureMultiResFilter(multires, args);
      multires->SetFixedImage( template_vol );
      multires->SetMovingImage( imageReader->GetOutput() );

      if ( args.verbosity > 0 )
      {      
        CommandIterationUpdate<PixelType, Dimension>::Pointer observer = CommandIterationUpdate<PixelType, Dimension>::New();
        CommandIterationUpdate<PixelType, Dimension>::Pointer multiresobserver = CommandIterationUpdate<PixelType, Dimension>::New();
        filter->AddObserver( itk::IterationEvent(), observer );
        multires->AddObserver( itk::IterationEvent(), multiresobserver );
      }

      try
      {
        multires->UpdateLargestPossibleRegion();
      }
      catch( itk::ExceptionObject& err )
      {
        std::cout << "Registration failed unexpectedly" << std::endl;
        std::cout << err << std::endl;
        exit( EXIT_FAILURE );
      }

      WriteDeformationField( multires->GetOutput(), DeformationName(args.volumeFileNames[i], j) );

      /*DEBUG*/
      std::stringstream deformation_name;
      deformation_name << args.volumeFileNames[i] << "_" << j << "_" "deformation-orig.nii.gz";
      WriteDeformationField(multires->GetOutput(), deformation_name.str());
      /*DEBUG*/
    } 

    NormalizeWarps(args.volumeFileNames, j);
    ComputeTemplate(args.volumeFileNames, j+1);
  }
  


  //}//end for mem allocations

}

}


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  for (unsigned int i = 0; i < volumeFileNames.size(); ++i)
    {
    if (!itksys::SystemTools::FileExists(volumeFileNames[i].c_str()))
      {
      std::cerr << "Error: volume file " << i << " does not exist." << std::endl;
      std::cerr << volumeFileNames[i] << std::endl;
      return EXIT_FAILURE;
      }
    }

  if (!itksys::SystemTools::FileExists(resultsDirectory.c_str()))
    {
    std::cerr << "Error: results directory does not exist." << std::endl;
    std::cerr << resultsDirectory << std::endl;
    return EXIT_FAILURE;
    }

  std::string volumeNames[volumeFileNames.size()];
  for (unsigned int i = 0; i < volumeFileNames.size(); ++i)
    {
    volumeNames[i] = volumeFileNames[i].c_str();
    }

 // set up the args structure from the parsed input
   arguments args;
   // args.fixedImageFile = demonsTargetVolume;
   // args.movingImageFile = demonsMovingVolume;
   // args.outputImageFile = demonsResampledVolume;
   // args.outputFieldFile = demonsDeformationVolume;
   // if (iteration.size() > 0)
   //   {
   //   args.numIterations.clear();
   //   }
   // for (unsigned int i = 0; i < iteration.size(); i++)
   //   {
   //   args.numIterations.push_back(iteration[i]);
   //   }
   // args.sigmaDef = smoothing;
   // args.sigmaUp = smoothingUp;
   // args.maxStepLength = maxStepLength;
   // args.useVanillaDem = turnOffDiffeomorph;
   // args.gradientType = gradientType;
   // args.useHistogramMatching = normalization;
   // args.verbosity = verbose;
   args.outputImageFile = outputVolume;
   args.numLevels = numMultiresLevels;
   args.numIterations = numInnerIterations;
   args.volumeFileNames = volumeFileNames;
   args.resultsDirectory = resultsDirectory;
   args.numOuterIterations = numOuterIterations;
   args.initialSigmaDiff = initialSigmaDiff;
   args.finalSigmaDiff = finalSigmaDiff;
   args.regWeight = regWeight;
   args.useJac = useJacFlag;

  DoGroupWiseRegistration(args, volumeNames, resultsDirectory.c_str(), outputVolume.c_str(), useJacFlag);
  //GroupWiseRegistration<3> gpreg;

  return EXIT_SUCCESS;
}


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


template <class TImage, class TDeformationFieldType>
void ComputeTemplate( arguments args, int iter, bool useWarps )
{
  typedef typename itk::ImageFileReader< TImage > ImageReaderType;
  typedef typename itk::ImageFileReader< TDeformationFieldType > DefFieldReaderType;
  typedef itk::AddImageFilter<TImage, TImage, TImage> AdderType;
  typedef itk::ImageFileWriter< TImage >    WriterType;
  typename ImageReaderType::Pointer         imageReader = ImageReaderType::New();
  typename DefFieldReaderType::Pointer      defFieldReader = DefFieldReaderType::New();
  typename WriterType::Pointer              writer =  WriterType::New();
  typename AdderType::Pointer               adder = AdderType::New();
  typename AdderType::Pointer               adder2 = AdderType::New();
  typename TImage::Pointer                  template_vol =  TImage::New();
  typedef itk::WarpImageFilter <TImage, TImage, TDeformationFieldType>      WarperType;
  typename WarperType::Pointer              warper = WarperType::New();
  typedef itk::DisplacementFieldJacobianDeterminantFilter <TDeformationFieldType, float>    JacobianDeterminantFilterType;
  typename JacobianDeterminantFilterType::Pointer jacdetfilter = JacobianDeterminantFilterType::New();
  typename TImage::Pointer 			            jacobianimage = TImage::New();
  typename TImage::Pointer                  jacobian = 0;
  typename TImage::Pointer                  image = 0;
  typename TImage::Pointer                  weight = TImage::New();
  typename TImage::RegionType               region;
  typename TImage::IndexType                start;


  /* For each image in the argument list */
  for (unsigned int i=0; i < args.volumeFileNames.size(); i++)
  {
    std::cout << args.volumeFileNames[i] << std::endl;

    /* Read in the image */
    imageReader->SetFileName( args.volumeFileNames[i] );
    try
    {
      imageReader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the input images." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }
    image = imageReader->GetOutput();


    /* Initialize template_vol as a zero volume */
    if (i==0)
    {
      region.SetSize( image->GetLargestPossibleRegion().GetSize() );
      start.Fill(0);
      region.SetIndex( start );
      template_vol->SetDirection( image->GetDirection() );
      template_vol->SetOrigin( image->GetOrigin() );
      template_vol->SetSpacing( image->GetSpacing());
      template_vol->SetRegions( region );
      template_vol->Allocate();
      template_vol->FillBuffer( 0.0 );
    }

    if (useWarps)
    {
      /* Initialize weight as a zero volume */
      if (i==0)
      {
        region.SetSize( image->GetLargestPossibleRegion().GetSize() );
        start.Fill(0);
        region.SetIndex( start );
        weight->SetDirection( image->GetDirection() );
        weight->SetOrigin( image->GetOrigin() );
        weight->SetSpacing( image->GetSpacing() );
        weight->SetRegions( region );
        weight->Allocate();
        weight->FillBuffer( 0.0 );
      }

      /* Read in the current image's warp */
      std::stringstream deformation_name;
      deformation_name << args.volumeFileNames[i] << "_" << iter-1 << "_deformation.nii.gz";
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

      /* Warp the image */
      warper->SetInput( image );
      warper->SetOutputSpacing( image->GetSpacing() );
      warper->SetOutputOrigin( image->GetOrigin() );
      warper->SetOutputDirection( image->GetDirection() );
      warper->SetDeformationField( defFieldReader->GetOutput() );
      image = warper->GetOutput();

      /* Write the warped image to disk */
      std::stringstream warped_image_name;
      warped_image_name << args.volumeFileNames[i] << "_" << iter-1 << "_warped.nii.gz";
      writer->SetFileName( warped_image_name.str() );
      writer->SetUseCompression( true );
      writer->SetInput( image );
      try
      {
        writer->Update();
      }
      catch( itk::ExceptionObject& err )
      {
        std::cout << "Could not write warped image to disk." << std::endl;
        std::cout << err << std::endl;
        exit( EXIT_FAILURE );
      }

      /* Compute the jacobian of the warp and add to the running total (weight) */
      jacdetfilter->SetInput( defFieldReader->GetOutput() );
      jacdetfilter->SetUseImageSpacing( false );
      jacdetfilter->UpdateLargestPossibleRegion();
      adder2->SetInput1( weight );
      adder2->SetInput2( jacdetfilter->GetOutput() );
      adder2->Update();
      weight = adder2->GetOutput();
    }

    /* Add the image to the running total (template_vol) */
    adder->SetInput1( template_vol );
    adder->SetInput2( image );
    adder->Update();
    template_vol = adder->GetOutput();
  }

  /* Divide template_vol by the weight */
  if (useWarps)
  {
    typedef itk::DivideImageFilter<TImage, TImage, TImage>   DivideByImageType;
    typename DivideByImageType::Pointer  divider = DivideByImageType::New();
    divider->SetInput1( template_vol );
    divider->SetInput2( weight );
    divider->Update();
    template_vol = divider->GetOutput();
  }
  else
  {
    typedef itk::MultiplyByConstantImageFilter< TImage, float, TImage >   MultiplyByConstantType;
    typename MultiplyByConstantType::Pointer  multiplier = MultiplyByConstantType::New();
    multiplier->SetConstant( 1.0/args.volumeFileNames.size() );
    multiplier->SetInput( template_vol );
    multiplier->Update();
    template_vol = multiplier->GetOutput();
  }

  /* Write the template_vol to disk (e.g. template0.nii.gz) */
  std::stringstream template_name;
  template_name << "template" << iter << ".nii.gz";
  writer->SetFileName( template_name.str() );
  writer->SetUseCompression( true );
  writer->SetInput( template_vol );
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


template <unsigned int Dimension>
void DoGroupWiseRegistration( arguments args, std::string p_arrVolumeNames[], std::string resultsDirectory, std::string outputVolume, bool useJacFlag )
{
  typedef float PixelType;
  typedef itk::OrientedImage< PixelType, Dimension >  OrientedImageType;// read and transform oriented images, operate on non-oriented
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef float VectorComponentType;
  typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
  typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;
  typedef itk::ESMInvConDemonsRegistrationFunction <ImageType, ImageType,DeformationFieldType>  DemonsRegistrationFunctionType;
  typedef typename itk::ImageFileReader< ImageType > ImageReaderType;
  typedef typename itk::ImageFileReader< DeformationFieldType > DeformationReaderType;
  typedef itk::ImageFileWriter< DeformationFieldType >  DeformationWriterType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typename WriterType::Pointer      writer =  WriterType::New();
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  std::stringstream template_name;
  std::stringstream deformation_name;

  //{//for mem allocations


 

  /* Compute the initial template */
  ComputeTemplate<ImageType, DeformationFieldType>(args, 0, false);


  /* For number of outer iterations */
  for (unsigned int j = 0; j < args.numOuterIterations; ++j)
  {

    

    /*  Load the template from disk */
    template_name.str("");
    template_name << "template" << j << ".nii.gz";
    typename ImageReaderType::Pointer imageReader2 = ImageReaderType::New();
    imageReader2->SetFileName(template_name.str());
    try
    {
      imageReader2->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the input images." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }
    typename ImageType::Pointer                  template_vol = 0;
    template_vol = imageReader2->GetOutput();
    template_vol->DisconnectPipeline();

    /* For each image in the argument list */
    for (unsigned int i = 0; i < args.volumeFileNames.size(); ++i)
    {
      /* Load the image */
      imageReader->SetFileName( args.volumeFileNames[i] );
      try
      {
        imageReader->Update();
      }
      catch( itk::ExceptionObject& err )
      {
        std::cout << "Could not read one of the input images." << std::endl;
        std::cout << err << std::endl;
        exit( EXIT_FAILURE );
      }


      /* Create and configure the diffeomorphic demons filter (filter) */
      typedef typename itk::DiffeomorphicDemonsRegistrationFilter < ImageType, ImageType, DeformationFieldType>   ActualRegistrationFilterType;
      typename         ActualRegistrationFilterType::Pointer          filter = ActualRegistrationFilterType::New();
      typename         DemonsRegistrationFunctionType::GradientType   gtype =  DemonsRegistrationFunctionType::Symmetric;
      filter->SetMaximumUpdateStepLength( 100.0 ); //TODO: check this //filter->SetMaximumUpdateStepLength( args.maxStepLength );
      // filter->SetMaximumUpdateStepLength( 0.1 ); //TODO: check this //filter->SetMaximumUpdateStepLength( args.maxStepLength );
      //filter->SetRegWeight( 10.0 );
      filter->SetUseGradientType( gtype ); //Symmetric is used in ESMInvConDemonsRegistrationFunction 
      // if ( args.sigmaDef > 0.1 )
      if ( args.initialSigmaDiff > 0.1 )
      {
        filter->SmoothDeformationFieldOn();
        filter->SetStandardDeviations( args.initialSigmaDiff );  //TODO Update with a range
      }
      else
        filter->SmoothDeformationFieldOff();

      if ( args.sigmaUp > 0.1 )
      {
        filter->SmoothUpdateFieldOn();
        filter->SetUpdateFieldStandardDeviations( args.sigmaUp );
      }
      else
        filter->SmoothUpdateFieldOff();

      // filter->SetInitialDeformationField( deformationReader->GetOutput() );

      /* Create and configure the multi-resolution filter and compute the deformation field */ 
      typedef typename itk::MultiResolutionPDEDeformableRegistration2< ImageType, ImageType, DeformationFieldType, PixelType >   MultiResRegistrationFilterType;
      typename MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();
      multires->SetRegistrationFilter( filter );
      multires->SetNumberOfLevels( args.numLevels );
      std::vector<unsigned int> numIterationsVector= std::vector<unsigned int>(args.numLevels, args.numIterations);
      multires->SetNumberOfIterations( &numIterationsVector[0]); // multires->SetNumberOfIterations( &args.numIterations[0] );

      if ( args.verbosity > 0 )
      {      
        typename CommandIterationUpdate<PixelType, Dimension>::Pointer observer = CommandIterationUpdate<PixelType, Dimension>::New();
        filter->AddObserver( itk::IterationEvent(), observer );
        typename CommandIterationUpdate<PixelType, Dimension>::Pointer multiresobserver = CommandIterationUpdate<PixelType, Dimension>::New();
        multires->AddObserver( itk::IterationEvent(), multiresobserver );
      }

      /* Compute the warp from the image to the template */
      multires->SetFixedImage( template_vol );
      multires->SetMovingImage( imageReader->GetOutput() );
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
      typename DeformationFieldType::Pointer defField = 0;
      defField = multires->GetOutput();
      defField->DisconnectPipeline();

      /* Save the warp, e.g. image1_0_deformation.nii.gz */
      deformation_name.str("");
      deformation_name << args.volumeFileNames[i] << "_" << j << "_" "deformation.nii.gz";
      typename DeformationWriterType::Pointer  defWriter =  DeformationWriterType::New();
      defWriter->SetFileName( deformation_name.str() );
      defWriter->SetInput( defField  );
      defWriter->SetUseCompression( true );
      try
      {
        defWriter->Update();
      }
      catch( itk::ExceptionObject& err )
      {
        std::cout << "Unexpected error." << std::endl;
        std::cout << err << std::endl;
        exit( EXIT_FAILURE );
      }

    } // end for each image

    /* Renormalize the warps TODO */
    // if (renormalize_warps_flag)
    //     renormalize_warps(SBJ_CELL, Output_dir, CurTemplateWarpName);
    // end

    /* Compute the new template */
    ComputeTemplate<ImageType, DeformationFieldType>(args, j+1, true);
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

  DoGroupWiseRegistration<3>(args, volumeNames, resultsDirectory.c_str(), outputVolume.c_str(), useJacFlag);

  return EXIT_SUCCESS;
}


  // std::stringstream warped_image_name;
  // warped_image_name << args.volumeFileNames[0] << "_" << "4" << "_warped.nii.gz";
  // typename ImageReaderType::Pointer reader = ImageReaderType::New();
  // reader->SetFileName(warped_image_name.str());
  // try
  // {
  //   reader->Update();
  // }
  // catch( itk::ExceptionObject& err )
  // {
  //   std::cout << "Could not read warped image." << std::endl;
  //   std::cout << err << std::endl;
  //   exit( EXIT_FAILURE );
  // }
  // warped_image_name.str("");
  // warped_image_name << args.volumeFileNames[1] << "_" << "4" << "_warped.nii.gz";
  // typename ImageReaderType::Pointer reader2 = ImageReaderType::New();
  // reader2->SetFileName(warped_image_name.str());
  // try
  // {
  //   reader2->Update();
  // }
  // catch( itk::ExceptionObject& err )
  // {
  //   std::cout << "Could not read warped image." << std::endl;
  //   std::cout << err << std::endl;
  //   exit( EXIT_FAILURE );
  // }

  // typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AdderType;
  // typename AdderType::Pointer               adder = AdderType::New();
  // adder->SetInput1( reader->GetOutput() );
  // adder->SetInput2( reader2->GetOutput() );
  // adder->Update();
  //     writer->SetFileName( "test_template.nii.gz" );
  //     writer->SetInput( adder->GetOutput()  );
  //     writer->SetUseCompression( true );
  //     try
  //     {
  //       writer->Update();
  //     }
  //     catch( itk::ExceptionObject& err )
  //     {
  //       std::cout << "Unexpected error." << std::endl;
  //       std::cout << err << std::endl;
  //       exit( EXIT_FAILURE );
  //     }
  //     exit(0);


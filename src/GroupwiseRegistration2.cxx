#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

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

// Declare the types of the images
typedef float PixelType;
typedef itk::OrientedImage< PixelType, Dimension >  OrientedImageType;// read and transform oriented images, operate on non-oriented
typedef itk::OrientImageFilter<OrientedImageType,OrientedImageType> OrientFilterType;//##
typedef itk::Image< PixelType, Dimension > ImageType;

template <unsigned int Dimension>
void DoGroupWiseRegistration( arguments args, std::string p_arrVolumeNames[], std::string resultsDirectory, std::string outputVolume, bool useJacFlag )
{


    //{//for mem allocations

  // Set up the file readers
    typedef typename itk::ImageFileReader< OrientedImageType > ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader  = ImageReaderType::New();
    typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( p_arrVolumeNames[0]);  //TODO: test array
    movingImageReader->SetFileName( p_arrVolumeNames[1] );
    try
    {
      fixedImageReader->Update();
      movingImageReader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the input images." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }

  // Set up pyramids
  
  // Begin inconvforces
     typedef float                                        VectorComponentType;
     typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
     typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;
     typedef itk::ESMInvConDemonsRegistrationFunction <ImageType,ImageType,DeformationFieldType>  DemonsRegistrationFunctionType;
     typedef itk::DisplacementFieldJacobianDeterminantFilter <DeformationFieldType, PixelType>    JacobianDeterminantFilterType;
     typedef itk::WarpImageFilter <ImageType, ImageType, DeformationFieldType>      WarperType;

     typename DeformationFieldType::Pointer  log_field = DeformationFieldType::New();
     typename DeformationFieldType::Pointer  log_inv_field = DeformationFieldType::New();
     typename DeformationFieldType::Pointer  update = DeformationFieldType::New();
     typename WarperType::Pointer            warper = WarperType::New();
     typename ImageType::Pointer  		       fixedimage = fixedImageReader->GetOutput();
     typename ImageType::Pointer 			       movingimage = movingImageReader->GetOutput();
     typename ImageType::Pointer 			       jacobianimage = ImageType::New();

     //invcondemonsforces<PixelType, ImageType, DeformationFieldType, 3>(fixedimage, movingimage, field, inv_field, 3, 0, 3);
     typename DeformationFieldType::SpacingType 	 spacing;
     spacing.Fill( 1.0 );
     typename DeformationFieldType::PointType   	 origin;
     //origin.Fill( 0.0 );
     //typename DeformationFieldType::SpacingType 	 spacing = movingimage->GetSpacing();
     //typename DeformationFieldType::PointType   	 origin  = movingimage->GetOrigin();
     typename DeformationFieldType::RegionType     region;
     typename DeformationFieldType::SizeType       size;
     typename DeformationFieldType::IndexType      start;

     std::cout<<"GroupWiseRegistration:Fixed image origin: "<<fixedimage->GetOrigin()<<std::endl;
     std::cout<<"GroupWiseRegistration:Fixed image size: "<<fixedimage->GetLargestPossibleRegion().GetSize()[0] << " " <<
         fixedimage->GetLargestPossibleRegion().GetSize()[1] << " " << fixedimage->GetLargestPossibleRegion().GetSize()[2] << std::endl;
     size = movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
     for (unsigned int d=0; d<Dimension; d++)
     {
       start[d] = 0;
       origin[d] = -127;
     }
     //start.Fill(0);
     region.SetSize( size );
     region.SetIndex( start );


     //todo: test this
     VectorPixelType initialvalue;
     initialvalue.Fill( 0.0 );
     //
     log_field->SetOrigin( origin );
     log_field->SetSpacing( spacing );
     log_field->SetRegions( region );
     log_field->Allocate();
     log_field->FillBuffer( 0.0 );
     update->SetOrigin( origin );
     update->SetSpacing( spacing );
     update->SetRegions( region );
     update->Allocate();
     update->FillBuffer( 0.0 );

  // Compute exponentials
     typedef itk::ExponentialDeformationFieldImageFilter <DeformationFieldType, DeformationFieldType>      FieldExponentiatorType;
     typename FieldExponentiatorType::Pointer exponentiatorOfLogField= FieldExponentiatorType::New();
     typename FieldExponentiatorType::Pointer exponentiatorOfNegativeLogField = FieldExponentiatorType::New();

     //exponentiatorOfLogField->SetInput( log_field );
     //exponentiatorOfLogField->UpdateLargestPossibleRegion();
     exponentiatorOfLogField->AutomaticNumberOfIterationsOn();
     exponentiatorOfLogField->SetMaximumNumberOfIterations( 2000u ); // Just set a high value so that automatic number of step // is not thresholded
     exponentiatorOfNegativeLogField->AutomaticNumberOfIterationsOn();
     exponentiatorOfNegativeLogField->SetMaximumNumberOfIterations( 2000u ); // Just set a high value so that automatic number of step // is not thresholded
     typename DeformationFieldType::Pointer field; //= exponentiatorOfLogField->GetOutput();
     typename DeformationFieldType::Pointer inv_field; //= exponentiatorOfNegativeLogField->GetOutput();

     typedef itk::MultiplyByConstantImageFilter< DeformationFieldType, signed int, DeformationFieldType > MultiplyByConstantType;
     typedef typename MultiplyByConstantType::Pointer   MultiplyByConstantPointer;
     MultiplyByConstantPointer m_Multiplier = MultiplyByConstantType::New();
     m_Multiplier->SetConstant( -1 );

  // Compute jacobian determinant
     typename JacobianDeterminantFilterType::Pointer jacdetfilter = JacobianDeterminantFilterType::New();
     typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxFilterType;
     typename MinMaxFilterType::Pointer minmaxfilter = MinMaxFilterType::New();
    // if (useJacFlag > 0)
    // {
    //   jacdetfilter->SetInput( field );
    //   jacdetfilter->SetUseImageSpacing( false );
    //   jacdetfilter->UpdateLargestPossibleRegion();
    //   jacobianimage = jacdetfilter->GetOutput();
    //   minmaxfilter->SetImage( jacdetfilter->GetOutput() );
    //   minmaxfilter->Compute();
    //   std::cout<<"Minimum of the determinant of the Jacobian of the warp: " <<minmaxfilter->GetMinimum()<<std::endl;
    //   std::cout<<"Maximum of the determinant of the Jacobian of the warp: " <<minmaxfilter->GetMaximum()<<std::endl;
    // }

  // Set up demons registration function
     typename DemonsRegistrationFunctionType::Pointer drfp = DemonsRegistrationFunctionType::New();
     typename DemonsRegistrationFunctionType::GradientType gtype = DemonsRegistrationFunctionType::Symmetric;
     drfp->SetUseGradientType( gtype );
     drfp->SetRegWeight( args.regWeight );

  // set up demons filter 
     typedef typename itk::DiffeomorphicDemonsRegistrationFilter < ImageType, ImageType, DeformationFieldType>   ActualRegistrationFilterType;
     //typedef typename ActualRegistrationFilterType::GradientType GradientType;
     typename ActualRegistrationFilterType::Pointer filter = ActualRegistrationFilterType::New();
     //filter->SetMaximumUpdateStepLength( args.maxStepLength );
     filter->SetMaximumUpdateStepLength( 1e1 ); //TODO: check this
     //filter->SetUseGradientType( static_cast<GradientType>(args.gradientType) );
     filter->SetUseGradientType( gtype );
    if ( args.sigmaDef > 0.1 )
       {
         filter->SmoothDeformationFieldOn();
         filter->SetStandardDeviations( args.sigmaDef );
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

       if ( args.verbosity > 0 )
       {      
         typename CommandIterationUpdate<PixelType, Dimension>::Pointer observer = CommandIterationUpdate<PixelType, Dimension>::New();
         filter->AddObserver( itk::IterationEvent(), observer );
       }

 // Set up the multi-resolution filter
     typedef typename itk::MultiResolutionPDEDeformableRegistration2< ImageType, ImageType, DeformationFieldType, PixelType >   MultiResRegistrationFilterType;
     typename MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();
     multires->SetRegistrationFilter( filter );
     multires->SetNumberOfLevels( args.numLevels );

     std::vector<unsigned int> numIterationsVector= std::vector<unsigned int>(args.numLevels, args.numIterations);
     multires->SetNumberOfIterations( &numIterationsVector[0]); // multires->SetNumberOfIterations( &args.numIterations[0] );
     multires->SetFixedImage( fixedimage );
     multires->SetMovingImage( movingimage );
     if ( args.verbosity > 0 )
     {
       typename CommandIterationUpdate<PixelType, Dimension>::Pointer multiresobserver = CommandIterationUpdate<PixelType, Dimension>::New();
       multires->AddObserver( itk::IterationEvent(), multiresobserver );
     }

 //Compute the deformation field
   try
   {
      multires->UpdateLargestPossibleRegion();
   }
   catch( itk::ExceptionObject& err )
   {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
   }
     std::cout << "numiterations: " << args.numIterations << std::endl;


   exit(0);

 // The outputs
   typename DeformationFieldType::Pointer defField = 0;
   defField = multires->GetOutput();
   defField->DisconnectPipeline();

   //}//end for mem allocations

   
   // warp the result
   typedef itk::WarpImageFilter < ImageType, ImageType, DeformationFieldType >  WarperType;
   warper->SetInput( movingimage );
   warper->SetInput( movingImageReader->GetOutput() );
   warper->SetOutputSpacing( spacing );
   warper->SetOutputOrigin( origin );
   // warper->SetOutputSpacing( fixedimage->GetSpacing() );
   // warper->SetOutputOrigin( fixedimage->GetOrigin() );
   warper->SetDeformationField( defField );
   
   // Write warped image out to file
   typedef PixelType OutputPixelType;
   typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
   typedef itk::CastImageFilter < ImageType, OutputImageType > CastFilterType;
   typedef itk::ImageFileWriter< OutputImageType >  WriterType;
   
   typename WriterType::Pointer      writer =  WriterType::New();
   typename CastFilterType::Pointer  caster =  CastFilterType::New();
   writer->SetFileName( args.outputImageFile.c_str() );
   caster->SetInput( warper->GetOutput() );
   writer->SetInput( caster->GetOutput()   );
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



  // Run demons registration 
     int numIter = 1;
     int strike = 0;
     float mse[numIter]; 
     for (int i = 0; i < numIter; i++) 
     {
       std::cout<<"i = "<<i<<std::endl;

       exponentiatorOfLogField->SetInput( log_field );
       exponentiatorOfLogField->UpdateLargestPossibleRegion();
       field = exponentiatorOfLogField->GetOutput();

       // exponentiatorOfNegativeLogField->SetInput( log_field ); //TODO: convert to negative of log_field
       // exponentiatorOfNegativeLogField->UpdateLargestPossibleRegion();
       m_Multiplier->SetInput( log_field );
       //m_Multiplier->GraftOutput( update ); //debugging
       m_Multiplier->Update();
       exponentiatorOfNegativeLogField->SetInput( m_Multiplier->GetOutput() ); //TODO: check this is negative of log_field
       inv_field = exponentiatorOfNegativeLogField->GetOutput();
       jacdetfilter->SetInput( field );
       jacdetfilter->SetUseImageSpacing( false );
       jacdetfilter->UpdateLargestPossibleRegion();
       jacobianimage = jacdetfilter->GetOutput();
       minmaxfilter->SetImage( jacobianimage );
       minmaxfilter->Compute();
       std::cout<<"Minimum of the determinant of the Jacobian of the warp: " <<minmaxfilter->GetMinimum()<<std::endl;
       std::cout<<"Maximum of the determinant of the Jacobian of the warp: " <<minmaxfilter->GetMaximum()<<std::endl;

      /*Debug*/
      // try
      // {
      //   typedef itk::ImageFileWriter< DeformationFieldType>  FieldWriterType;
      //   typename FieldWriterType::Pointer      fieldwriter =  FieldWriterType::New();
      //   fieldwriter->SetUseCompression( true );
      //   fieldwriter->SetFileName( "field.nii.gz" );
      //   fieldwriter->SetInput( field  );
      //   fieldwriter->Update();
      //   fieldwriter->SetFileName( "inv_field.nii.gz" );
      //   fieldwriter->SetInput( inv_field );
      //   fieldwriter->Update();
      // }
      // catch( itk::ExceptionObject& err )
      // {
      //    std::cout << "Unexpected error." << std::endl;
      //    std::cout << err << std::endl;
      //    exit( EXIT_FAILURE );
      // }
      // exit(0);
     /*End Debug*/


       /**Debug**/
       //typedef itk::MinimumMaximumImageCalculator<DeformationFieldType> MinMaxFilterType2;
       //typename MinMaxFilterType2::Pointer minmaxfilter2 = MinMaxFilterType2::New();
       //minmaxfilter2->SetImage( m_Multiplier->GetOutput() );
       //minmaxfilter2->Compute();
       //std::cout<<"Minimum of the inv_field: " <<minmaxfilter2->GetMinimum()<<std::endl;
       //std::cout<<"Maximum of the inv_field: " <<minmaxfilter2->GetMaximum()<<std::endl;

       //minmaxfilter2->SetImage( log_field );
       //minmaxfilter2->Compute();
       //std::cout<<"Minimum of the log_field: " <<minmaxfilter2->GetMinimum()<<std::endl;
       //std::cout<<"Maximum of the log_field: " <<minmaxfilter2->GetMaximum()<<std::endl;

       /**End Debug**/

       drfp->SetDeformationField( field );
       drfp->SetInvDeformationField( inv_field );
       drfp->SetFixedImage( fixedImageReader->GetOutput() );
       drfp->SetMovingImage( movingImageReader->GetOutput() );
       if (useJacFlag > 0)
       {
         drfp->SetUseJacobian(true);
         drfp->SetJacobianDetImage(jacobianimage);
       }
       else
       {
         drfp->SetUseJacobian(false);
       }
       drfp->InitializeIteration();

       const itk::Size<Dimension> radius = drfp->GetRadius();

       // Break the input into a series of regions.  The first region is free
       // of boundary conditions, the rest with boundary conditions.  We operate
       // on the output region because input has been copied to output.
       typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator <DeformationFieldType> FaceCalculatorType;
       typedef typename FaceCalculatorType::FaceListType FaceListType;
       typedef typename DemonsRegistrationFunctionType::NeighborhoodType NeighborhoodIteratorType;
       typedef itk::ImageRegionIterator<DeformationFieldType> UpdateIteratorType;

       FaceCalculatorType faceCalculator;

       FaceListType faceList = faceCalculator(field, region, radius);
       typename FaceListType::iterator fIt = faceList.begin();

       // Ask the function object for a pointer to a data structure it
       // will use to manage any global values it needs.  We'll pass this
       // back to the function object at each calculation and then
       // again so that the function object can use it to determine a
       // time step for this iteration.
       void * globalData = drfp->GetGlobalDataPointer();

       // Process the non-boundary region.
       NeighborhoodIteratorType nD(radius,field, *fIt);
       UpdateIteratorType       nU(update,  *fIt);
       nD.GoToBegin();
       while( !nD.IsAtEnd() )
       {
         nU.Value() = drfp->ComputeUpdate(nD, globalData);
         ++nD;
         ++nU;
       }

       // Process each of the boundary faces.
       NeighborhoodIteratorType bD;
       UpdateIteratorType   bU;
       for (++fIt; fIt != faceList.end(); ++fIt)
       {
         bD = NeighborhoodIteratorType(radius,field, *fIt);
         bU = UpdateIteratorType(update, *fIt);

         bD.GoToBegin();
         bU.GoToBegin();
         while ( !bD.IsAtEnd() )
         {
           bU.Value() = drfp->ComputeUpdate(bD, globalData);
           ++bD;
           ++bU;
         }
       }
       // Ask the finite difference function to compute the time step for
       // this iteration.  We give it the global data pointer to use, then
       // ask it to free the global data memory.
       //timeStep = df->ComputeGlobalTimeStep(globalData);
       drfp->ReleaseGlobalDataPointer(globalData);


       // Warp the image  (In MATLAB: warped_mov_im = warpimage(mov_im, def_x, def_y, def_z);)
       warper->SetInput( movingImageReader->GetOutput() );
       warper->SetOutputSpacing( spacing );
       warper->SetOutputOrigin( origin );
       //warper->SetOutputSpacing( fixedimage->GetSpacing() );
       //warper->SetOutputOrigin( fixedimage->GetOrigin() );
       //warper->SetDeformationField( exponentiatorOfLogField->GetOutput() );
       //warper->SetDeformationField( update );
       warper->SetDeformationField( field );
       if ( std::numeric_limits<PixelType>::has_quiet_NaN )
       {
         warper->SetEdgePaddingValue( std::numeric_limits<PixelType>::quiet_NaN() );
       }
       warper->UpdateLargestPossibleRegion();
       warper->Update();

       // Calculate MSE
       double finalSSD = 0.0;
       typedef itk::ImageRegionConstIterator<ImageType> ImageConstIterator;
       ImageConstIterator iterfix = ImageConstIterator( fixedimage, fixedimage->GetRequestedRegion() );
       ImageConstIterator itermovwarp = ImageConstIterator( warper->GetOutput(), fixedimage->GetRequestedRegion() );
       for (iterfix.GoToBegin(), itermovwarp.GoToBegin(); !iterfix.IsAtEnd(); ++iterfix, ++itermovwarp)
       {
         if (iterfix.Get() == iterfix.Get() && itermovwarp.Get() == itermovwarp.Get())
           finalSSD += vnl_math_sqr( iterfix.Get() - itermovwarp.Get() );
       }
       const double finalMSE = finalSSD / static_cast<double>( fixedimage->GetRequestedRegion().GetNumberOfPixels() );
       mse[i] = finalMSE;
       std::cout<<"****MSE fixed image vs. warped moving image: "<<finalMSE<<std::endl;


       // Smooth update
       //double var[Dimension];
       //for (unsigned int i=0; i<Dimension; ++i)
       //  var[i] = 0.6;
       //typedef itk::SmoothingRecursiveGaussianImageFilter<DeformationFieldType,DeformationFieldType> GaussianFilterType;
       //typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
       //smoother->SetInput( update );
       //smoother->SetSigma( 0.6 ); //TODO change this
       //smoother->Update();


       // Update deformation field
       typedef itk::AddImageFilter< DeformationFieldType, DeformationFieldType, DeformationFieldType > AddFilterType;
       typename AddFilterType::Pointer addFilter = AddFilterType::New();
       addFilter->SetInput1( log_field );
       //addFilter->SetInput2( smoother->GetOutput() );
       addFilter->SetInput2( update );
       try
       {
         addFilter->Update();
       }
       catch( itk::ExceptionObject & err )
       {
         std::cout << "ExceptionObject caught !" << std::endl;
         std::cout << err << std::endl;
         exit( EXIT_FAILURE );
       }
       log_field = addFilter->GetOutput();

       // Test add images (works)
       /*typedef itk::AddImageFilter< ImageType, ImageType, ImageType > AddImageFilterType;
         typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
         addImageFilter->SetInput1( movingImageReader->GetOutput() );
       //addFilter->SetInput2( smoother->GetOutput() );
       addImageFilter->SetInput2( warper->GetOutput() );
       try
       {
       addImageFilter->Update();
       }
       catch( itk::ExceptionObject & err )
       {
       std::cout << "ExceptionObject caught !" << std::endl;
       std::cout << err << std::endl;
       exit( EXIT_FAILURE );
       }*/

       if ( i > 1)
       {
         if (mse[i] > mse[i - 1])
           strike++;
       }

       if (strike > 3)
         break;
     }
     //end for options.numiter
     
  // Write image
     // typedef itk::ImageFileWriter< ImageType >  WriterType;
    // typename WriterType::Pointer      writer =  WriterType::New();
    writer->SetFileName( outputVolume );
    //writer->SetInput( movingImageReader->GetOutput()  );
    writer->SetInput( warper->GetOutput()  );
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

#ifndef _itkDiffeomorphicDemonsRegistrationFilter_txx
#define _itkDiffeomorphicDemonsRegistrationFilter_txx
#include "itkDiffeomorphicDemonsRegistrationFilter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include <itkImageDuplicator.h>

namespace itk {

/*
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::DiffeomorphicDemonsRegistrationFilter()
{
  typename DemonsRegistrationFunctionType::Pointer drfp;
  drfp = DemonsRegistrationFunctionType::New();
  // drfp->SetRegWeight(this->GetRegWeight());
  // drfp->SetRegWeight(0.1); //inverse error
  drfp->SetRegWeight(1.0);
  // std::cout << "Reg Weight is " << this->GetRegWeight() << std::endl;

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>( drfp.GetPointer() ) );

  m_Multiplier = MultiplyByConstantType::New();
  // m_Multiplier->InPlaceOn();  //This was breaking the code in ApplyUpdate(...)

  m_Exponentiator = FieldExponentiatorType::New();
  
  m_Warper = VectorWarperType::New();
  FieldInterpolatorPointer VectorInterpolator =
     FieldInterpolatorType::New();
  m_Warper->SetInterpolator(VectorInterpolator);

  m_Adder = AdderType::New();
  m_Adder->InPlaceOn();

  /* Set the inverse deformation field to null */
  this->SetInvDeformationField(NULL);
}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{

  std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration(): " << std::endl;

  // update variables in the equation object
  DemonsRegistrationFunctionType *f = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());

  if ( !f )
  {
    itkExceptionMacro(<<"FiniteDifferenceFunction not of type DemonsRegistrationFunctionType");
  }

  if (this->GetInvDeformationField() == NULL)
  {
    typedef ImageDuplicator<TDeformationField> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(this->GetDeformationField());
    duplicator->Update();
    m_invDeformationField = duplicator->GetOutput();
    m_invDeformationField->FillBuffer(itk::NumericTraits< typename TDeformationField::PixelType >::Zero);
    //typename TDeformationField::RegionType            region;
    //typename TDeformationField::IndexType             start;
    //region.SetSize( this->GetDeformationField()->GetLargestPossibleRegion().GetSize() );
    //start.Fill(0);
    //region.SetIndex( start );
    //m_invDeformationField->SetDirection( this->GetDeformationField()->GetDirection() );  //this causes a seg fault 
    //m_invDeformationField->SetOrigin( this->GetDeformationField()->GetOrigin() );
    //m_invDeformationField->SetSpacing( this->GetDeformationField()->GetSpacing());

    //m_invDeformationField->SetRegions( region );
    //m_invDeformationField->Allocate();
    //m_invDeformationField->FillBuffer( 0.0 );

    //this->SetInvDeformationField( );
    //this->SetInvDeformationField( this->GetDeformationField() ); //should be zeros
  }

  f->SetDeformationField( this->GetDeformationField() ); 
  f->SetInvDeformationField( this->GetInvDeformationField() );

  /*Debug*/
  //try
  //{
    //typedef itk::ImageFileWriter< DeformationFieldType>  FieldWriterType;
    //typename FieldWriterType::Pointer      fieldwriter =  FieldWriterType::New();
    //fieldwriter->SetUseCompression( true );
    //fieldwriter->SetFileName( "invdeformation.nii.gz" );
    //fieldwriter->SetInput( this->GetInvDeformationField()  );
    //fieldwriter->Update();
    //fieldwriter->SetFileName( "deformation.nii.gz" );
    //fieldwriter->SetInput( this->GetDeformationField()  );
    //fieldwriter->Update();
  //}
  //catch( itk::ExceptionObject& err )
  //{
    //std::cout << "Unexpected error." << std::endl;
    //std::cout << err << std::endl;
    //exit( EXIT_FAILURE );
  //}
  /*End Debug*/


  // call the superclass  implementation ( initializes f )
  Superclass::InitializeIteration();


  std::cout << "Finished InitializeIteration" << std::endl;
}


/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GetMetric() const
{

  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());

  if( !drfp )
  {
    itkExceptionMacro( << 
        "Could not cast difference function to DiffeomorphicDemonsRegistrationFunction" );
  }

  return drfp->GetMetric();
}

/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GetIntensityDifferenceThreshold() const
{
 
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to DemonsRegistrationFunction" );
   }
   
  return drfp->GetIntensityDifferenceThreshold();
}

/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::SetIntensityDifferenceThreshold(double threshold) 
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
   }
   
  drfp->SetIntensityDifferenceThreshold(threshold);
}


/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GetMaximumUpdateStepLength() const
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to DemonsRegistrationFunction" );
   }
  
  return drfp->GetMaximumUpdateStepLength();
}

/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::SetMaximumUpdateStepLength(double threshold) 
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
   }
  
  drfp->SetMaximumUpdateStepLength(threshold);
}


/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
const double &
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GetRMSChange() const
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to DiffeomorphicDemonsRegistrationFunction" );
   }
   
  return drfp->GetRMSChange();
}


/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
typename DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GradientType
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::GetUseGradientType() const
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to DemonsRegistrationFunction" );
   }
  
  return drfp->GetUseGradientType();
}

/*
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::SetUseGradientType(GradientType gtype) 
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
   }
  
  drfp->SetUseGradientType(gtype);
}


template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::AllocateUpdateBuffer()
{
    std::cout << "In itkDiffeomorphicDemonsRegistrationFilter::AllocateUpdateBuffer" << std::endl;
  // The update buffer looks just like the output.
  DeformationFieldPointer output = this->GetOutput();
  DeformationFieldPointer upbuf = this->GetUpdateBuffer();

  upbuf->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  upbuf->SetRequestedRegion(output->GetRequestedRegion());
  upbuf->SetBufferedRegion(output->GetBufferedRegion());
  upbuf->SetSpacing(output->GetSpacing());
  upbuf->SetOrigin(output->GetOrigin());
  upbuf->Allocate();
    std::cout << "Exiting itkDiffeomorphicDemonsRegistrationFilter::AllocateUpdateBuffer" << std::endl;
}


/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::ApplyUpdate(TimeStepType dt)
{
  std::cout << "Entering ApplyUpdate" << std::endl;
   this->GetUpdateBuffer()->Modified();
   
  // If we smooth the update buffer before applying it, then the are
  // approximating a viscuous problem as opposed to an elastic problem
  if ( this->GetSmoothUpdateField() )
    {
    this->SmoothUpdateField();
    }

  /* use time step if necessary */
  if ( fabs(dt - 1.0)>1.0e-4 )
  {
     std::cout<<"Using timestep: "<<dt<<std::endl;
     m_Multiplier->SetConstant( dt );
     m_Multiplier->SetInput( this->GetUpdateBuffer() );
     m_Multiplier->GraftOutput( this->GetUpdateBuffer() );
     // in place update //m_Multiplier->UpdateLargestPossibleRegion();
     m_Multiplier->Update();
     this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
  }


  /* Configure the exponentiator */
  const double imposedMaxUpStep = this->GetMaximumUpdateStepLength();
  if ( imposedMaxUpStep > 0.0 )
  {
     /* max(norm(Phi))/2^N < 0.25*pixelspacing */
     const double numiterfloat = 2.0 + vcl_log(imposedMaxUpStep)/vnl_math::ln2;
     unsigned int numiter = 0;
     if ( numiterfloat > 0.0 )
        numiter = static_cast<unsigned int>( 1.0 + numiterfloat );
     
     m_Exponentiator->AutomaticNumberOfIterationsOff();
     m_Exponentiator->SetMaximumNumberOfIterations( numiter );
  }
  else
  {
     m_Exponentiator->AutomaticNumberOfIterationsOn();
     m_Exponentiator->SetMaximumNumberOfIterations( 2000u ); // just a high value
  }


  /*DEBUG*/
  // m_Multiplier->SetConstant( 14.0 );
  // m_Multiplier->SetInput( this->GetUpdateBuffer() );
  // m_Multiplier->UpdateLargestPossibleRegion();
  // this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
  /*END DEBUG*/


  /* Multiply the update buffer by -1 */ 
  m_Multiplier->SetConstant( -1 );
  m_Multiplier->SetInput( this->GetUpdateBuffer() ); 
  m_Multiplier->Update();// m_Multiplier->UpdateLargestPossibleRegion();


  /* Compose the vector fields to update the inverse warp field */
  m_Exponentiator->SetInput( m_Multiplier->GetOutput() );
  m_Warper->SetOutputSpacing( this->GetDeformationField()->GetSpacing() );
  m_Warper->SetOutputOrigin( this->GetDeformationField()->GetOrigin() );
  m_Warper->SetInput( this->GetInvDeformationField() );
  m_Warper->SetDeformationField( m_Exponentiator->GetOutput() );
  m_Adder->SetInput1( m_Warper->GetOutput() );
  m_Adder->SetInput2( m_Exponentiator->GetOutput() );
  m_Adder->GetOutput()->SetRequestedRegion( this->GetInvDeformationField()->GetRequestedRegion() );
  m_Adder->Update();
  this->SetInvDeformationField( m_Adder->GetOutput() );
  this->GetInvDeformationField()->DisconnectPipeline();


  /* Compose the vector fields to update the warp field */
  std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate: Composing the vector fields" << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
  m_Exponentiator->SetInput( this->GetUpdateBuffer() );
  m_Warper->SetOutputSpacing( this->GetUpdateBuffer()->GetSpacing() );
  m_Warper->SetOutputOrigin( this->GetUpdateBuffer()->GetOrigin() );
  m_Warper->SetInput( this->GetOutput() );
  m_Warper->SetDeformationField( m_Exponentiator->GetOutput() );
  m_Adder->SetInput1( m_Warper->GetOutput() );
  m_Adder->SetInput2( m_Exponentiator->GetOutput() );
  //m_Adder->UpdateLargestPossibleRegion();
  m_Adder->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
  m_Adder->Update();
  //std::cout<<"out buff spac: "<<this->GetOutput()->GetSpacing()<<std::endl;
  //std::cout<<"up buff spac: "<<this->GetUpdateBuffer()->GetSpacing()<<std::endl;
  //std::cout<<"exp out spac: "<<m_Exponentiator->GetOutput()->GetSpacing()<<std::endl;
  //std::cout<<"warp out spac: "<<m_Warper->GetOutput()->GetSpacing()<<std::endl;
  //std::cout<<"add out spac: "<<m_Adder->GetOutput()->GetSpacing()<<std::endl;
  

  /* Before we update the output we use that buffer to smooth the inverse warp field */
  if ( this->GetSmoothDeformationField() )
  {
    this->GraftOutput( this->GetInvDeformationField() );
    this->SmoothDeformationField();
    this->SetInvDeformationField( this->GetOutput() );
    this->GetInvDeformationField()->DisconnectPipeline();
  }


  /* Update the output with the new warp field */
  this->GraftOutput( m_Adder->GetOutput() );


  /* Set the RMS change */
  DemonsRegistrationFunctionType *drfp = dynamic_cast<DemonsRegistrationFunctionType *> (this->GetDifferenceFunction().GetPointer());
  if( !drfp )
   {
   itkExceptionMacro( << "Could not cast difference function to DemonsRegistrationFunction" );
   }
  this->SetRMSChange( drfp->GetRMSChange() );


  /* Smooth the output (the warp field) */
  if ( this->GetSmoothDeformationField() )
  {
     this->SmoothDeformationField();
  }

}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
DiffeomorphicDemonsRegistrationFilter<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{ 
  Superclass::PrintSelf( os, indent );
  os << indent << "Intensity difference threshold: " <<
    this->GetIntensityDifferenceThreshold() << std::endl;
}


} // end namespace itk

#endif

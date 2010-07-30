#ifndef _itkDiffeomorphicDemonsRegistrationFilter_txx
#define _itkDiffeomorphicDemonsRegistrationFilter_txx
#include "itkDiffeomorphicDemonsRegistrationFilter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

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

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfp.GetPointer() ) );

  m_Multiplier = MultiplyByConstantType::New();
  // m_Multiplier->InPlaceOn();

  m_Exponentiator = FieldExponentiatorType::New();
  
  m_Warper = VectorWarperType::New();
  FieldInterpolatorPointer VectorInterpolator =
     FieldInterpolatorType::New();
  m_Warper->SetInterpolator(VectorInterpolator);

  m_Adder = AdderType::New();
  m_Adder->InPlaceOn();

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

  f->SetDeformationField( this->GetDeformationField() );
  f->SetInvDeformationField( this->GetDeformationField() );

  // if (this->GetInvDeformationField() == NULL)
  //   f->SetInvDeformationField( this->GetDeformationField() );
  // else
  // {
  //   f->SetInvDeformationField( this->GetDeformationField() );
  //   // f->SetInvDeformationField( this->GetInvDeformationField() );
  //   // std::cout << "inv deformation origin: " << this->GetInvDeformationField()->GetOrigin() << std::endl;
  //   // std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
  //   // std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
  //   // std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
  // }

  // call the superclass  implementation ( initializes f )
  Superclass::InitializeIteration();
  if (this->GetInvDeformationField() != NULL)
  {
    std::cout << "---------------------" << std::endl;
    std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration(): " << std::endl;
    std::cout << "After superclass InitializeIteration call:" << std::endl;
    std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "inv deformation size: " << this->GetInvDeformationField()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    std::cout << "---------------------" << std::endl;
  }
    std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration: output size: " << this->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration: output size: " << this->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration: output size: " << this->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    std::cout << "DiffeomorphicDemonsRegistrationFilter:InitializeIteration(): leaving" << std::endl;
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
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate" << std::endl;
   this->GetUpdateBuffer()->Modified();
   
  // If we smooth the update buffer before applying it, then the are
  // approximating a viscuous problem as opposed to an elastic problem
  if ( this->GetSmoothUpdateField() )
    {
    this->SmoothUpdateField();
    }

  // use time step if necessary
  if ( fabs(dt - 1.0)>1.0e-4 )
  {
     std::cout<<"Using timestep: "<<dt<<std::endl;
     m_Multiplier->SetConstant( dt );
     m_Multiplier->SetInput( this->GetUpdateBuffer() );
     m_Multiplier->GraftOutput( this->GetUpdateBuffer() );
     // in place update
     //m_Multiplier->UpdateLargestPossibleRegion();
     m_Multiplier->Update();
     // graft output back to this->GetUpdateBuffer()
     this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
  }



  // compute the exponential
  m_Exponentiator->SetInput( this->GetUpdateBuffer() );
  const double imposedMaxUpStep = this->GetMaximumUpdateStepLength();
  if ( imposedMaxUpStep > 0.0 )
  {
     // max(norm(Phi))/2^N < 0.25*pixelspacing
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


  std::cout << "In itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate: Going to warp/compose the vector fields" << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[0] << std::endl;

  // compose the vector fields
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
  
  // Region passing stuff
  this->GraftOutput( m_Adder->GetOutput() );

  

  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
   {
   itkExceptionMacro( << 
     "Could not cast difference function to DemonsRegistrationFunction" );
   }

  this->SetRMSChange( drfp->GetRMSChange() );

  /*
   * Smooth the deformation field
   */
  if ( this->GetSmoothDeformationField() )
  {
     this->SmoothDeformationField();

     //double var[DeformationFieldType::ImageDimension];
     //for (unsigned int i=0; i<DeformationFieldType::ImageDimension; ++i)
     //   var[i] = vnl_math_sqr( this->GetUpdateFieldStandardDeviations()[i] );
        
     //typedef SmoothingRecursiveGaussianImageFilter<DeformationFieldType,DeformationFieldType> GaussianFilterType;
     //typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
     //smoother->SetInput( this->GetOutput() );
     //smoother->SetVariance( var );
     //smoother->SetSigma( this->GetUpdateFieldStandardDeviations()[0] );
     //smoother->Update();

     //this->GraftOutput( smoother->GetOutput() );

  }

  // Update the inverse deformation field
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate: update buffer size: " << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate: update buffer size: " << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate: update buffer size: " << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    m_Multiplier->SetConstant( 1 );
    m_Multiplier->SetInput( this->GetUpdateBuffer() );
    // m_Multiplier->SetInput( this->GetDeformationField() );
    // m_Multiplier->UpdateLargestPossibleRegion();
    m_Multiplier->Update();
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate:multiplier output size: " << m_Multiplier->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate:multiplier output size: " << m_Multiplier->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "itkDiffeomorphicDemonsRegistrationFilter::ApplyUpdate:multiplier output size: " << m_Multiplier->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
  // m_Exponentiator->SetInput( m_Multiplier->GetOutput() );
  // if (this->GetInvDeformationField() == NULL)
  // {
  //   // std::cout << "apply update: Inv deformation is null" << std::endl;
  //   // this->SetInvDeformationField(this->GetDeformationField());
  //   this->SetInvDeformationField(m_Exponentiator->GetOutput());
  // }
  // else
  // {
  // // m_Warper->SetOutputSpacing( this->GetDeformationField()->GetSpacing() );
  // // m_Warper->SetOutputOrigin( this->GetDeformationField()->GetOrigin() );
  //   // std::cout << "in here" << std::endl;
  //   // m_Warper->SetInput( this->GetInvDeformationField() );
  //   // m_Warper->SetDeformationField( m_Exponentiator->GetOutput() );
  //   // m_Adder->SetInput1( m_Warper->GetOutput() );
  //   // m_Adder->SetInput2( m_Exponentiator->GetOutput() );
  //   // this->SetInvDeformationField( m_Adder->GetOutput() );

  // }

  // m_Adder->GetOutput()->SetRequestedRegion( this->GetInvDeformationField()->GetRequestedRegion() );
  //
  // m_Adder->Update();

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

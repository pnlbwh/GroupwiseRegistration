/*=========================================================================
 
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkESMInvConDemonsRegistrationFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-07-11 19:02:04 $
  Version:   $Revision: 1.4 $
 
  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.
 
=========================================================================*/

#ifndef __itkESMInvConDemonsRegistrationFunction_txx
#define __itkESMInvConDemonsRegistrationFunction_txx

#include "itkESMInvConDemonsRegistrationFunction.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_trace.h"
#include "vnl/algo/vnl_determinant.h"

namespace itk {
    
/**
 * Default constructor
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ESMInvConDemonsRegistrationFunction()
{
    
    RadiusType r;
    unsigned int j;
    for( j = 0; j < ImageDimension; j++ )
    {
        r[j] = 0;
    }
    this->SetRadius(r);
    
    m_TimeStep = 1.0;
    m_DenominatorThreshold = 1e-9;
    m_IntensityDifferenceThreshold = 0.001;
    m_MaximumUpdateStepLength = 0.5;
    m_RegWeight = 1000;
    
    this->SetMovingImage(NULL);
    this->SetFixedImage(NULL);
    this->SetInvDeformationField(NULL);
    this->SetJacobianDetImage(NULL);
    this->SetFwWeightImage(NULL);
    
    m_FixedImageSpacing.Fill( 1.0 );
    m_FixedImageOrigin.Fill( 0.0 );
    m_Normalizer = 0.0;
    m_FixedImageGradientCalculator = GradientCalculatorType::New();
    m_MappedMovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
    
    this->m_UseGradientType = Symmetric;
    
    typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();
    
    m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
    
    m_MovingImageWarper = WarperType::New();
    m_MovingImageWarper->SetInterpolator( m_MovingImageInterpolator );
    m_MovingImageWarper->SetEdgePaddingValue( NumericTraits<MovingPixelType>::max() );
    
    typename DefaultInterpolatorType::Pointer interp_fixed =
    DefaultInterpolatorType::New();
    
    m_FixedImageInterpolator = static_cast<InterpolatorType*>(
    interp_fixed.GetPointer() );
    
    m_FixedImageWarper = WarperType::New();
    m_FixedImageWarper->SetInterpolator( m_FixedImageInterpolator );
    m_FixedImageWarper->SetEdgePaddingValue( NumericTraits<FixedPixelType>::max() );
    
    m_Metric = NumericTraits<double>::max();
    m_SumOfSquaredDifference = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_RMSChange = NumericTraits<double>::max();
    m_SumOfSquaredChange = 0.0;
    m_UseJacobian = false;
    m_UseFwWeight = false;
    }
    
    
/*
 * Standard "PrintSelf" method.
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    
    os << indent << "UseGradientType: ";
    os << m_UseGradientType << std::endl;
    os << indent << "MaximumUpdateStepLength: ";
    os << m_MaximumUpdateStepLength << std::endl;
    
    os << indent << "MovingImageIterpolator: ";
    os << m_MovingImageInterpolator.GetPointer() << std::endl;
    os << indent << "FixedImageGradientCalculator: ";
    os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
    os << indent << "MappedMovingImageGradientCalculator: ";
    os << m_MappedMovingImageGradientCalculator.GetPointer() << std::endl;
    os << indent << "DenominatorThreshold: ";
    os << m_DenominatorThreshold << std::endl;
    os << indent << "IntensityDifferenceThreshold: ";
    os << m_IntensityDifferenceThreshold << std::endl;
    
    os << indent << "Metric: ";
    os << m_Metric << std::endl;
    os << indent << "SumOfSquaredDifference: ";
    os << m_SumOfSquaredDifference << std::endl;
    os << indent << "NumberOfPixelsProcessed: ";
    os << m_NumberOfPixelsProcessed << std::endl;
    os << indent << "RMSChange: ";
    os << m_RMSChange << std::endl;
    os << indent << "SumOfSquaredChange: ";
    os << m_SumOfSquaredChange << std::endl;
    
    }
    
/**
 *
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SetIntensityDifferenceThreshold(double threshold)
{
    m_IntensityDifferenceThreshold = threshold;
    }
    
/**
 *
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    double
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::GetIntensityDifferenceThreshold() const
{
    return m_IntensityDifferenceThreshold;
    }
    
/**
 * Set the function state values before each iteration
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::InitializeIteration()
{
    if( !this->GetMovingImage() || !this->GetFixedImage()
    || !m_MovingImageInterpolator )
    {
        itkExceptionMacro(
        << "MovingImage, FixedImage and/or Interpolator not set" );
    }
    
    // cache fixed image information
    m_FixedImageSpacing = this->GetFixedImage()->GetSpacing();
    m_FixedImageOrigin  = this->GetFixedImage()->GetOrigin();
    
    // compute the normalizer
    if( m_MaximumUpdateStepLength > 0.0 )
    {
        m_Normalizer = 0.0;
        for( unsigned int k = 0; k < ImageDimension; k++ )
        {
            m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
        }
        m_Normalizer *= m_MaximumUpdateStepLength * m_MaximumUpdateStepLength /
        static_cast<double>( ImageDimension );
    }
    else
    {
        // set it to minus one to denote a special case
        // ( unrestricted update length )
        m_Normalizer = -1.0;
    }
    
    /*if ( this->m_UseJacobian && !this->GetJacobianDetImage()){
      itkExceptionMacro(
        << "Jacobian Image needs to be set ! ");  
    }*/
    // setup gradient calculator
    m_FixedImageGradientCalculator->SetInputImage( this->GetFixedImage() );
    m_MappedMovingImageGradientCalculator->SetInputImage( this->GetMovingImage() );
    
    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration " << std::endl;
    // Compute warped moving image
    m_MovingImageWarper->SetOutputSpacing( this->GetFixedImage()->GetSpacing() );
    m_MovingImageWarper->SetOutputOrigin( this->GetFixedImage()->GetOrigin() );
    m_MovingImageWarper->SetOutputDirection( this->GetFixedImage()->GetDirection() );
    m_MovingImageWarper->SetInput( this->GetMovingImage() );
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image size" << this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image size" << this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image size" << this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image requested size" << this->GetMovingImage()->GetRequestedRegion().GetSize()[0] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image requested size" << this->GetMovingImage()->GetRequestedRegion().GetSize()[1] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: moving image requested size" << this->GetMovingImage()->GetRequestedRegion().GetSize()[2] << std::endl;
    m_MovingImageWarper->SetDeformationField( this->GetDeformationField() );
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: deformation size" << this->GetDeformationField()->GetLargestPossibleRegion().GetSize()[0] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: deformation size" << this->GetDeformationField()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: deformation size" << this->GetDeformationField()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: requested deformation size" << this->GetDeformationField()->GetRequestedRegion().GetSize()[0] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: requested deformation size" << this->GetDeformationField()->GetRequestedRegion().GetSize()[1] << std::endl;
    std::cout << "itkESMInvConDemonsRegistrationFunction::InitializeIteration: requested deformation size" << this->GetDeformationField()->GetRequestedRegion().GetSize()[2] << std::endl;
    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration: about to set region for moving image warper" << std::endl;
    m_MovingImageWarper->GetOutput()->SetRequestedRegion( this->GetDeformationField()->GetRequestedRegion() );
    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration: about to warp moving image" << std::endl;
    m_MovingImageWarper->Update();
    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration:after warped moving image " << std::endl;
    
    // Compute warped moving image
    m_FixedImageWarper->SetOutputSpacing( this->GetFixedImage()->GetSpacing() );
    m_FixedImageWarper->SetOutputOrigin( this->GetFixedImage()->GetOrigin() );
    m_FixedImageWarper->SetOutputDirection( this->GetFixedImage()->GetDirection() );
    m_FixedImageWarper->SetInput( this->GetFixedImage() );
    m_FixedImageWarper->SetDeformationField( this->GetInvDeformationField() );
    m_FixedImageWarper->GetOutput()->SetRequestedRegion( this->GetInvDeformationField()->GetRequestedRegion() );
    m_FixedImageWarper->Update();
    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration:after warped fixed image " << std::endl;
    
    // setup moving image interpolator for further access
    m_MovingImageInterpolator->SetInputImage( this->GetMovingImage() );
    m_FixedImageInterpolator->SetInputImage( this->GetFixedImage() );
    
    // initialize metric computation variables
    m_SumOfSquaredDifference  = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_SumOfSquaredChange      = 0.0;

    std::cout << "In itkESMInvConDemonsRegistrationFunction:InitializeIteration: leaving" << std::endl;
}
    
    
/**
 * Compute update at a non boundary neighbourhood
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    typename ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::PixelType
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ComputeUpdate(const NeighborhoodType &it, void * gd,
    const FloatOffsetType& itkNotUsed(offset))
{
    // std::cout << "In itkESMInvConDemonsRegistrationFunction:ComputeUpdate: start" << std::endl;
    
    GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
    PixelType update;
    IndexType FirstIndex = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex();
    IndexType LastIndex = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex() +
    this->GetFixedImage()->GetLargestPossibleRegion().GetSize();
    
    const IndexType index = it.GetIndex();
    
    // Get fixed image related information
    // Note: no need to check if the index is within
    // fixed image buffer. This is done by the external filter.
    const double fixedValue = static_cast<double>(this->GetFixedImage()->GetPixel( index ) );
    double jacobian = 1.0f;
    double fw_weight = 1.0f;
    if (m_UseJacobian)
    {
        jacobian = static_cast<double>(this->GetJacobianDetImage()->GetPixel(index));
    }
    
    if (m_UseFwWeight)
    {
        fw_weight = static_cast<double>(this->GetFwWeightImage()->GetPixel(index));
    }
    MovingPixelType warpedFixedPixValue
    = m_FixedImageWarper->GetOutput()->GetPixel( index );
    
       
    if( warpedFixedPixValue == NumericTraits <FixedPixelType>::max() )
    {
        update.Fill( 0.0 );
        return update;
    }
    double warpedFixedValue = static_cast<double>( warpedFixedPixValue );
    
    const double movingValue = static_cast<double>(this->GetMovingImage()->GetPixel( index ) );
    // Get moving image related information
    // check if the point was mapped outside of the moving image using
    // the "special value" NumericTraits<MovingPixelType>::max()
    
    MovingPixelType warpedMovingPixValue
    = m_MovingImageWarper->GetOutput()->GetPixel( index );
    
    if( warpedMovingPixValue == NumericTraits <MovingPixelType>::max() )
    {
        update.Fill( 0.0 );
        return update;
    }
    
    
    double warpedMovingValue = static_cast<double>( warpedMovingPixValue );
    //const double movingValue = static_cast<double>( movingPixValue );
    
    CovariantVectorType usedGradientTimes2;
    
    // we don't use a CentralDifferenceImageFunction here to be able to
    // check for NumericTraits<MovingPixelType>::max()
    CovariantVectorType warpedMovingGradient;
    IndexType tmpIndex = index;
    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
        // bounds checking
        if( FirstIndex[dim]==LastIndex[dim]  || index[dim]<FirstIndex[dim] || index[dim]>=LastIndex[dim] )
        {
            warpedMovingGradient[dim] = 0.0;
            continue;
        }
        else if ( index[dim] == FirstIndex[dim] )
        {
            // compute derivative
            tmpIndex[dim] += 1;
            warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits <MovingPixelType>::max() )
            {
                // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // forward difference
                warpedMovingGradient[dim] = static_cast<double>( warpedMovingPixValue ) - warpedMovingValue;
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex[dim] -= 1;
            continue;
        }
        else if ( index[dim] == (LastIndex[dim]-1) )
        {
            // compute derivative
            tmpIndex[dim] -= 1;
            warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
                // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedMovingGradient[dim] = warpedMovingValue - static_cast<double>( warpedMovingPixValue );
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex[dim] += 1;
            continue;
        }
        
        
        // compute derivative
        tmpIndex[dim] += 1;
        warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
        if ( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
        {
            // backward difference
            warpedMovingGradient[dim] = warpedMovingValue;
    
            tmpIndex[dim] -= 2;
            warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
            // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedMovingGradient[dim] -= static_cast<double>(
                m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex ) );
        
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
        }
        else
        {
            warpedMovingGradient[dim] = static_cast<double>( warpedMovingPixValue );
            
            tmpIndex[dim] -= 2;
            warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            if ( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
                // forward difference
                warpedMovingGradient[dim] -= warpedMovingValue;
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            else
            {
                // normal case, central difference
                warpedMovingGradient[dim] -= static_cast<double>( warpedMovingPixValue );
                warpedMovingGradient[dim] *= 0.5 / m_FixedImageSpacing[dim];
            }
        }
        tmpIndex[dim] += 1;
    }
    
    CovariantVectorType warpedFixedGradient;
    IndexType tmpIndex2 = index;
    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
        // bounds checking
        if( FirstIndex[dim]==LastIndex[dim] || index[dim]<FirstIndex[dim] || index[dim]>=LastIndex[dim] )
        {
            warpedFixedGradient[dim] = 0.0;
            continue;
        }
        else if ( index[dim] == FirstIndex[dim] )
        {
            // compute derivative
            tmpIndex2[dim] += 1;
            warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits <FixedPixelType>::max() )
            {
                // weird crunched border case
                warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // forward difference
                warpedFixedGradient[dim] = static_cast<double>( warpedFixedPixValue ) - warpedFixedValue;
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex2[dim] -= 1;
            continue;
        }
        else if ( index[dim] == (LastIndex[dim]-1) )
        {
            // compute derivative
            tmpIndex2[dim] -= 1;
            warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // weird crunched border case
                warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedFixedGradient[dim] = warpedFixedValue - static_cast<double>( warpedFixedPixValue );
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex2[dim] += 1;
            continue;
        }
        
        
        // compute derivative
        tmpIndex2[dim] += 1;
        warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
        if ( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
        {
            // backward difference
            warpedFixedGradient[dim] = warpedFixedValue;
    
            tmpIndex2[dim] -= 2;
            warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // weird crunched border case
                 warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedFixedGradient[dim] -= static_cast<double>(
                m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 ) );
        
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
        }
        else
        {
            warpedFixedGradient[dim] = static_cast<double>( warpedFixedPixValue );
            
            tmpIndex2[dim] -= 2;
            warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            if ( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // forward difference
                warpedFixedGradient[dim] -= warpedFixedValue;
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            else
            {
                // normal case, central difference
                warpedFixedGradient[dim] -= static_cast<double>( warpedFixedPixValue );
                warpedFixedGradient[dim] *= 0.5 / m_FixedImageSpacing[dim];
            }
        }
        tmpIndex2[dim] += 1;
    }
    
    vnl_vector<double> warpedFixedGradient_vec(ImageDimension);
    vnl_vector<double> warpedMovingGradient_vec(ImageDimension);
    
    vnl_matrix<double> Hess(ImageDimension, ImageDimension);
    vnl_matrix<double> invHess(ImageDimension, ImageDimension);
    
    
    const double speedValue = fixedValue - warpedMovingValue;
    const double speedValue2 = warpedFixedValue - movingValue;
    
    //std::cout << "Reg Weight (inside) " << this->m_RegWeight << std::endl;
    
    vnl_diag_matrix<double> Reg(ImageDimension, this->m_RegWeight);
    warpedFixedGradient_vec = warpedFixedGradient.GetVnlVector();
    warpedMovingGradient_vec = warpedMovingGradient.GetVnlVector();
    
     
    //std::cout << "Jacobian " << jacobian << "Fw weight: " << fw_weight << std::endl;
    
    Hess = 0.25*fw_weight*outer_product(warpedMovingGradient_vec, warpedMovingGradient_vec) + 
        jacobian * 0.25*outer_product(warpedFixedGradient_vec, warpedFixedGradient_vec);
    
    //Hess = 0.25*outer_product(warpedMovingGradient_vec, warpedMovingGradient_vec) + 
    //    0.25*outer_product(warpedFixedGradient_vec, warpedFixedGradient_vec);
    
    if (jacobian == 0)
    {
       Hess = 4.0*Hess;
    }
    
    invHess = vnl_inverse(Hess + Reg);
    vnl_matrix<double> temp_mat(ImageDimension,1);
    
    temp_mat.set_column(0,0.25*fw_weight*(fixedValue - warpedMovingValue)*warpedMovingGradient_vec + 
            jacobian * 0.25*(warpedFixedValue - movingValue)*warpedFixedGradient_vec);
     
    if (jacobian == 0)
    {
        temp_mat = 4.0*temp_mat;
    }
    
    vnl_matrix<double> update_mat(ImageDimension,1);
    
    update_mat = invHess * temp_mat;
    
           
  /**
   * Compute Update.
   * We avoid the mismatch in units between the two terms.
   * and avoid large step using a normalization term.
   */
        
    if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold || vnl_math_abs(speedValue2) < m_IntensityDifferenceThreshold)
    {
        update.Fill( 0.0 );
    }
    else
    {
        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
            update[dim] = update_mat[dim][0];
            
        }
    }
    
        
    // WARNING!! We compute the global data without taking into account the current update step.
    // There are several reasons for that: If an exponential, a smoothing or any other operation
    // is applied on the update field, we cannot compute the newMappedCenterPoint here; and even
    // if we could, this would be an often unnecessary time-consuming task.
    if ( globalData )
    {
        globalData->m_SumOfSquaredDifference += vnl_math_sqr( speedValue );
        globalData->m_NumberOfPixelsProcessed += 1;
        globalData->m_SumOfSquaredChange += update.GetSquaredNorm();
    }
    
    return update;
}
    
    
/**
 * Update the metric and release the per-thread-global data.
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ReleaseGlobalDataPointer( void *gd ) const
{
    GlobalDataStruct * globalData = (GlobalDataStruct *) gd;
    
    m_MetricCalculationLock.Lock();
    m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
    m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
    m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
    if( m_NumberOfPixelsProcessed )
    {
        m_Metric = m_SumOfSquaredDifference /
        static_cast<double>( m_NumberOfPixelsProcessed );
        m_RMSChange = vcl_sqrt( m_SumOfSquaredChange /
        static_cast<double>( m_NumberOfPixelsProcessed ) );
    }
    m_MetricCalculationLock.Unlock();
    
    delete globalData;
    }
    
} // end namespace itk

#endif

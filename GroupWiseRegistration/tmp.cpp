typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  typename ImageType::Pointer template_vol = 0;
  typename ImageType::Pointer template_vol2 = 0;

  for (int i=0; i < args.volumeFileNames.size(); i++)
  {
    std::cout << args.volumeFileNames[i] << std::endl;

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

    if (i==0)
    {
      imageReader2->SetFileName( args.volumeFileNames[0] );
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
      template_vol = imageReader2->GetOutput();
      std::cout << "setting intitial template_vol to " << args.volumeFileNames[i] << std::endl;
    }
    else
    {
      adder->SetInput1( template_vol );
      adder->SetInput2( imageReader->GetOutput() );
      // adder->InPlaceOn();
      adder->Update();
      template_vol = adder->GetOutput();
    }
  }
  typedef itk::MultiplyByConstantImageFilter< ImageType, float, ImageType > MultiplyByConstantType;
  typedef typename MultiplyByConstantType::Pointer   MultiplyByConstantPointer;
  MultiplyByConstantPointer m_Multiplier = MultiplyByConstantType::New();
  m_Multiplier->SetConstant( 1.0/args.volumeFileNames.size() );
  m_Multiplier->SetInput( template_vol );
  m_Multiplier->Update();
  template_vol = m_Multiplier->GetOutput();


// typename TDeformationField::RegionType     region;
    // typename TDeformationField::IndexType      start;
    // // region.SetSize( this->GetOutput()->GetLargestPossibleRegion().GetSize() );
    // region.SetSize( this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize() );
    // // start.Fill(0.0);
    // // region.SetIndex( start );
    // m_logDeformationField->SetDirection( this->GetUpdateBuffer()->GetDirection() );
    // m_logDeformationField->SetOrigin( this->GetUpdateBuffer()->GetOrigin() );
    // m_logDeformationField->SetSpacing( this->GetUpdateBuffer()->GetSpacing());
    // m_logDeformationField->SetRegions( region );
    // m_logDeformationField->Allocate();
    // m_logDeformationField->FillBuffer( 0.0 );


/*
   * Update the log field
   */
  typedef typename itk::ImageFileReader< TDeformationField > DeformationReaderType;
  typename DeformationReaderType::Pointer                    reader =  DeformationReaderType::New();
  typedef itk::ImageFileWriter< TDeformationField >  DeformationWriterType;
  typename DeformationWriterType::Pointer            writer =  DeformationWriterType::New();

  reader->SetFileName( "log_field.nii.gz" );
  try
  {
    reader->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Could not read one of the input images." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }
  m_logDeformationField = reader->GetOutput();
  m_Adder->SetInput1( m_logDeformationField );
  m_Adder->SetInput2( this->GetUpdateBuffer() );
  // m_Adder->GetOutput()->SetRequestedRegion( this->GetDeformationField()->GetRequestedRegion() );
  // m_Adder->GetOutput()->SetRequestedRegion( this->GetUpdateBuffer()->GetRequestedRegion() );
  // if ( this->GetSmoothDeformationField() )
  // {
  //   this->GraftOutput( m_Adder->GetOutput() );
  //   this->SmoothDeformationField();
  // }
  writer->SetFileName( "log_field.nii.gz" );
  writer->SetUseCompression( true );
  // writer->SetInput( this->GetOutput() );
  writer->SetInput( m_Adder->GetOutput() );
std::cout<<"log field image size: "<<m_logDeformationField->GetLargestPossibleRegion().GetSize()[0] << " " <<
         m_logDeformationField->GetLargestPossibleRegion().GetSize()[1] << " " << m_logDeformationField->GetLargestPossibleRegion().GetSize()[2] << std::endl;
std::cout<<"update buffer image size: "<<this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[0] << " " <<
         this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[1] << " " << this->GetUpdateBuffer()->GetLargestPossibleRegion().GetSize()[2] << std::endl;

  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Could not write log field to disk." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }


/*
   * Set the inverse deformation field
   */
  if (isFirstIteration)
  {
    typedef itk::ImageFileWriter< TDeformationField >  DeformationWriterType;
    typename DeformationWriterType::Pointer            writer =  DeformationWriterType::New();
    writer->SetFileName( "log_field.nii.gz" );
    writer->SetUseCompression( true );
    writer->SetInput( this->GetDeformationField() );  //should be zeros
    try
    {
      writer->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not write log field to disk." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }
    f->SetInvDeformationField( this->GetDeformationField() );  //should be zeros
    isFirstIteration = false;
  }
  else
  {
    typedef typename itk::ImageFileReader< TDeformationField > DeformationReaderType;
    typename DeformationReaderType::Pointer      reader =  DeformationReaderType::New();
    reader->SetFileName( "log_field.nii.gz" );
    try
    {
      reader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cout << "Could not read the log field." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
    }
      std::cout << "Setting the function's inverse deformation field" << std::endl;
    m_Multiplier->SetConstant( -1 );
    m_Multiplier->SetInput( reader->GetOutput() );
    m_Exponentiator->SetInput( m_Multiplier->GetOutput() );
    m_Exponentiator->Update();
    f->SetInvDeformationField( m_Exponentiator->GetOutput() );  //should be zeros
  }

//DEBUG
  // if (this->GetInvDeformationField() == NULL)
  // {
  //   this->SetInvDeformationField( this->GetDeformationField() );
  //   // f->SetInvDeformationField( this->GetInvDeformationField() );  //should be zeros

  //   if (this->GetDeformationField() == NULL)
  //     std::cout << "  InitializeIteration: GetDeformationField is NULL!" << std::endl;

  //   std::cout << "  ************************************************" << std::endl;
  //   std::cout << "  InitializeIteration: GetInvDeformationField is NULL" << std::endl;
  //   std::cout << "  ************************************************" << std::endl;
  //   /*Debug*/
  //   try
  //   {
  //     typedef itk::ImageFileWriter< DeformationFieldType>  FieldWriterType;
  //     typename FieldWriterType::Pointer      fieldwriter =  FieldWriterType::New();
  //     fieldwriter->SetUseCompression( true );
  //     fieldwriter->SetFileName( "invdeformation-initial.nii.gz" );
  //     fieldwriter->SetInput( this->GetInvDeformationField()  );
  //     fieldwriter->Update();
  //   }
  //   catch( itk::ExceptionObject& err )
  //   {
  //     std::cout << "Unexpected error." << std::endl;
  //     std::cout << err << std::endl;
  //     exit( EXIT_FAILURE );
  //   }
  //   // exit(0);
  //   /*End Debug*/

  // }
  // else
  // {
  //   f->SetInvDeformationField( this->GetInvDeformationField() );
  //   /*Debug*/
  //   try
  //   {
  //     typedef itk::ImageFileWriter< DeformationFieldType>  FieldWriterType;
  //     typename FieldWriterType::Pointer      fieldwriter =  FieldWriterType::New();
  //     fieldwriter->SetUseCompression( true );
  //     fieldwriter->SetFileName( "invdeformation.nii.gz" );
  //     fieldwriter->SetInput( this->GetInvDeformationField()  );
  //     fieldwriter->Update();
  //     fieldwriter->SetFileName( "deformation.nii.gz" );
  //     fieldwriter->SetInput( this->GetDeformationField()  );
  //     fieldwriter->Update();
  //   }
  //   catch( itk::ExceptionObject& err )
  //   {
  //     std::cout << "Unexpected error." << std::endl;
  //     std::cout << err << std::endl;
  //     exit( EXIT_FAILURE );
  //   }
  //   /*End Debug*/

  // }
  // f->SetDeformationField( this->GetDeformationField() );


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

 //template_vol = imageReader2->GetOutput();
    //template_vol->DisconnectPipeline();
    
    //ImageReaderType::Pointer imageReader2 = ImageReaderType::New();
    //imageReader2->SetFileName(template_name.str());
    //try
    //{
      //imageReader2->Update();
    //}
    //catch( itk::ExceptionObject& err )
    //{
      //std::cout << "Could not read one of the input images." << std::endl;
      //std::cout << err << std::endl;
      //exit( EXIT_FAILURE );
    //}
    //ImageType::Pointer                  template_vol = 0;
    //template_vol = imageReader2->GetOutput();
    //template_vol->DisconnectPipeline();

  DemonsRegistrationFunctionType::GradientType   gtype;
//if ( args.initialSigmaDiff > 0.1 ) // if ( args.sigmaDef > 0.1 )
      //{
        //filter->SmoothDeformationFieldOn();
        //filter->SetStandardDeviations( args.initialSigmaDiff );  //TODO Update with a range
      //}
      //else
        //filter->SmoothDeformationFieldOff();

      //if ( args.sigmaUp > 0.1 )
      //{
        //filter->SmoothUpdateFieldOn();
        //filter->SetUpdateFieldStandardDeviations( args.sigmaUp );
      //}
      //else
        //filter->SmoothUpdateFieldOff();

      // filter->SetInitialDeformationField( deformationReader->GetOutput() );


 DeformationWriterType::Pointer  defWriter =  DeformationWriterType::New();
      defWriter->SetFileName( deformation_name.str() );
      defWriter->SetInput( field  );
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

defWriter->SetFileName( deformation_name.str() );
      defWriter->SetInput( field  );
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
      /* Save the warp, e.g. image1_0_deformation.nii.gz */
      deformation_name.str("");
      deformation_name << args.volumeFileNames[i] << "_" << j << "_" "deformation.nii.gz";
      WriteDeformationField(field, deformation_name.str());

      /* DEBUG: Save the warp as, e.g. image1_0_deformation-orig.nii.gz */
      deformation_name.str("");
      deformation_name << args.volumeFileNames[i] << "_" << j << "_" "deformation-orig.nii.gz";


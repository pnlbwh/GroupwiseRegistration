#include "itkESMInvConDemonsRegistrationFunction.h"

template<unsigned int Dimension>
class GroupWiseRegistration
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  //typedef itk::OrientedImage< PixelType, Dimension >  OrientedImageType;// read and transform oriented images, operate on non-oriented
  typedef float VectorComponentType;
  typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
  typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

  typedef typename itk::ImageFileReader< ImageType > ReaderType;
  typedef typename itk::ImageFileReader< DeformationFieldType > DeformationReaderType;

  typedef itk::ImageFileWriter< DeformationFieldType >  DeformationWriterType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  typedef itk::ESMInvConDemonsRegistrationFunction <ImageType, ImageType, DeformationFieldType>  DemonsRegistrationFunctionType;


  void runGroupWiseRegistration();

private:
  std::string outputImageFile;
  unsigned int numLevels;
  unsigned int numIterations;
  std::vector<std::string> volumeFileNames;   
  std::string resultsDirectory;
  unsigned int numOuterIterations;
  float initialSigmaDiff;
  float finalSigmaDiff;
  float regWeight;
  bool  useJac;
  float sigmaDef;        
  float sigmaUp;          
  float maxStepLength;     
  unsigned int verbosity;

  typename WriterType::Pointer      m_writer;
  typename ReaderType::Pointer      m_reader;

};

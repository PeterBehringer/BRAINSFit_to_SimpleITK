import os
import SimpleITK as sitk

%run TestFilesCommonCode.ipynb

def command_iteration(method) :
 print("{0:3} = {1:10.5f} : {2}".format(method.GetOptimizerIteration(),method.GetMetricValue(),method.GetOptimizerPosition()))

# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax.nrrd'

fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)

# read input masks

fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'

fixedMask=sitk.ReadImage(fixedMaskFilename)
movingMask=sitk.ReadImage(movingMaskFilename)

# set output file

outputRegistrationVolume = '/Users/peterbehringer/MyTesting/outputRegistrationVolume_1.h5'


# perform Rigid Transformation

# First transform is passed through an initializer
tx1 = sitk.CenteredVersorTransformInitializer(fixedVolume, movingVolume, sitk.VersorRigid3DTransform())

ctx = sitk.Transform(tx1) # Set composite transform
                          # This composite transform is update at each stage
                          # and finally will be considered as the output of
                          # the registration process since the InPlace is TRUE.
ctx.SetFixedParameters(ctx.GetFixedParameters()) # hack to force deep copy, as registion is done in place..
# Set the registration filter
R = sitk.ImageRegistrationMethod()
R.SetMetricAsMattesMutualInformation( 200 )
R.SetOptimizerAsRegularStepGradientDescent(learningRate =0.5,
                                           minStep=1e-2,
                                           numberOfIterations=250,
                                           gradientMagnitudeTolerance=1e-4,
                                           estimateLearningRate=R.Never)
R.SetOptimizerScales([1, 1, 1, 1.0/250, 1.0/250, 1.0/250])
R.SetShrinkFactorsPerLevel([1])
R.SetSmoothingSigmasPerLevel([0])
R.SetInitialTransform(ctx)
R.SetInterpolator(sitk.sitkLinear)
R.SetMetricSamplingPercentage(0.5)
R.SetMetricSamplingStrategy(R.RANDOM)
R.RemoveAllCommands()
R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )
outTx = R.Execute(sitk.Cast(fixedVolume,sitk.sitkFloat32), sitk.Cast(movingVolume,sitk.sitkFloat32))

sitk.WriteTransform(outTx,outputRegistrationVolume)

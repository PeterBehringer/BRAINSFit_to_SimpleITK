import os
import SimpleITK as sitk


def command_iteration(method) :
 print("{0:3} = {1:10.5f} : {2}".format(method.GetOptimizerIteration(),method.GetMetricValue(),method.GetOptimizerPosition()))

# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax.nrrd'

fixedVolume=sitk.ReadImage(fixedImageFilename)
movingVolume=sitk.ReadImage(movingImageFilename)

# read input masks

fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'

fixedMask=sitk.ReadImage(fixedMaskFilename)
movingMask=sitk.ReadImage(movingMaskFilename)

# set output file

outputRegistrationVolume = '/Users/peterbehringer/MyTesting/outputRegistrationVolume.h5'

# create fixedVolume
tx1 = sitk.CenteredVersorTransformInitializer(fixedVolume, movingVolume, sitk.VersorRigid3DTransform())
ctx = sitk.Transform(tx1)

reg = sitk.ImageRegistrationMethod()
reg.SetMetricAsMattesMutualInformation()
reg.SetOptimizerAsRegularStepGradientDescent(learningRate =0.5,
                                           minStep=1e-2,
                                           numberOfIterations=250,
                                           gradientMagnitudeTolerance=1e-4,
                                           estimateLearningRate=reg.Never)
reg.SetOptimizerScales([1, 1, 1, 1.0/250, 1.0/250, 1.0/250])
reg.SetShrinkFactorsPerLevel([1])
reg.SetSmoothingSigmasPerLevel([0])
reg.SetInitialTransform(ctx)
reg.SetInterpolator(sitk.sitkLinear)
reg.SetMetricSamplingPercentage(0.5)
reg.SetMetricSamplingStrategy(reg.RANDOM)
reg.RemoveAllCommands()
reg.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(reg) )
outTx = reg.Execute(sitk.Cast(fixedVolume,sitk.sitkFloat32), sitk.Cast(movingVolume,sitk.sitkFloat32))

sitk.WriteTransform(outTx,outputRegistrationVolume)


# choose optimizer

# choose interpolator

# update

# get output image
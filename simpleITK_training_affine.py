import os
import SimpleITK as sitk
import slicer


def command_iteration(method) :
 print("{0:3} = {1:10.5f} : {2}".format(method.GetOptimizerIteration(),method.GetMetricValue(),method.GetOptimizerPosition()))

def prepare_resampled_values(fixed,moving,transform):
    """ Resample moving into space of fixed, and then rescale the
    intensity values of both to fit in a UInt8 data type.
    """
    f = sitk.ResampleImageFilter()
    f.SetTransform(transform)
    f.SetReferenceImage(fixed)
    print(f.GetTransform)
    deformed_moving = f.Execute(moving)
    fixed = sitk.Cast(sitk.RescaleIntensity(fixed),sitk.sitkUInt8)
    deformed_moving = sitk.Cast(sitk.RescaleIntensity(deformed_moving),sitk.sitkUInt8)
    return fixed, deformed_moving


# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-intraop.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'

fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)

# read input masks

fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'

fixedMask=sitk.ReadImage(fixedMaskFilename, sitk.sitkFloat32)
movingMask=sitk.ReadImage(movingMaskFilename, sitk.sitkFloat32)

# set output file

outputTransform = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutTrans_Affine.h5'
outputVolume = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutVol_Affine.nrrd'


# perform Affine Transformation
# _______________________________

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

# Execute
R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )
affineOut = R.Execute(sitk.Cast(fixedVolume,sitk.sitkFloat32), sitk.Cast(movingVolume,sitk.sitkFloat32))

# Write the transform
sitk.WriteTransform(affineOut,outputTransform)


# perform volume transformation
print ('perform volume transformation')

slicer.mrmlScene.Clear(0)
tfmLogic = slicer.modules.transforms.logic()
vlmLogic=slicer.modules.volumes.logic()


volume=vlmLogic.AddArchetypeVolume('/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd','volume',4)
volumeNode=slicer.mrmlScene.GetNodesByClass('vtkMRMLScalarVolumeNode').GetItemAsObject(0)

transform=tfmLogic.AddTransform(outputTransform,slicer.mrmlScene)

volumeNode.SetAndObserveTransformNodeID(transform.GetID())
tfm = tfmLogic.AddTransform(outputTransform, slicer.mrmlScene)

tfmLogic.hardenTransform(volumeNode)

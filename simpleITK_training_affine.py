import os
import SimpleITK as sitk
import Slicer


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

outputRegistrationTransform = '/Users/peterbehringer/MyTesting/OutTrans_Affine.h5'


# perform Affine Transformation
# _______________________________


tx1 = sitk.AffineTransform(fixedVolume.GetDimension())

ctx = sitk.Transform(tx1) # Set composite transform
                          # This composite transform is update at each stage
                          # and finally will be considered as the output of
                          # the registration process since the InPlace is TRUE.
ctx.SetFixedParameters(ctx.GetFixedParameters()) # hack to force deep copy, as registion is done in place..


# Set the registration filter
R = sitk.ImageRegistrationMethod()
R.SetMetricAsJointHistogramMutualInformation()
R.SetOptimizerAsGradientDescent(learningRate=1.0,
                                numberOfIterations=100,
                                convergenceMinimumValue = 1e-6,
                                convergenceWindowSize = 10,
                                estimateLearningRate = R.EachIteration,
                                maximumStepSizeInPhysicalUnits = 1.0)
R.SetOptimizerScalesFromPhysicalShift()
R.SetShrinkFactorsPerLevel([3,2,1])
R.SetSmoothingSigmasPerLevel([2,1,1])
R.SetInitialTransform(ctx)
R.SetInterpolator(sitk.sitkLinear)
R.SetMetricSamplingPercentage(0.5)
R.SetMetricSamplingStrategy(R.RANDOM)
R.RemoveAllCommands()
R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )
affineOut = R.Execute(fixedVolume,movingVolume)

sitk.WriteTransform(affineOut,outputRegistrationTransform)


"""
def transformVolumes(fiducialsIn, transform, fiducialsOut):
  ### done in Slicer!
  fidLogic = slicer.modules.markups.logic()
  tfmLogic = slicer.modules.transforms.logic()
  fidId = fidLogic.LoadMarkupsFiducials(fiducialsIn, 'na')
  print 'Fiducials loaded:',fidId
  fid = slicer.mrmlScene.GetNodeByID(fidId)
  tfm = tfmLogic.AddTransform(transform, slicer.mrmlScene)
  fid.SetAndObserveTransformNodeID(tfm.GetID())
  tfmLogic.hardenTransform(fid)

  fidStorage = fid.GetStorageNode()
  fidStorage.SetFileName(fiducialsOut)
  fidStorage.WriteData(fid)
  #slicer.mrmlScene.Clear()
"""
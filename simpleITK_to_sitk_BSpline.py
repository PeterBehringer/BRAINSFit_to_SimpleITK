
#=========================================================================
#
#  BRAINSFit to SimpleITK Prostate MR Registration
#
#  BSpline Registration Phase
#
#=========================================================================

import os
import SimpleITK as sitk
import slicer
from math import pi
import numpy as np


def command_iteration(method) :
    if (method.GetOptimizerIteration()==0):
        print("Scales: ", method.GetOptimizerScales())
    print("{0:3} = {1:7.5f} : {2}".format(method.GetOptimizerIteration(),
                                           method.GetMetricValue(),
                                           method.GetOptimizerPosition()))


# mac paths
# =========================================================================

# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-intraop.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
movingImageAffineFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/RegReference_Case1/reg-Affine.nrrd'
fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)
movingVolumeAffine=sitk.ReadImage(movingImageAffineFilename, sitk.sitkFloat32)

# read input masks
fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'
fixedMask=sitk.ReadImage(fixedMaskFilename, sitk.sitkFloat32)
movingMask=sitk.ReadImage(movingMaskFilename, sitk.sitkFloat32)

# set output file paths
eulerTransPath = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/initial_exhaustive.h5'
eulerTransPathAfterRotation = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTransAfterRotation.h5'
rotatedImage='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Rotated_image_1.nrrd'
bestEulerTransPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/bestEulerTrans.h5'
outTxPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTx_exhaustive.h5'
outTxPath1='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTx_exhaustive1.h5'
trans_rigid_beforeRegistrationPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/trans_rigid_beforeRegistration.h5'
outTxAfterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTxAfterRigid.h5'

output_volume_after_init_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_init.nrrd'
fixed_volume_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/fixed_volume.nrrd'
rigid_versor_trans_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_transform_rigid_versor.h5'
output_volume_after_paramPassing_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_paramPassing.nrrd'
rigid_versor_trans_after_rigid_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/rigid_versor_trans_after_rigid.h5'
output_volume_after_rigid_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_rigid.nrrd'
output_volume_after_affine_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_affine.nrrd'
affine_trans_after_affine_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/affine_trans_after_affine.h5'
BSplineOut_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/bspline_trans.h5'
output_volume_after_bspline_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_bspline.nrrd'

FilteredFixedImagePATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/filteredFixedImage.nrrd'
FilteredMovingImagePATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/filteredMovingImage.nrrd'



# BSPLINE REGISTRATION
# =========================================================================
R = sitk.ImageRegistrationMethod()
# Define the second stage as BSpline transform
BSplineOrder=3
# Set fixed parameters of the BSpline transform using transform initializer
transfromDomainMeshSize=[8-BSplineOrder, 8-BSplineOrder, 8-BSplineOrder]
bspline_trans = sitk.BSplineTransformInitializer(fixedVolume, transfromDomainMeshSize, BSplineOrder )
ctx = sitk.Transform()
ctx.AddTransform(bspline_trans)

R.SetInitialTransform(bspline_trans)
R.SetMetricAsJointHistogramMutualInformation()
R.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-7,
                       maximumNumberOfIterations=1500,
                       maximumNumberOfCorrections=25,
                       maximumNumberOfFunctionEvaluations=900,
                       costFunctionConvergenceFactor=1.00E+09)
R.SetShrinkFactorsPerLevel([1])
R.SetSmoothingSigmasPerLevel([0])
R.SetInterpolator(sitk.sitkLinear)
R.SetOptimizerScalesFromPhysicalShift()
R.SetMetricSamplingPercentage(0.5)
R.SetMetricSamplingStrategy(R.RANDOM)
R.SetMetricMovingMask(movingMask)
R.SetMetricFixedMask(fixedMask)

R.RemoveAllCommands()
R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )
out1 = R.Execute(fixedVolume,movingVolumeAffine)

print("-------")
print(out1)
print("Optimizer stop condition: {0}".format(R.GetOptimizerStopConditionDescription()))
print(" Iteration: {0}".format(R.GetOptimizerIteration()))
print(" Metric value: {0}".format(R.GetMetricValue()))


# WRITE TRANSFORM
# =========================================================================

sitk.WriteTransform(out1,BSplineOut_PATH)

# WRITE IMAGES
# =========================================================================

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(out1)

out = resampler.Execute(movingVolumeAffine)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,output_volume_after_bspline_PATH)

simg2 = sitk.Cast(sitk.RescaleIntensity(fixedVolume), sitk.sitkUInt8)
sitk.WriteImage(simg2,fixed_volume_PATH)
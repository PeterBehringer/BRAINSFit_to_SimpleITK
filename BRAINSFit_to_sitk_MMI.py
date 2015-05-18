
#=========================================================================
#
#  BRAINSFit to SimpleITK Prostate MR Registration
#
#=========================================================================

import os
import SimpleITK as sitk
from math import pi
import numpy as np


def command_iteration(method) :
    if (method.GetOptimizerIteration()==0):
        print("Scales: ", method.GetOptimizerScales())
    print("{0:3} = {1:7.5f} : {2}".format(method.GetOptimizerIteration(),
                                           method.GetMetricValue(),
                                           method.GetOptimizerPosition()))


# file paths
# =========================================================================

# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-intraop.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
fixedImage_cropped_Filename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-intraop.nrrd'
movingImage_cropped_Filename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)
movingVolume_cropped=sitk.ReadImage(movingImage_cropped_Filename, sitk.sitkFloat32)
fixedVolume_cropped=sitk.ReadImage(fixedImage_cropped_Filename, sitk.sitkFloat32)


# read input masks
fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'
fixedMask=sitk.ReadImage(fixedMaskFilename, sitk.sitkFloat32)
movingMask=sitk.ReadImage(movingMaskFilename, sitk.sitkFloat32)

# set output file paths
output_volume_after_init_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_init.nrrd'
fixed_volume_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/fixed_volume.nrrd'
rigid_versor_trans_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_transform_rigid_versor.h5'
output_volume_after_paramPassing_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_paramPassing.nrrd'
rigid_versor_trans_after_rigid_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/rigid_versor_trans_after_rigid.h5'
output_volume_after_rigid_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_rigid.nrrd'
output_volume_after_affine_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_after_affine.nrrd'
affine_trans_after_affine_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/affine_trans_after_affine.h5'
euler_trans_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/euler_trans.h5'
resampled_moving_mask_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/resampled_moving_mask.nrrd'
roi_mask_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/roi_mask.nrrd'
bounding_box_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/bounding_box.nrrd'
roi_mask_thresholded_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/roi_mask_thresholded.nrrd'
croppedImage_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/croppedImage.nrrd'
rigid_versor_trans_before_rigid_PATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/rigid_versor_trans_before_rigid.h5'

# INITIALIZATION
# =========================================================================

# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
Reg.SetMetricFixedMask(fixedMask)
Reg.SetMetricMovingMask(movingMask)
Reg.SetMetricSamplingPercentage(.9)
Reg.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.REGULAR)

# Reg.SetMetricAsCorrelation()
Reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

Reg.SetInterpolator(sitk.sitkLinear)

# Set Rotation Parameters for ExhaustiveOptimizer
# =========================================================================

# def 1 degree in rad
one_degree=1.0*pi/180.0

# search in neighbourhood from -12degree to +12degree
angleRange = 12.0* one_degree

# set step size to 3degree
stepSize = 3.0 * one_degree

# samples per axis
sample_per_axis=round((angleRange*2)//stepSize)

# set the number of samples from zero in +- degree direction in each dimension, -> int not double
# [rotationX, rotationY, rotationZ, translationX, translationY, translationZ]
Reg.SetOptimizerAsExhaustive([4,0,4,0,0,0])

# set the step size for each dimension
# [rotationScaleX, rotationScaleY, rotationScaleZ, translationScaleX, translationScaleY, translationScaleZ]
Reg.SetOptimizerScales([stepSize,stepSize,stepSize,1.0,1.0,1.0])

# Initialize the transform with a translation and the center of rotation from the moments of intensity.
# =========================================================================

# Set the Euler3DTransform
euler_trans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))

# Set, Execute & write
Reg.SetInitialTransform(euler_trans,inPlace=True)
Reg.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(Reg))
Reg.Execute(fixedVolume_cropped, movingVolume_cropped)

sitk.WriteTransform(euler_trans,euler_trans_PATH)

# get image volume
resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(euler_trans)

output_volume_after_init = resampler.Execute(movingVolume_cropped)

simg1 = sitk.Cast(sitk.RescaleIntensity(fixedVolume_cropped), sitk.sitkUInt8)
simg2 = sitk.Cast(sitk.RescaleIntensity(output_volume_after_init), sitk.sitkUInt8)

sitk.WriteImage(simg1,fixed_volume_PATH)
sitk.WriteImage(simg2,output_volume_after_init_PATH)

# Enhance searchSpace by Passing Parameters from Euler3DTransform to RigidVersor3DTransform
# =========================================================================

rigid_versor_trans = sitk.VersorRigid3DTransform()
rigid_versor_trans.SetCenter(euler_trans.GetCenter())
rigid_versor_trans.SetTranslation(euler_trans.GetTranslation())
rigid_versor_trans.SetMatrix(euler_trans.GetMatrix())
sitk.WriteTransform(rigid_versor_trans,rigid_versor_trans_before_rigid_PATH)
# make sure params are passed correctly:

print ('euler_trans before parameter passing :')
print ('_____________________________')
print ('')
print euler_trans
print ('')
print ('euler_trans.GetParameters() :')
print ('_____________________________')
print euler_trans.GetParameters()
print ('')
print ('')

print ('rigid_versor_trans after parameter passing :')
print ('_____________________________')
print ('')
print rigid_versor_trans
print ('')
print ('rigid_versor_trans.GetParameters() :')
print ('_____________________________')
print rigid_versor_trans.GetParameters()
print ('')
print ('')

# -> transforms applied to moving image look exactly the same


# CREATE BOUNDING BOX
# =========================================================================
"""
mask_resampler = sitk.ResampleImageFilter()
mask_resampler.SetTransform(euler_trans)
mask_resampler.SetReferenceImage(fixedVolume)
mask_resampler.SetOutputSpacing(fixedVolume.GetSpacing())
mask_resampler.SetOutputDirection(fixedVolume.GetDirection())
mask_resampler.SetOutputOrigin(fixedVolume.GetOrigin())
mask_resampler.SetOutputPixelType(sitk.sitkFloat32)

resampled_moving_mask=mask_resampler.Execute(movingMask)
sitk.WriteImage(resampled_moving_mask,resampled_moving_mask_PATH)

# add transformed mask and fixed mask
add_image_filter=sitk.AddImageFilter()
roi_mask=add_image_filter.Execute(resampled_moving_mask,fixedMask)
sitk.WriteImage(roi_mask,roi_mask_PATH)

# apply threshold
binary_threshold=sitk.BinaryThresholdImageFilter()
binary_threshold.SetLowerThreshold(0.0001)
binary_threshold.SetUpperThreshold(100000)
binary_threshold.SetOutsideValue(0)
binary_threshold.SetInsideValue(1)

roi_mask_thresholded=binary_threshold.Execute(roi_mask)

# get the bounding Box
cast = sitk.CastImageFilter()
cast.SetOutputPixelType(2)
roi_mask_casted = cast.Execute(roi_mask_thresholded)

label_statistics_filter=sitk.LabelStatisticsImageFilter()
label_statistics_filter.Execute(roi_mask_casted,roi_mask_casted)
bounding_box = label_statistics_filter.GetBoundingBox(1)

sitk.WriteImage(roi_mask_thresholded,roi_mask_thresholded_PATH)

size = roi_mask_casted.GetSize()

bbMin = (max(0,bounding_box[0]-30),max(0,bounding_box[2]-30),max(0,bounding_box[4]-5))
bbMax = (size[0]-min(size[0],bounding_box[1]+30),size[1]-min(size[1],bounding_box[3]+30),size[2]-(min(size[2],bounding_box[5]+5)))

lower_boundaries=(bounding_box[0],bounding_box[2],bounding_box[4])
upper_boundaries=(bounding_box[1],bounding_box[3],bounding_box[5])

print ('size : '+str(size))
print ('bbMin :'+str(bbMin))
print ('bbMax :'+str(bbMax))
print ('lower_boundaries :'+str(lower_boundaries))
print ('upper_boundaries :'+str(upper_boundaries))

# CROP MOVING IMAGE
# =========================================================================

crop = sitk.CropImageFilter()
crop.SetLowerBoundaryCropSize(bbMin)
crop.SetUpperBoundaryCropSize(bbMax)
movingVolume_cropped = crop.Execute(movingVolume)

simg1 = sitk.Cast(sitk.RescaleIntensity(movingVolume_cropped), sitk.sitkUInt8)
sitk.WriteImage(simg1,croppedImage_PATH)
"""

# RIGID REGISTRATION
# =========================================================================

# set up registration method
Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(rigid_versor_trans,inPlace=True)

# Reg2.SetMetricAsCorrelation()
Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)
Reg2.SetMetricFixedMask(fixedMask)
Reg2.SetMetricMovingMask(movingMask)
Reg2.SetInterpolator(sitk.sitkLinear)

Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=0.2,
                                          numberOfIterations=1500,
                                          minStep=0.005,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


Reg2.SetOptimizerScales([1.0,1.0,1.0,1.0/1000,1.0/1000,1.0/1000])

smoothingSigmas=[0]
Reg2.SetSmoothingSigmasPerLevel(smoothingSigmas)

Reg2.SetMetricSamplingStrategy(Reg2.RANDOM)
Reg2.SetMetricSamplingPercentage(0.002)

Reg2.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(True)

shrinkFactors=[1]
Reg2.SetShrinkFactorsPerLevel(shrinkFactors)

# execute
print ('Now lets try Rigid')
print ('')
Reg2.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(Reg2))
Reg2.Execute(fixedVolume_cropped, movingVolume_cropped)

sitk.WriteTransform(rigid_versor_trans,rigid_versor_trans_after_rigid_PATH)

# resample image
resampler.SetReferenceImage(fixedVolume_cropped)
resampler.SetTransform(rigid_versor_trans)
out = resampler.Execute(movingVolume_cropped)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,output_volume_after_rigid_PATH)
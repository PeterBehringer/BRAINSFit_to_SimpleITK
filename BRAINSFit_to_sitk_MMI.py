
#=========================================================================
#
#  BRAINSFit to SimpleITK Prostate MR Registration
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


# file paths
# =========================================================================

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


# INITIALIZATION
# =========================================================================

# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
# Reg.SetMetricFixedMask(fixedMask)
# Reg.SetMetricMovingMask(movingMask)

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
Reg.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(euler_trans,euler_trans_PATH)

# get image volume
resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(euler_trans)

output_volume_after_init = resampler.Execute(movingVolume)

simg1 = sitk.Cast(sitk.RescaleIntensity(fixedVolume), sitk.sitkUInt8)
simg2 = sitk.Cast(sitk.RescaleIntensity(output_volume_after_init), sitk.sitkUInt8)

sitk.WriteImage(simg1,fixed_volume_PATH)
sitk.WriteImage(simg2,output_volume_after_init_PATH)

# Enhance searchSpace by Passing Parameters from Euler3DTransform to RigidVersor3DTransform
# =========================================================================

rigid_versor_trans = sitk.VersorRigid3DTransform()
rigid_versor_trans.SetCenter(euler_trans.GetCenter())
rigid_versor_trans.SetTranslation(euler_trans.GetTranslation())
rigid_versor_trans.SetMatrix(euler_trans.GetMatrix())

# make sure params are passed correctly:

print ('euler_trans before parameter passing :')
print ('_____________________________')
print ()
print euler_trans
print ()
print ('euler_trans.GetParameters() :')
print ('_____________________________')
print euler_trans.GetParameters()
print ()
print ()

print ('rigid_versor_trans after parameter passing :')
print ('_____________________________')
print ()
print rigid_versor_trans
print ()
print ('rigid_versor_trans.GetParameters() :')
print ('_____________________________')
print rigid_versor_trans.GetParameters()
print ()
print ()

# -> checked/transforms applied to moving image look exactly the same


# RIGID REGISTRATION
# =========================================================================


# set up registration method
Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(rigid_versor_trans,inPlace=True)

# Reg2.SetMetricAsCorrelation()
Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

# Reg2.SetMetricFixedMask(fixedMask)
# Reg2.SetMetricMovingMask(movingMask)
Reg2.SetInterpolator(sitk.sitkLinear)

Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=0.2,
                                          numberOfIterations=1500,
                                          minStep=0.005,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


Reg2.SetOptimizerScales([1.0,1.0,1.0,1.0/1000,1.0/1000,1.0/1000])
# Reg2.SetOptimizerScalesFromJacobian()

smoothingSigmas=[0]
Reg2.SetSmoothingSigmasPerLevel(smoothingSigmas)
Reg2.SetMetricSamplingStrategy(Reg2.RANDOM)
Reg2.SetMetricSamplingPercentage(1)
Reg2.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(True)

shrinkFactors=[1]
Reg2.SetShrinkFactorsPerLevel(shrinkFactors)



# execute
print ('Now lets try Rigid')
Reg2.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(Reg2))
Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(rigid_versor_trans,rigid_versor_trans_after_rigid_PATH)

# resample image
resampler.SetTransform(rigid_versor_trans)
out = resampler.Execute(movingVolume)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,output_volume_after_rigid_PATH)

print ('rigid_versor_trans after rigid optimization :')
print ('_____________________________')
print ()
print rigid_versor_trans
print ()
print ('rigid_versor_trans.GetParameters() :')
print ('_____________________________')
print rigid_versor_trans.GetParameters()
print ()
print ()
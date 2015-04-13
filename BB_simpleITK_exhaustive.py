
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

"""

# windows paths
# =========================================================================

# read input volumes
fixedImageFilename = 'C:\Users\pb691\Desktop\MyImageData\Images\Case1-t2ax-intraop.nrrd'
movingImageFilename= 'C:\Users\pb691\Desktop\MyImageData\Images\Case1-t2ax-N4.nrrd'
fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)

# read input masks
fixedMaskFilename = 'C:\Users\pb691\Desktop\MyImageData\Segmentations\Rater1\Case1-t2ax-intraop-TG-rater1.nrrd'
movingMaskFilename= 'C:\Users\pb691\Desktop\MyImageData\Segmentations\Rater1\Case1-t2ax-TG-rater1.nrrd'
fixedMask=sitk.ReadImage(fixedMaskFilename, sitk.sitkFloat32)
movingMask=sitk.ReadImage(movingMaskFilename, sitk.sitkFloat32)

# set output file paths
outputTransformInit='C:\Users\pb691\Desktop\MyTesting\output_transform_afterInit.h5'
outTxAfterRigidPath='C:\Users\pb691\Desktop\MyTesting\output_transform_afterRigid.h5'
fixedVolumePath='C:\Users\pb691\Desktop\MyTesting\FixedVolume.nrrd'
outVolumePath='C:\Users\pb691\Desktop\MyTesting\output_volume_afterInit.nrrd'
afterRigidPath='C:\Users\pb691\Desktop\MyTesting\output_volume_afterRigid.nrrd'
tx2outputPath='C:\Users\pb691\Desktop\MyTesting\output_transform_afterParamsPassing.h5'
rigid_versorPath='C:\Users\pb691\Desktop\MyTesting\Test__rigid_versor.h5'
"""

# mac paths
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

FilteredFixedImagePATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/filteredFixedImage.nrrd'
FilteredMovingImagePATH='/Users/peterbehringer/MyTesting/SimpleITK_Tests/filteredMovingImage.nrrd'

"""
# NORMALIZATION OF INPUT IMAGES
# tried to fix Mutual Information by Normalization and Histogram Matching

ImFilter=sitk.NormalizeImageFilter()
HistMatchFilter=sitk.HistogramMatchingImageFilter()

fImg1=ImFilter.Execute(fixedVolume)
fImg2=ImFilter.Execute(movingVolume)

FilteredFixedImage=HistMatchFilter.Execute(fImg1, fImg1)
FilteredMovingImage=HistMatchFilter.Execute(fImg2, fImg1)


sitk.WriteImage(FilteredFixedImage,FilteredFixedImagePATH)
sitk.WriteImage(FilteredMovingImage,FilteredMovingImagePATH)
"""

# INITIALIZATION
# =========================================================================

# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
Reg.SetMetricFixedMask(fixedMask)
Reg.SetMetricMovingMask(movingMask)

Reg.SetMetricAsCorrelation()
# Reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

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
eulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))

# Set, Execute & write
Reg.SetInitialTransform(eulerTrans,inPlace=True)
Reg.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(Reg))
Reg.Execute(fixedVolume, movingVolume)

print eulerTrans
sitk.WriteTransform(eulerTrans,outTxPath1)

# get image volume
resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(eulerTrans)

output_volume_after_init = resampler.Execute(movingVolume)

simg1 = sitk.Cast(sitk.RescaleIntensity(fixedVolume), sitk.sitkUInt8)
simg2 = sitk.Cast(sitk.RescaleIntensity(output_volume_after_init), sitk.sitkUInt8)

sitk.WriteImage(simg1,fixed_volume_PATH)
sitk.WriteImage(simg2,output_volume_after_init_PATH)

# Enhance searchSpace by Passing Parameters from Euler3DTransform to RigidVersor3DTransform
# =========================================================================

rigid_versor_trans = sitk.VersorRigid3DTransform()
rigid_versor_trans.SetCenter(eulerTrans.GetCenter())
rigid_versor_trans.SetTranslation(eulerTrans.GetTranslation())
rigid_versor_trans.SetMatrix(eulerTrans.GetMatrix())

# make sure params are passed correctly:
print ('Transform after parameter passing :')
print ('_____________________________')
print ()
print rigid_versor_trans

print ('eulerTrans')
print ('______________')
print (eulerTrans)

# -> checked, transforms applied to moving image look exactly the same


# RIGID REGISTRATION
# =========================================================================


# set up registration method
Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(rigid_versor_trans,inPlace=True)

Reg2.SetMetricAsCorrelation()
# Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

Reg2.SetMetricFixedMask(fixedMask)
Reg2.SetMetricMovingMask(movingMask)
Reg2.SetInterpolator(sitk.sitkLinear)

# BRAINSFIT IGT SLICER 3.6 PARAMS
# --minimumStepLength	0.005
# --numberOfIterations	1500
# --translationScale 1000

# BRAINSFitHelperTemplate.hxx PARAMS
# m_MaximumStepLength(0.2)
# m_MinimumStepLength(1, 0.005)
# m_RelaxationFactor(0.5)
# m_ProjectedGradientTolerance(1e-5)

# PARAMETER SetOptimizerAsRegularStepGradientDescent
# double learningRate,                                    m_MaximumStepLength = 0.2
# double minStep,                                         0.005
# unsigned int numberOfIterations,                        1500
# double relaxationFactor                                 0.5
# double gradientMagnitudeTolerance                       1e-5
# EstimateLearningRateType estimateLearningRate=Never     ?
# double maximumStepSizeInPhysicalUnits                   0.2

# learningRate -> m_MaximumStepLength -> 0.2000

# parameter, where I dont know where to set yet:
# samplingPercentage = 0.002
# reproportionScale = 1
# skewScale = 1




Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=0.2,
                                          numberOfIterations=1500,
                                          minStep=0.005,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)



# set Optimizer Scales
# as done in BRAINSCommonLib/genericRegistrationHelper.hxx #L255 ff for VersorRigid3D
#
# optimizerScales[0] : 1.0
# optimizerScales[1] : 1.0
# optimizerScales[2] : 1.0
# optimizerScales[3] : 1.0 / m_translationScale
# optimizerScales[4] : 1.0 / m_translationScale
# optimizerScales[5] : 1.0 / m_translationScale

Reg2.SetOptimizerScales([1.0,1.0,1.0,1.0/1000,1.0/1000,1.0/1000])
Reg2.SetOptimizerScalesFromJacobian()


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



"""


# AFFINE REGISTRATION

# --useScaleVersor3D

# --useScaleSkewVersor3D

# =========================================================================


# Enhance searchSpace by Passing Parameters from RigidVersor3DTransform to AffineTransform()
# =========================================================================

affine_trans = sitk.AffineTransform(fixedVolume.GetDimension())
affine_trans.SetCenter(eulerTrans.GetCenter())
affine_trans.SetTranslation(eulerTrans.GetTranslation())
affine_trans.SetMatrix(eulerTrans.GetMatrix())

# make sure params are passed correctly:
print ('Affine Transform after parameter passing :')
print ('_____________________________')
print ()
print affine_trans

print ('Reference: eulerTrans')
print ('______________')
print (eulerTrans)


# set up registration method
Reg3=sitk.ImageRegistrationMethod()
Reg3.SetInitialTransform(affine_trans,inPlace=True)

Reg3.SetMetricAsCorrelation()
# Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

Reg3.SetMetricFixedMask(fixedMask)
Reg3.SetMetricMovingMask(movingMask)
Reg3.SetInterpolator(sitk.sitkLinear)

# BRAINSFIT IGT SLICER 3.6 PARAMS
# --minimumStepLength	0.005
# --numberOfIterations	1500
# --translationScale 1000

# BRAINSFitHelperTemplate.hxx PARAMS
# m_MaximumStepLength(0.2)
# m_MinimumStepLength(1, 0.005)
# m_RelaxationFactor(0.5)
# m_ProjectedGradientTolerance(1e-5)

# PARAMETER SetOptimizerAsRegularStepGradientDescent
# double learningRate,                                    ?
# double minStep,                                         0.005
# unsigned int numberOfIterations,                        1500
# double relaxationFactor                                 0.5
# double gradientMagnitudeTolerance                       1e-5
# EstimateLearningRateType estimateLearningRate=Never     ?
# double maximumStepSizeInPhysicalUnits                   0.2

Reg3.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0,
                                          minStep=0.005,
                                          numberOfIterations=1500,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)
# set Optimizer Scales
# as done in BRAINSCommonLib/genericRegistrationHelper.hxx #L255 ff for VersorRigid3D
#
# optimizerScales[0] : 1.0
# optimizerScales[1] : 1.0
# optimizerScales[2] : 1.0
# optimizerScales[3] : 1.0 / m_translationScale
# optimizerScales[4] : 1.0 / m_translationScale
# optimizerScales[5] : 1.0 / m_translationScale

Reg3.SetOptimizerScales([1.0,1.0,1.0,1.0/1000,1.0/1000,1.0/1000])


# execute
print ('Now lets try Affine')
Reg3.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(Reg3) )
Reg3.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(affine_trans,affine_trans_after_affine_PATH)

# resample image
resampler.SetTransform(affine_trans)
out = resampler.Execute(movingVolume)
simg3 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg3,output_volume_after_affine_PATH)



# BSPLINE REGISTRATION
# =========================================================================
"""
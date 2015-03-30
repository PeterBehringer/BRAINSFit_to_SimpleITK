from __future__ import print_function
from __future__ import division


import SimpleITK as sitk
import sys
import os
from math import pi


def command_iteration(method) :
    if (method.GetOptimizerIteration()==0):
        print("Scales: ", method.GetOptimizerScales())
    print("{0:3} = {1:7.5f} : {2}".format(method.GetOptimizerIteration(),
                                           method.GetMetricValue(),
                                           method.GetOptimizerPosition()))



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


# Initialization
# _______________________________
R = sitk.ImageRegistrationMethod()

R.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

sample_per_axis=12
if fixedVolume.GetDimension() == 2:
    tx = sitk.Euler2DTransform()
    # Set the number of samples (radius) in each dimension, with a
    # default step size of 1.0
    R.SetOptimizerAsExhaustive([sample_per_axis//2,0,0])
    # Utilize the scale to set the step size for each dimension
    R.SetOptimizerScales([2.0*pi/sample_per_axis, 1.0,1.0])
elif fixedVolume.GetDimension() == 3:
    tx = sitk.Euler3DTransform()
    R.SetOptimizerAsExhaustive([sample_per_axis//2,sample_per_axis//2,sample_per_axis//4,0,0,0])
    R.SetOptimizerScales([2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,1.0,1.0,1.0])

# Initialize the transform with a translation and the center of
# rotation from the moments of intensity.
tx = sitk.CenteredTransformInitializer(fixedMask, movingMask, tx)

R.SetInitialTransform(tx)

R.SetInterpolator(sitk.sitkLinear)

R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )

outTx = R.Execute(fixedVolume, movingVolume)

print("-------")
print(outTx)
print("Optimizer stop condition: {0}".format(R.GetOptimizerStopConditionDescription()))
print(" Iteration: {0}".format(R.GetOptimizerIteration()))
print(" Metric value: {0}".format(R.GetMetricValue()))

sitk.WriteTransform(outTx,  outputTransformInit)

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(outTx)

out = resampler.Execute(movingVolume)

simg1 = sitk.Cast(sitk.RescaleIntensity(fixedVolume), sitk.sitkUInt8)
simg2 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)

sitk.WriteImage(fixedVolume,fixedVolumePath)
sitk.WriteImage(out,outVolumePath)


# RIGID REGISTRATION PHASE
# _______________________________

tx2=sitk.VersorRigid3DTransform()
outTx.AddTransform(tx2)

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(outTx)
Reg2.SetMetricAsCorrelation()
Reg2.SetMetricFixedMask(fixedMask)
Reg2.SetMetricMovingMask(movingMask)
Reg2.SetInterpolator(sitk.sitkLinear)


# BRAINSFIT IGT SLICER 3.6 PARAMETER
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
# unsigned int numberOfIterations,                        100 (1500 actually)
# double relaxationFactor                                 0.5
# double gradientMagnitudeTolerance                       1e-5
# EstimateLearningRateType estimateLearningRate=Never     ?
# double maximumStepSizeInPhysicalUnits                   0.2

Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0,
                                          minStep=0.005,
                                          numberOfIterations=1500,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


outTxAfterRigid = Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(outTxAfterRigid,outTxAfterRigidPath)

# resample image
resampler.SetTransform(outTxAfterRigid)
out = resampler.Execute(movingVolume)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,afterRigidPath)
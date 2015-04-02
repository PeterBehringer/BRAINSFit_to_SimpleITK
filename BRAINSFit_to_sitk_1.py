
#=========================================================================
#
#  BRAINSFit to SimpleITK Prostate MR Registration
#
#=========================================================================

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


"""
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
trans_rigid_beforeRegistrationPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/trans_rigid_beforeRegistration.h5'
outTxAfterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTxAfterRigid.h5'

# set output file paths
outputTransformInit='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_transform_afterInit.h5'
outTxAfterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_transform_afterRigid.h5'
fixedVolumePath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/FixedVolume.nrrd'
outVolumePath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_afterInit.nrrd'
afterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_volume_afterRigid.nrrd'
tx2outputPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/output_transform_afterParamsPassing.h5'
rigid_versorPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Test__rigid_versor.h5'

######  REGISTRATION
######  Initialization and Rotation

Reg = sitk.ImageRegistrationMethod()

Reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

# def 1 degree in rad
one_degree=1.0*pi/180.0

# search in neighbourhood from -12degree to +12degree
angleRange = 12.0* one_degree

# set step size to 3degree
stepSize = 3.0 * one_degree

# samples per axis
sample_per_axis=round((angleRange*2)//stepSize)

print ('angle range : '+str(angleRange))
print ('sample per axis : '+str(sample_per_axis))
print ('step size : '+str(stepSize))

tx = sitk.Euler3DTransform()

sample_per_axis = 8

# set the number of samples (radius) in each dimension
Reg.SetOptimizerAsExhaustive([sample_per_axis/2,0,sample_per_axis/2,0,0,0])

# set the step size for each dimension
Reg.SetOptimizerScales([0.066*2*pi/sample_per_axis,0.0,0.066*2*pi/sample_per_axis,1.0,1.0,1.0])

# Initialize the transform with a translation and the center of
# rotation from the moments of intensity.
tx = sitk.CenteredTransformInitializer(fixedMask, movingMask, tx)

Reg.SetInitialTransform(tx,inPlace=True)

Reg.SetInterpolator(sitk.sitkLinear)

Reg.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )

Reg.Execute(fixedVolume, movingVolume)

print("-------")
print(tx)
print("Optimizer stop condition: {0}".format(R.GetOptimizerStopConditionDescription()))
print(" Iteration: {0}".format(R.GetOptimizerIteration()))
print(" Metric value: {0}".format(R.GetMetricValue()))


sitk.WriteTransform(tx,outputTransformInit)

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(bestEulerTrans)

out = resampler.Execute(movingVolume)

simg1 = sitk.Cast(sitk.RescaleIntensity(fixedVolume), sitk.sitkUInt8)
simg2 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)

sitk.WriteImage(fixedVolume,fixedVolumePath)
sitk.WriteImage(out,outVolumePath)


# RIGID REGISTRATION PHASE
# _______________________________

print ('tx Initial after Initialization and Euler Rotation :')
print ('_____________________________')
print bestEulerTrans
print ()
print ('type of bestEulerTrans :')
print type(bestEulerTrans)
print ()

# (1) pass parameters
"""
rigid_versor = sitk.VersorRigid3DTransform()

rotationCenter = bestEulerTrans.GetCenter()
rotation = sitk.VersorTransform(bestEulerTrans.GetMatrix(), rotationCenter)

rigid_versor.SetCenter(rotation.GetCenter())
rigid_versor.SetRotation(rotation.GetVersor())
rigid_versor.SetMatrix(rotation.GetMatrix())


print ('rigid_versor')
print ('______________')
print (rigid_versor)

sitk.WriteTransform(rigid_versor,  rigid_versorPath)

# rigid_euler = sitk.Euler3DTransform()
# rigid_euler.SetMatrix(rotation.GetMatrix())
# rigid_euler.SetCenter(rotation.GetCenter())




tx2=sitk.VersorRigid3DTransform()
tx2.SetParameters(outTx.GetParameters())

tx2=sitk.VersorRigid3DTransform()
outTx.AddTransform(tx2)

sitk.WriteTransform(tx2,  tx2outputPath)

"""

# (2) set up registration method

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(bestEulerTrans)


Reg2.SetMetricAsCorrelation()
#Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

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
# set Optimizer Scales
# as done in BRAINSCommonLib/genericRegistrationHelper.hxx L 255 ff for VersorRigid3D
#
# optimizerScales[0] : 1.0
# optimizerScales[1] : 1.0
# optimizerScales[2] : 1.0
# optimizerScales[3] : 1.0 / m_translationScale
# optimizerScales[4] : 1.0 / m_translationScale
# optimizerScales[5] : 1.0 / m_translationScale

Reg2.SetOptimizerScales([1.0,1.0,1.0,1.0/1000,1.0/1000,1.0/1000])

print ('Now lets try Rigid')
Reg2.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(Reg2) )

# (3) execute

outTxAfterRigid = Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(outTxAfterRigid,outTxAfterRigidPath)

# (4) resample image

resampler.SetTransform(outTxAfterRigid)
out = resampler.Execute(movingVolume)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,afterRigidPath)


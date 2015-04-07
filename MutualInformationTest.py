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

# read input volumes

fixedImageFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
movingImageFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd'
fixedVolume=sitk.ReadImage(fixedImageFilename, sitk.sitkFloat32)
movingVolume=sitk.ReadImage(movingImageFilename, sitk.sitkFloat32)

# read input masks
fixedMaskFilename = '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'
movingMaskFilename= '/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Segmentations/Rater1/Case1-t2ax-TG-rater1.nrrd'
fixedMask=sitk.ReadImage(fixedMaskFilename, sitk.sitkFloat32)
movingMask=sitk.ReadImage(movingMaskFilename, sitk.sitkFloat32)

# set output file paths
outTxAfterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTxAfterRigid.h5'
afterRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outVolAfterRigid.nrrd'


# RIGID REGISTRATION PHASE
# _______________________________



rigid_versor=sitk.VersorRigid3DTransform()

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(rigid_versor,inPlace=True)

#Reg2.SetMetricAsCorrelation()
Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

Reg2.SetMetricFixedMask(fixedMask)
Reg2.SetMetricMovingMask(movingMask)
Reg2.SetInterpolator(sitk.sitkLinear)
Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=0.2000,
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

Reg2.Execute(fixedVolume, movingVolume)
sitk.WriteTransform(rigid_versor,outTxAfterRigidPath)

# (4) resample image

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(outTxAfterRigid)

out = resampler.Execute(movingVolume)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,afterRigidPath)



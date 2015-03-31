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



# INITIALIZATION (Translation & Rotation)
# _______________________________

# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
Reg.SetMetricFixedMask(fixedMask)
Reg.SetMetricMovingMask(movingMask)
Reg.SetMetricAsCorrelation()
Reg.SetInterpolator(sitk.sitkLinear)

# use exhaustive Optimizer
sample_per_axis=12
Reg.SetOptimizerAsExhaustive([sample_per_axis/2,sample_per_axis/2,sample_per_axis/4,0,0,0])
Reg.SetOptimizerScales([2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,1.0,1.0,1.0])

# Set the Euler3DTransform
eulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))
Reg.SetInitialTransform(eulerTrans)

# write the transform
sitk.WriteTransform(eulerTrans,eulerTransPath)

# Execute & write
Reg.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(Reg))
outTx = Reg.Execute(fixedVolume, movingVolume)
sitk.WriteTransform(outTx,outTxPath)

# RIGID REGISTRATION PHASE
# _______________________________

# trans_afterInit=sitk.ReadTransform('/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTx_exhaustive.h5')

quickSetVersor=sitk.VersorRigid3DTransform()
quickSetVersor.SetCenter(outTx.GetCenter())
quickSetVersor.SetTranslation(outTx.GetTranslation())
quickSetVersor.SetMatrix(outTx.GetMatrix())

sitk.WriteTransform(quickSetVersor,trans_rigid_beforeRegistrationPath)

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(quickSetVersor)
Reg2.SetMetricAsCorrelation()
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

# PARAMS SetOptimizerAsRegularStepGradientDescent
# double learningRate,                                    ?
# double minStep,                                         0.005
# unsigned int numberOfIterations,                        100 (1500 actually)
# double relaxationFactor                                 0.5
# double gradientMagnitudeTolerance                       1e-5
# EstimateLearningRateType estimateLearningRate=Never     ?
# double maximumStepSizeInPhysicalUnits                   0.2

Reg2.SetOptimizerAsRegularStepGradientDescent(minStep=0.005,
                                          numberOfIterations=100,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


outTxAfterRigid = Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(outTxAfterRigid,outTxAfterRigidPath)


# AFFINE REGISTRATION PHASE
# _______________________________




# BSPLINE REGISTRATION PHASE
# _______________________________

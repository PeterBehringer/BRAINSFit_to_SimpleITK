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

outputTransform = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutTrans_Rigid_1.h5'
ctxInitialPath ='/Users/peterbehringer/MyTesting/SimpleITK_Tests/PetersTry_ctxInitial.h5'
outTxPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Brad_afterInitAndRotation.h5'
txBeforeRigidPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Brad_beforeRigidWithInitialParams.h5'
outTxAfterAffinePath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Brad_outTxAfterRigidPath.h5'

# INITIALIZATION AND EULER ROTATION
# _______________________________

R = sitk.ImageRegistrationMethod()

R.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

sample_per_axis=12
tx = sitk.Euler3DTransform()
R.SetOptimizerAsExhaustive([sample_per_axis/2,sample_per_axis/2,sample_per_axis/4,0,0,0])
R.SetOptimizerScales([2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,2.0*pi/sample_per_axis,1.0,1.0,1.0])

# Initialize the transform with a translation and the center of
# rotation from the moments of intensity.

tx = sitk.CenteredTransformInitializer(fixedMask, movingMask, tx)

R.SetInitialTransform(tx)
R.SetInterpolator(sitk.sitkLinear)
R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )

outTx = R.Execute(fixedVolume, movingVolume)


print ('tx Initial after Initialization and Euler Rotation :')
print ('_____________________________')
print ()
print tx
print ()
sitk.WriteTransform(outTx,outTxPath)

# RIGID REGISTRATION PHASE
# _______________________________

tx2=sitk.VersorRigid3DTransform()
tx.AddTransform(tx2)

# make sure params are passed correctly:
print ('txBeforeRigidPath :')
print ('_____________________________')
print ()
print tx
print ()
sitk.WriteTransform(tx,txBeforeRigidPath)

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(tx)
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

Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0,
                                          minStep=0.005,
                                          numberOfIterations=100,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


outTxAfterAffine = Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(outTxAfterAffine,outTxAfterAffinePath)

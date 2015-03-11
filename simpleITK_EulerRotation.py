import os
import SimpleITK as sitk
import slicer
from math import pi
import numpy as np

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
outputVolume = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutVol_Rigid_1.nrrd'
outputTransform_Initializer = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutTrans_Initializer_Rigid__1.h5'
ctx1Data = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/ctx1.h5'
ctx2Data = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/ctx2.h5'
eulerTransPath = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTrans.h5'
eulerTransPathAfterRotation = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTransAfterRotation.h5'
rotatedImage='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Rotated_image_1.nrrd'
bestEulerTransPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/bestEulerTrans.h5'

# INITIALIZATION - ROTATION - METRIC MEASURE
# _______________________________

# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
Reg.SetMetricFixedMask(fixedMask)
Reg.SetMetricMovingMask(movingMask)
Reg.SetMetricAsCorrelation()
Reg.SetInterpolator(sitk.sitkLinear)
Reg.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0,
                                          minStep=1e-4,
                                          numberOfIterations=1,
                                          gradientMagnitudeTolerance=1e-8 )


print ('>> Step 0: default >> Metricvalue : '+str(Reg.MetricEvaluate(fixedVolume,movingVolume)))

# Set the Euler3DTransform
eulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))
Reg.SetInitialTransform(eulerTrans)

# write the transform
sitk.WriteTransform(eulerTrans,eulerTransPath)

print ('>> Step 1: after centered Initialization >> Metricvalue : '+str(Reg.MetricEvaluate(fixedVolume,movingVolume)))

# print out center
# center=eulerTrans.GetCenter()
# print ('Rotation Center : '+str(center))

# set Rotation
eulerTrans.SetRotation(0,0,45)
# print ('Parameters eulerTrans : '+str(eulerTrans.GetParameters()))

# Update Transform
Reg.SetInitialTransform(eulerTrans)

# write the transform
sitk.WriteTransform(eulerTrans,eulerTransPathAfterRotation)

print ('>> Step 2: after 45 Rotation axial >> Metricvalue : '+str(Reg.MetricEvaluate(fixedVolume,movingVolume)))

# set Rotation
eulerTrans.SetRotation(23,0,0)
# print ('Parameters eulerTrans : '+str(eulerTrans.GetParameters()))

# Update Transform
Reg.SetInitialTransform(eulerTrans)

print ('>> Step 3: after 23 Rotation axial >> Metricvalue : '+str(Reg.MetricEvaluate(fixedVolume,movingVolume)))

one_degree=1.0*pi/180.0

axAngleRange = 2.0
sagAngleRange = 2.0

axStepSize = 3.0 * one_degree
sagStepSize = 3.0 * one_degree

# set current Metric value
initialMetricValue= Reg.MetricEvaluate(fixedVolume,movingVolume)

# initialize output Transform with Translation from CenteredTransformInitializer
bestEulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))


for axAngle in np.arange(0.0,axAngleRange,axStepSize):


        eulerTrans.SetRotation(0,0,axAngle)
        Reg.SetInitialTransform(eulerTrans)
        currentMetricValue = Reg.MetricEvaluate(fixedVolume,movingVolume)

        if abs(currentMetricValue) < abs(initialMetricValue):
            bestEulerTrans.SetRotation(0,0,axAngle)
            # print ('new best Euler Trans found at axAngle = '+str(axAngle))

sitk.WriteTransform(bestEulerTrans,bestEulerTransPath)

"""
for axAngle in np.arange(0.0,axAngleRange,axStepSize):
    for sagAngle in np.arange(0.0,sagAngleRange,sagStepSize):


        eulerTrans.SetRotation(sagAngle,0,axAngle)
        Reg.SetInitialTransform(eulerTrans)
        currentMetricValue = Reg.MetricEvaluate(fixedVolume,movingVolume)

        if abs(currentMetricValue) < abs(initialMetricValue):
            bestEulerTrans.SetRotation(sagAngle,0,axAngle)
            print ('new best Euler Trans found at sagAngle = '+str(sagAngle)+' and axAngle = '+str(axAngle))

sitk.WriteTransform(bestEulerTrans,bestEulerTransPath)
"""
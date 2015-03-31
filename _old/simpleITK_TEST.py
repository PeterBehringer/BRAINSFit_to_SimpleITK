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
eulerTransPath = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTrans_TEST.h5'
eulerTransPathAfterRotation = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTransAfterRotation.h5'
rotatedImage='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Rotated_image_1.nrrd'
bestEulerTransPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/bestEulerTrans.h5'
outTxPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/outTx.h5'
quickSetVersorPath = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/quickSetVersor_TEST.h5'

# INITIALIZATION
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

# Set the Initial Euler3DTransform
# ___________________________

eulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))
Reg.SetInitialTransform(eulerTrans)

# Set Rotation
eulerTrans.SetRotation(0.0,0.0,45.0)

# Print out Parameters
print ('')
print ('')
print ('PARAMETER eulerTrans')
print ('__________________________')
print ('Center = '+str(eulerTrans.GetCenter()))
print ('Translation = ' +str(eulerTrans.GetTranslation()))
print ('Parameters = '+str(eulerTrans.GetParameters()))
print ('Matrix = '+str(eulerTrans.GetMatrix()))

# write the transform
sitk.WriteTransform(eulerTrans,eulerTransPath)

# Set the VersorRigid3DTransform
# ___________________________

# set parameter
quickSetVersor=sitk.VersorRigid3DTransform()
quickSetVersor.SetCenter(eulerTrans.GetCenter())
quickSetVersor.SetTranslation(eulerTrans.GetTranslation())
quickSetVersor.SetMatrix(eulerTrans.GetMatrix())

# Print out Parameters
print ('')
print ('')
print ('PARAMETER quickSetVersor')
print ('__________________________')
print ('Center = '+str(quickSetVersor.GetCenter()))
print ('Translation = ' +str(quickSetVersor.GetTranslation()))
print ('Matrix = '+str(quickSetVersor.GetMatrix()))
print ('Parameters = '+str(quickSetVersor.GetParameters()))
print ('Versor = '+str(quickSetVersor.GetVersor()))

# write the transform
sitk.WriteTransform(quickSetVersor,quickSetVersorPath)

###
# -> both transforms look the same

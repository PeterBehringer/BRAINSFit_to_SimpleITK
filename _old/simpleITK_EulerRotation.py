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
eulerTransPath = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/PetersTry_InitialTrans.h5'
eulerTransPathAfterRotation = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/eulerTransAfterRotation.h5'
rotatedImage='/Users/peterbehringer/MyTesting/SimpleITK_Tests/Rotated_image_1.nrrd'
bestEulerTransPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/PetersTry_AfterEulerRotation.h5'
outTxPath='/Users/peterbehringer/MyTesting/SimpleITK_Tests/PetersTry_AfterRigid.h5'
quickSetVersorInitial='/Users/peterbehringer/MyTesting/SimpleITK_Tests/quickVersorInitial.h5'


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

# Set the Euler3DTransform
eulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))
Reg.SetInitialTransform(eulerTrans)

# write the transform
sitk.WriteTransform(eulerTrans,eulerTransPath)


# ROTATE & MEASURE METRIC
# like here https://github.com/BRAINSia/BRAINSTools/blob/19fa37dfbdee37deff4ccee412bb601f7a787bda/BRAINSCommonLib/BRAINSFitHelperTemplate.hxx#L325-L339
# _______________________________

one_degree=1.0*pi/180.0

axAngleRange = 1.0
# sagAngleRange = 12.0

axStepSize = 3.0 * one_degree
# sagStepSize = 3.0 * one_degree

# set current Metric value
initialMetricValue= Reg.MetricEvaluate(fixedVolume,movingVolume)

# initialize output Transform with Translation from CenteredTransformInitializer
bestEulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))

for axAngle in np.arange(-axAngleRange,axAngleRange,axStepSize):
    #for sagAngle in np.arange(-sagAngleRange,sagAngleRange,sagStepSize):

        eulerTrans.SetRotation(0,0,axAngle)
        print ('current axAngle : '+str(axAngle))
        Reg.SetInitialTransform(eulerTrans)
        currentMetricValue = Reg.MetricEvaluate(fixedVolume,movingVolume)

        if currentMetricValue < initialMetricValue:
            print ('new best axAngle : '+str(axAngle))
            bestEulerTrans.SetRotation(0,0,axAngle)
            # print ('new best Euler Trans found at sagAngle = '+str(sagAngle)+' and axAngle = '+str(axAngle))

sitk.WriteTransform(bestEulerTrans,bestEulerTransPath)

# RIGID REGISTRATION PHASE
# _______________________________

quickSetVersor=sitk.VersorRigid3DTransform()



# quickSetVersor.SetCenter(bestEulerTrans.GetCenter())
# quickSetVersor.SetTranslation(bestEulerTrans.GetTranslation())
# quickSetVersor.SetMatrix(bestEulerTrans.GetMatrix())

sitk.WriteTransform(quickSetVersor,quickSetVersorInitial)

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

Reg2.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0,
                                          minStep=0.005,
                                          numberOfIterations=100,
                                          relaxationFactor = 0.5,
                                          gradientMagnitudeTolerance=1e-5,
                                          maximumStepSizeInPhysicalUnits=0.2)


outTx = Reg2.Execute(fixedVolume, movingVolume)

sitk.WriteTransform(outTx,outTxPath)

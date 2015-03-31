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

outputTransform = 'C:\Users\pb691\Desktop\MyTesting\outputTransform.h5'
outputVolume = 'C:\Users\pb691\Desktop\MyTesting\outputVolume.h5'
outputTransform_Initializer = 'C:\Users\pb691\Desktop\MyTesting\outputTransform_Initializer.h5'
ctx1Data = 'C:\Users\pb691\Desktop\MyTesting\ctx1Data.h5'
ctx2Data = 'C:\Users\pb691\Desktop\MyTesting\ctx2Data.h5'
eulerTransPath = 'C:\Users\pb691\Desktop\MyTesting\eulerTransPath.h5'
eulerTransPathAfterRotation = 'C:\Users\pb691\Desktop\MyTesting\eulerTransPathAfterRotation.h5'
rotatedImage='C:\Users\pb691\Desktop\MyTesting\sdrotatedImage.h5'
bestEulerTransPath='C:\Users\pb691\Desktop\MyTesting\BestEulerTransPath.h5'
outTxPath='C:\Users\pb691\Desktop\MyTesting\outTxAfterRigidPath.h5'
quickSetVersorInitial='C:\Users\pb691\Desktop\MyTesting\quickSetVersorInitial.h5'
txiBeforeRigidPath='C:\Users\pb691\Desktop\MyTesting\TxBeforeRigid.h5'

# set output file paths
outputTransformInit='C:\Users\pb691\Desktop\MyTesting\output_transform_afterInit.h5'
outTxAfterRigidPath='C:\Users\pb691\Desktop\MyTesting\output_transform_afterRigid.h5'
fixedVolumePath='C:\Users\pb691\Desktop\MyTesting\FixedVolume.nrrd'
outVolumePath='C:\Users\pb691\Desktop\MyTesting\output_volume_afterInit.nrrd'
afterRigidPath='C:\Users\pb691\Desktop\MyTesting\output_volume_afterRigid.nrrd'
tx2outputPath='C:\Users\pb691\Desktop\MyTesting\output_transform_afterParamsPassing.h5'
# ___
rigid_versorPath='C:\Users\pb691\Desktop\MyTesting\Test__rigid_versor.h5'


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

axAngleRange = 12.0* one_degree
sagAngleRange = 12.0* one_degree

axStepSize = 3.0 * one_degree
sagStepSize = 3.0 * one_degree

# set current Metric value
initialMetricValue= Reg.MetricEvaluate(fixedVolume,movingVolume)

# initialize output Transform with Translation from CenteredTransformInitializer
bestEulerTrans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixedMask,movingMask,sitk.Euler3DTransform()))

for axAngle in np.arange(-axAngleRange,axAngleRange,axStepSize):
     for sagAngle in np.arange(-sagAngleRange,sagAngleRange,sagStepSize):

        eulerTrans.SetRotation(0,0,axAngle)
        print ('current axAngle : '+str(axAngle))
        Reg.SetInitialTransform(eulerTrans)
        currentMetricValue = Reg.MetricEvaluate(fixedVolume,movingVolume)

        if currentMetricValue < initialMetricValue:
            print ('new best axAngle : '+str(axAngle))
            bestEulerTrans.SetRotation(sagAngle,0,axAngle)
            print ('new best Euler Trans found at sagAngle '+str(sagAngle) +' and axAngle = '+str(axAngle))

sitk.WriteTransform(bestEulerTrans,bestEulerTransPath)

print ('tx Initial after Initialization and Euler Rotation :')
print ('_____________________________')
print ()
print bestEulerTrans
print ()

# PASS PARAMETERS
# _______________________________

quickSetVersor=sitk.VersorRigid3DTransform()
quickSetVersor.SetCenter(bestEulerTrans.GetCenter())
quickSetVersor.SetTranslation(bestEulerTrans.GetTranslation())
quickSetVersor.SetMatrix(bestEulerTrans.GetMatrix())

# make sure params are passed correctly:
print ('Transform after parameter passing :')
print ('_____________________________')
print ()
print quickSetVersor
print ()
sitk.WriteTransform(quickSetVersor,txiBeforeRigidPath)


# RIGID REGISTRATION PHASE
# _______________________________

Reg2=sitk.ImageRegistrationMethod()
Reg2.SetInitialTransform(quickSetVersor)

#Reg2.SetMetricAsCorrelation()
Reg2.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

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

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixedVolume)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(outTxAfterRigid)

out = resampler.Execute(movingVolume)
simg1 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
sitk.WriteImage(simg1,afterRigidPath)



# AFFINE REGISTRATION PHASE
# _______________________________


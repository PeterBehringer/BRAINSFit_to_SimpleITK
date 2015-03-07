import os
import SimpleITK as sitk
import slicer


def cshow3d(fixed, moving, transform=sitk.Transform(), title=""):
    """Show a tile of slices inline in the notebooks.
    """
    rescaled_fixed, deformed_rescaled_moving = prepare_resampled_values( fixed, moving, transform)
    myshow3d(sitk.Compose(deformed_rescaled_moving, deformed_rescaled_moving, sitk.Maximum(deformed_rescaled_moving,rescaled_fixed)), zslices=range(0, deformed_rescaled_moving.GetSize()[2],5), dpi=8,title=title)

def command_iteration(method) :
 print("{0:3} = {1:10.5f} : {2}".format(method.GetOptimizerIteration(),method.GetMetricValue(),method.GetOptimizerPosition()))

def prepare_resampled_values(fixed,moving,transform):
    """ Resample moving into space of fixed, and then rescale the
    intensity values of both to fit in a UInt8 data type.
    """
    f = sitk.ResampleImageFilter()
    f.SetTransform(transform)
    f.SetReferenceImage(fixed)
    print(f.GetTransform)
    deformed_moving = f.Execute(moving)
    fixed = sitk.Cast(sitk.RescaleIntensity(fixed),sitk.sitkUInt8)
    deformed_moving = sitk.Cast(sitk.RescaleIntensity(deformed_moving),sitk.sitkUInt8)
    return fixed, deformed_moving


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

# set output file

outputTransform = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutTrans_Rigid_1.h5'
outputVolume = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutVol_Rigid_1.nrrd'
outputTransform_Initializer = '/Users/peterbehringer/MyTesting/SimpleITK_Tests/OutTrans_Initializer_Rigid__1.h5'


# COMPUTE RIGID REGISTRATION
# _______________________________

# Initialize Rigid Transform by
transformForInitializer=sitk.VersorRigid3DTransform()
initializer=sitk.CenteredTransformInitializer(fixedMask,movingMask,transformForInitializer)
initialTrans = sitk.Transform(initializer)
initialTrans.SetFixedParameters(initializer.GetFixedParameters())

# write initial Transform
sitk.WriteTransform(initialTrans,outputTransform_Initializer)

# create ImageRegistrationMethod
Reg = sitk.ImageRegistrationMethod()

# set initial Transform for initial moving Transform
Reg.SetMovingInitialTransform(initialTrans)

# create rigid composite transform
compositeTrans=sitk.VersorRigid3DTransform()
Reg.SetInitialTransform(compositeTrans)

# Sparse Search
# Calculate similarity metric - how to do in simpleITK?

# Set the Metric as MMI
Reg.SetMetricAsMattesMutualInformation( 200 )
Reg.SetOptimizerAsRegularStepGradientDescent(learningRate =1.0,
                                           minStep=1e-3,
                                           numberOfIterations=250,
                                           gradientMagnitudeTolerance=1e-4,
                                           estimateLearningRate=R.Never)
Reg.SetOptimizerScales([1, 1, 1, 1.0/250, 1.0/250, 1.0/250])

# set linear Interpolator
Reg.SetInterpolator(sitk.sitkLinear)
Reg.SetMetricSamplingPercentage(0.5)
Reg.SetMetricSamplingStrategy(R.RANDOM)
Reg.SetShrinkFactorsPerLevel([1])
Reg.SetSmoothingSigmasPerLevel([0])

# execute
Reg.RemoveAllCommands()
RegAddCommand( sitk.sitkIterationEvent, lambda: command_iteration(Reg ))
outTx = R.Execute(sitk.Cast(fixedVolume,sitk.sitkFloat32), sitk.Cast(movingVolume,sitk.sitkFloat32))

# Write the transform
sitk.WriteTransform(outTx,outputTransform)

# perform volume transformation
print ('perform volume transformation')

slicer.mrmlScene.Clear(0)
tfmLogic = slicer.modules.transforms.logic()
vlmLogic=slicer.modules.volumes.logic()

volume=vlmLogic.AddArchetypeVolume('/Users/peterbehringer/MyImageData/ProstateRegistrationValidation/Images/Case1-t2ax-N4.nrrd','outputVolume1',4)
volumeNode=slicer.mrmlScene.GetNodesByClass('vtkMRMLScalarVolumeNode').GetItemAsObject(0)

transform=tfmLogic.AddTransform(outputTransform,slicer.mrmlScene)

volumeNode.SetAndObserveTransformNodeID(transform.GetID())
tfm = tfmLogic.AddTransform(outputTransform, slicer.mrmlScene)

tfmLogic.hardenTransform(volumeNode)

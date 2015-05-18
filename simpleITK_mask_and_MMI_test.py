import SimpleITK as sitk
from math import pi


def command_iteration(method) :
    if (method.GetOptimizerIteration()==0):
        print("Scales: ", method.GetOptimizerScales())
    print("{0:3} = {1:7.5f} : {2}".format(method.GetOptimizerIteration(),
                                           method.GetMetricValue(),
                                           method.GetOptimizerPosition()))


fixed_image=sitk.Image(200,200,10,sitk.sitkFloat32)
moving_image=sitk.Image(200,200,10,sitk.sitkFloat32)
fixed_image_mask=sitk.Image(200,200,10,sitk.sitkFloat32)
moving_image_mask=sitk.Image(200,200,10,sitk.sitkFloat32)

# _____________________________________________________________________________________________________________________ #

# create gaussian source for fixed_image
source=sitk.GaussianImageSource()
source.SetSigma(10.0)
origin=[0.0,0.0,0.0]
source.SetOrigin(origin)
fixed_image=source.Execute()

# _____________________________________________________________________________________________________________________ #

# create gaussian source for moving_image
source2=sitk.GaussianImageSource()
source2.SetSigma(10.0)
origin=[0.0,0.0,0.0]
source2.SetOrigin(origin)
moving_image=source2.Execute()

# set a translation by x=y=10
translation=[10.0,10.0,0]
transform=sitk.TranslationTransform(moving_image.GetDimension(),translation)

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixed_image)
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(1)
resampler.SetTransform(transform)
moving_image=resampler.Execute(moving_image)

# _____________________________________________________________________________________________________________________ #

# create spherical moving_image_mask

threshold=sitk.BinaryThresholdImageFilter()
threshold.SetLowerThreshold(150)
moving_image_mask=threshold.Execute(moving_image)

# _____________________________________________________________________________________________________________________ #

# create spherical fixed_image_mask
source3=sitk.GaussianImageSource()
source3.SetSigma(10.0)
origin=[0.0,0.0,0.0]
source3.SetOrigin(origin)
fixed_image_mask=source3.Execute()

threshold=sitk.BinaryThresholdImageFilter()
threshold.SetLowerThreshold(150)
fixed_image_mask=threshold.Execute(fixed_image_mask)

# _____________________________________________________________________________________________________________________ #
# perform Registration


# Initialize ImageRegistrationMethod()
Reg=sitk.ImageRegistrationMethod()
Reg.SetMetricFixedMask(fixed_image_mask)
Reg.SetMetricMovingMask(moving_image_mask)

# Reg.SetMetricAsCorrelation()
Reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins = 50)

Reg.SetInterpolator(sitk.sitkLinear)
Reg.SetMetricSamplingPercentage(.9)
Reg.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.RANDOM)


# Set Rotation Parameters for ExhaustiveOptimizer
# =========================================================================

# def 1 degree in rad
one_degree=1.0*pi/180.0

# search in neighbourhood from -12degree to +12degree
angleRange = 12.0* one_degree

# set step size to 3degree
stepSize = 3.0 * one_degree

# samples per axis
sample_per_axis=round((angleRange*2)//stepSize)

# set the number of samples from zero in +- degree direction in each dimension, -> int not double
# [rotationX, rotationY, rotationZ, translationX, translationY, translationZ]
Reg.SetOptimizerAsExhaustive([4,0,4,0,0,0])

# set the step size for each dimension
# [rotationScaleX, rotationScaleY, rotationScaleZ, translationScaleX, translationScaleY, translationScaleZ]
Reg.SetOptimizerScales([stepSize,stepSize,stepSize,1.0,1.0,1.0])

# Initialize the transform with a translation and the center of rotation from the moments of intensity.
# =========================================================================

# Set the Euler3DTransform
euler_trans=sitk.Euler3DTransform(sitk.CenteredTransformInitializer(fixed_image_mask,moving_image_mask,sitk.Euler3DTransform()))

# Set, Execute & write
Reg.SetInitialTransform(euler_trans,inPlace=True)
Reg.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(Reg))
Reg.Execute(fixed_image, moving_image)

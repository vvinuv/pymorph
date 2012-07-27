import pyfits
import ndimage
import math
import numpy 
import os
import _ni_support

def _minmax(coor, minc, maxc):
    if coor[0] < minc[0]:
        minc[0] = coor[0]
    if coor[0] > maxc[0]:
        maxc[0] = coor[0]
    if coor[1] < minc[1]:
        minc[1] = coor[1]
    if coor[1] > maxc[1]:
        maxc[1] = coor[1]
    return minc, maxc


def rotate_modi(input, angle, xcntr, ycntr, axes = (1, 0), reshape = False,
           output_type = None, output = None, order = 1,
           mode = 'nearest', cval = 0.0, prefilter = False):
    """Rotate an array.

    The array is rotated in the plane defined by the two axes given by the
    axes parameter using spline interpolation of the requested order. The
    angle is given in degrees. Points outside the boundaries of the input
    are filled according to the given mode. If reshape is true, the output
    shape is adapted so that the input array is contained completely in
    the output. The parameter prefilter determines if the input is pre-
    filtered before interpolation, if False it is assumed that the input
    is already filtered.
    """
    xcntr = xcntr + 1.0 #This function needs iraf like center
    ycntr = ycntr + 1.0
    input = numpy.asarray(input)
    axes = list(axes)
    rank = input.ndim
    if axes[0] < 0:
        axes[0] += rank
    if axes[1] < 0:
        axes[1] += rank
    if axes[0] < 0 or axes[1] < 0 or axes[0] > rank or axes[1] > rank:
        raise RuntimeError, 'invalid rotation plane specified'
    if axes[0] > axes[1]:
        axes = axes[1], axes[0]
    angle = numpy.pi / 180 * angle
    m11 = math.cos(angle)
    m12 = math.sin(angle)
    m21 = -math.sin(angle)
    m22 = math.cos(angle)
    matrix = numpy.array([[m11, m12],
                             [m21, m22]], dtype = numpy.float64)
    iy = input.shape[axes[0]]
    ix = input.shape[axes[1]]
    if reshape:
        mtrx = numpy.array([[ m11, -m21],
                               [-m12,  m22]], dtype = numpy.float64)
        minc = [0, 0]
        maxc = [0, 0]
        coor = numpy.dot(mtrx, [0, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, 0])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        oy = int(maxc[0] - minc[0] + 0.5)
        ox = int(maxc[1] - minc[1] + 0.5)
    else:
        oy = input.shape[axes[0]]
        ox = input.shape[axes[1]]
    offset = numpy.zeros((2,), dtype = numpy.float64)
#    offset[0] = float(oy) / 2.0 - 0.5
#    offset[1] = float(ox) / 2.0 - 0.5
    offset[1] = xcntr - 1.0 #x center
    offset[0] = ycntr - 1.0 #y center
    offset = numpy.dot(matrix, offset)
    tmp = numpy.zeros((2,), dtype = numpy.float64)
#    tmp[0] = float(iy) / 2.0 - 0.5
#    tmp[1] = float(ix) / 2.0 - 0.5
    tmp[1] = xcntr - 1.0
    tmp[0] = ycntr - 1.0
    offset = tmp - offset
#    offset = [406, 566]
    output_shape = list(input.shape)
    output_shape[axes[0]] = oy
    output_shape[axes[1]] = ox
    output_shape = tuple(output_shape)
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    ndimage.affine_transform(input, matrix, offset, output_shape, None, output,
                         order, mode, cval, prefilter)
    return return_value

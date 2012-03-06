from skymaps import SkyDir
import numpy as np

def rotate_north(skydir,target,anti=False):
    """ Transformation that will rotate target to celestial north """

    axis=SkyDir(target.ra()-90,0)
    theta=np.radians(90-target.dec())
    if anti: theta*=-1

    newdir=SkyDir(skydir.ra(),skydir.dec())
    newdir().rotate(axis(),theta)
    return newdir

def anti_rotate_north(skydir,target):
    return rotate_north(skydir,target,anti=True)

def rotate_equator(skydir,target,anti=False):
    """ Rotate skydir such that target would be rotated 
        to the celestial equator. 
        
        A few simple tests:

            >>> a = SkyDir(0,0)
            >>> b = SkyDir(30,30)
            >>> rotate_equator(a, b)
            SkyDir(330.000,-30.000)
            >>> anti_rotate_equator(a, b)
            SkyDir(30.000,30.000)

        There was previously a bug when the 'target' was the equator.
        I think it is fixed now:

            >>> rotate_equator(b, a)
            SkyDir(30.000,30.000)
            >>> anti_rotate_equator(b, a)
            SkyDir(30.000,30.000)


        Another test: 

            >>> sd = SkyDir(-.2,-.2)
            >>> target = SkyDir(5,5)
            >>> l=rotate_equator(sd, target)
            >>> print '%.2f, %.2f' % (l.ra(), l.dec())
            354.80, -5.20
            >>> l=anti_rotate_equator(sd, target)
            >>> print '%.2f, %.2f' % (l.ra(), l.dec())
            4.80, 4.80
    """
    if np.allclose([target.ra(), target.dec()], [0,0]):
        return skydir


    equator=SkyDir(0,0)

    axis=target.cross(equator)

    theta=equator.difference(target)
    if anti: theta*=-1

    newdir=SkyDir(skydir.ra(),skydir.dec())
    newdir().rotate(axis(),theta)
    return newdir

def anti_rotate_equator(skydir,target):
    """ Performs the opposite rotation of the rotate_equator function. """
    return rotate_equator(skydir,target,anti=True)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

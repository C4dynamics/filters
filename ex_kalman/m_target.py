# import math 
from m_params import params
import numpy as np
from builtins import object

class target(object):
    """ 
        data point 
        vertically falling target
        ballistic coeficient 
    """
    z = 0
    vz = 0
    az = 0
    beta = 0
    data = np.array([[0, 0, 0, 0]]) # r, v, a, beta 
    
    def __init__(obj, x0, beta):
        obj.z, obj.vz, obj.az = x0
        obj.beta = beta
        obj.data = np.array([[obj.z, obj.vz, obj.az, obj.beta]]) #


    @staticmethod
    def deriv(y, t, tgt): 
        """ 
        run the target equations of motion. 
        """
        if y.ndim == 1:
            zT = y[0]
            vzT = y[1]
        else:
            zT = y[:, 0]
            vzT = y[:, 1]
            
        dzT = vzT
        dvzT = .0034 * np.exp(-zT.astype('float') / 22000 / params.ft2m) * params.g * vzT**2 / 2 / tgt.beta - params.g
            
        return dzT, dvzT

 

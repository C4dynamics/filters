import importlib
import numpy as np
# import sys, os
# sys.path.append(os.getcwd() + '\..\ex_kalman')
# sys.path.append(os.getcwd() + '\..\luenberger')

from m_params import c_params

import m_filtertype 
importlib.reload(m_filtertype)
from m_filtertype import e_filtertype

import m_luenberger 
importlib.reload(m_luenberger)
from m_luenberger import c_luenberger

import m_e_kalman
importlib.reload(m_e_kalman)
from m_e_kalman import c_e_kalman


# np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)                 

class c_radar:
  """ 
    an electromagnetic radar 
    meausres target vertical position.
    vk: noise variance  
  """
  
  # 
  # state vector
  #   r position (in one dimension range)
  #   v velocity 
  #   b ballistic coefficient 
  ##
  r = 0 # measured range 
  # vx = 0 # velocity 
  
  
  ts = 50e-3                        # sample time 
  q33 = 0
  data = np.array([[0, 0, 0, 0]])   # time, z target, v target, beta target 
  
  # luenberger gains 
  L = 0
  
  # 
  # seeker errors 
  ##
  bias = 0 # deg
  sf = 1 # scale factor error -10%
  noise = np.sqrt(500 * c_params.ft2m) # 1 sigma 
  # misalignment = 0.01 # deg

  def __init__(obj, x0, filtype):
      '''
       initial estimate: 
       100025 ft    25ft error
       6150 ft/s    150ft/s error
       800 lb/ft^2  300lb/ft^2 
      '''
      if hasattr(filtype, 'value'):
        if filtype.value == e_filtertype.ex_kalman.value:
          obj.ifilter = c_e_kalman(x0, [obj.noise, 141.42 * c_params.ft2m, 300 * c_params.lbft2kgm], obj.ts)
        elif filtype.value == e_filtertype.luenberger.value:
          # linear model
          beta0 = x0[2]
          A = np.array([[0, 1, 0], [0, -np.sqrt(2 * 0.0034 * c_params.g / beta0), -c_params.g / beta0], [0, 0, 0]])
          b = np.array([[0], [0], [0]])
          c = np.array([1, 0, 0])
          obj.ifilter = c_luenberger(A, b, c)
      
      
      obj.data[0, :] =  np.insert(x0, 0, 0)  
      
      
  def measure(obj, x):
    # 
    # apply errors
    ##
    obj.r = x * obj.sf + obj.bias + obj.noise * np.random.randn(1) 
  
      
  def filter(obj, t):
    
    # check that t is at time of operation 
    
    if (1000 * t) % (1000 * obj.ts) > 1e-6 or t <= 0:
      return
    
    
    # 
    # prepare kalman matrices 
    ##
    
    # print(t)
    rho = .0034 * np.exp(-obj.ifilter.x[0, 0] / 22000 / c_params.ft2m)
    f21 = -rho * c_params.g * obj.ifilter.x[1, 0]**2 / 44000 / obj.ifilter.x[2, 0] 
    f22 =  rho * c_params.g * obj.ifilter.x[1, 0] / obj.ifilter.x[2, 0]
    f23 = -rho * c_params.g * obj.ifilter.x[1, 0]**2 / 2 / obj.ifilter.x[2, 0]**2 
    
    Phi = np.array([[1, obj.ifilter.tau, 0], 
                    [f21 * obj.ifilter.tau, 1 + f22 * obj.ifilter.tau, f23 * obj.ifilter.tau], 
                    [0, 0, 1]])
    
    Q   = np.array([[0, 0, 0], 
                    [0, obj.q33 * f23**2 * obj.ifilter.tau**3 / 3, obj.q33 * f23 * obj.ifilter.tau**2 / 2], 
                    [0, obj.q33 * f23 * obj.ifilter.tau**2 / 2, obj.q33 * obj.ifilter.tau]])

    f = lambda w: c_radar.system_model(w)

    #
    # prepare luenberger matrices
    ##

    ''' predict '''
    obj.ifilter.predict(f, Phi, Q)
    ''' correct '''    
    x = obj.ifilter.correct(f, obj.r)  
 
 
    #
    # store results 
    ##
 
    # obj.data = np.concatenate((obj.data, np.expand_dims(np.insert(x, 0, t), axis = 0)), axis = 0)  
    obj.data = np.vstack((obj.data, np.insert(x, 0, t))).copy()
    
    return x

    

  @staticmethod
  def system_model(x):
    dx = np.zeros((len(x), 1))
    dx[0, 0] = x[1, 0]
    dx[1, 0] = .0034 * np.exp(-x[0, 0].astype('float') / 22000 / c_params.ft2m) * c_params.g * x[1, 0]**2 / 2 / x[2, 0] - c_params.g
    dx[2, 0] = 0
    return dx

 

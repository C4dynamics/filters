from m_e_kalman import e_kalman
from m_params import params
import numpy as np
# np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)                 

class radar:
  """ 
    an electromagnetic radar 
    meausres target position.
    vk: measure's error variance  
  """
  

  vk = (22.361 * params.ft2m)      # measurement variance:
  ts = 50e-3                # sample time 
  q33 = 0
  data = np.array([[0, 0, 0, 0]]) # time, z target, v target, beta target 

  def __init__(obj, x0):
      '''
       initial estimate: 
       100025 ft    25ft error
       6150 ft/s    150ft/s error
       800 lb/ft^2  300lb/ft^2 
      '''
      obj.filter = e_kalman(x0, [obj.vk, 141.42 * params.ft2m, 300 * params.lbft2kgm], obj.ts)
      obj.data[0, :] =  np.insert(x0, 0, 0)  # np.concatenate(([0], x0)) # np.reshape(x0, (-1, len(x0)))
      
  def measure(obj, t, rng_in):
    
    # check that t is at time of operation 
    
    if (1000 * t) % (1000 * obj.ts) > 1e-6 or t <= 0:
      return
    
    # print(t)
    rho = .0034 * np.exp(-obj.filter.x[0, 0] / 22000 / params.ft2m)
    f21 = -rho * params.g * obj.filter.x[1, 0]**2 / 44000 / obj.filter.x[2, 0] 
    f22 =  rho * params.g * obj.filter.x[1, 0] / obj.filter.x[2, 0]
    f23 = -rho * params.g * obj.filter.x[1, 0]**2 / 2 / obj.filter.x[2, 0]**2 
    
    Phi = np.array([[1, obj.filter.tau, 0], 
                    [f21 * obj.filter.tau, 1 + f22 * obj.filter.tau, f23 * obj.filter.tau], 
                    [0, 0, 1]])
    
    Q   = np.array([[0, 0, 0], 
                    [0, obj.q33 * f23**2 * obj.filter.tau**3 / 3, obj.q33 * f23 * obj.filter.tau**2 / 2], 
                    [0, obj.q33 * f23 * obj.filter.tau**2 / 2, obj.q33 * obj.filter.tau]])


    ''' predict '''
    f = lambda w: radar.system_model(w)
    obj.filter.predict(f, Phi, Q)


    ''' correct '''    
    x = obj.filter.correct(f, rng_in)  
    obj.data = np.concatenate((obj.data, np.expand_dims(np.insert(x, 0, t), axis = 0)), axis = 0)  
    return x

  @staticmethod
  def system_model(x):
    dx = np.zeros((len(x), 1))
    dx[0, 0] = x[1, 0]
    dx[1, 0] = .0034 * np.exp(-x[0, 0].astype('float') / 22000 / params.ft2m) * params.g * x[1, 0]**2 / 2 / x[2, 0] - params.g
    dx[2, 0] = 0
    return dx

 

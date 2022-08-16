import math
from skimage import measure

class microstrip_sub:
  def __init__(sub, t_metal, metalCond, t_sub, relPerm, freq):
    """ define the microstrip substrate
    Args:
    t_metal = conductor thickness in mils (e.g., 1.4 mil)
    metalCond = conductivity of the metal in S/m
    t_sub = substrate height between conductor and ground plane in mils (e.g., 19 mils)
    relPerm = relative permittivity of the substrate material
    freq = desired operating frequency in Hz (e.g., 2.4e9)
    """    
    sub.t_metal = t_metal
    sub.sigMetal = metalCond
    sub.t_sub = t_sub
    sub.relPerm = relPerm
    sub.freq = freq

class microstrip_calc:

  def ee_HandJ(sub,
               u: float):
    A  = 1.0 + (1.0/49.0)*math.log((math.pow(u,4.0) + math.pow((u/52.0),2.0))/(math.pow(u,4.0) + 0.432)) \
            + (1.0/18.7)*math.log(1.0 + math.pow((u/18.1),3.0));
    B  = 0.564*math.pow(((sub.relPerm-0.9)/(sub.relPerm+3.0)),0.053);
    ee = (sub.relPerm+1.0)/2.0 + ((sub.relPerm-1.0)/2.0)*math.pow((1.0 + 10.0/u),(-A*B));
    return ee

  def Z0_HandJ(u: float):
    LIGHTSPEED = 2.99792458e8;
    FREESPACEZ0 = 4.0*math.pi*1.0e-7*LIGHTSPEED;
    F = 6.0 + (2.0*math.pi - 6.0)*math.exp(-1*math.pow((30.666/u),0.7528));
    z01 = (FREESPACEZ0/(2*math.pi))*math.log(F/u + math.sqrt(1.0 + math.pow((2/u),2.0)));
    return z01
  
  def calcMicrostrip(sub,
                     width: float,
                     length: float):
    """ calculate the impedance of a microstrip line based on dimensions
    Args:
    sub --> defined by class microstrip_sub, includes, metal thickness, substrate 
            thickness, relative permittivity and frequency of operation
    width = conductor width in mils (e.g., 38 mils)
    length = physical length of conductor (e.g., 1000 mils)
    """                
    u = width/sub.t_sub # ratio of the conductor width to height
    if sub.t_metal > 0:
      t = sub.t_metal/sub.t_sub
      du1 = (t*math.log(1.0+4.0*math.exp(1)/t/math.pow(1.0/math.tanh(math.sqrt(6.517*u)),2.0)))/math.pi  # from Hammerstad and Jensen: ACCURATE 
           #MODELS FOR MICROSTRIP COMPUTER-AIDED DESIGN
      u1 = u + du1
      dur = (du1/2)*(1.0+1.0/\
           math.cosh(math.sqrt(sub.relPerm-1))); # from Hammerstad and Jensen
      ur = u + dur
    else:
      u1 = u
      ur = u
    Y = microstrip_calc.ee_HandJ(sub, ur)
    Z0 = microstrip_calc.Z0_HandJ(ur)/math.sqrt(Y)    
    ereff0 = Y*math.pow(microstrip_calc.Z0_HandJ(u1)/microstrip_calc.Z0_HandJ(ur),2.0);
    fn = 1e-9*.0254*sub.freq*sub.t_sub#/1e7 # 1e-9*.0254 converts to GHz*mm
    P1 = 0.27488 + (0.6315 + 0.525 / math.pow(1 + 0.157*fn,20) )*u \
         - 0.065683*math.exp(-8.7513*u);
    P2 = 0.33622*(1 - math.exp(-0.03442*sub.relPerm));
    P3 = 0.0363*math.exp(-4.6*u)*(1 - math.exp(-math.pow((fn / 38.7),4.97)));
    P4 = 1 + 2.751*( 1 -  math.exp(-math.pow((sub.relPerm/15.916),8)));
    P = P1*P2*math.pow((0.1844 + P3*P4)*fn,1.5763);
    ereff = (sub.relPerm*P+ereff0)/(1+P); # equavlent relative dielectric constant
    #fn = 1e-9*.00254*sub.freq*sub.t_sub#/1e6 # 1e-9*.00254 converts to GHz*cm
    R1 = 0.03891*(math.pow(sub.relPerm,1.4));
    R2 = 0.267*(math.pow(u,7.0));
    R3 = 4.766*math.exp(-3.228*(math.pow(u,0.641)));
    R4 = 0.016 + math.pow((0.0514*sub.relPerm),4.524);
    R5 = math.pow(fn/28.843,12.0);
    R6 = 22.20*(math.pow(u,1.92));
    R7 = 1.206 - 0.3144*math.exp(-R1)*(1 - math.exp(-R2));
    R8 = 1.0 + 1.275*(1.0 -  math.exp(-0.004625*R3*math.pow(sub.relPerm,1.674)\
         *math.pow(fn/18.365,2.745)));
    R9 = (5.086*R4*R5/(0.3838 + 0.386*R4))*(math.exp(-R6)/(1 + 1.2992*R5));
    R9 = R9 * math.pow(sub.relPerm-1,6)/(1 + 10*math.pow(sub.relPerm-1,6));
    R10 = 0.00044*math.pow(sub.relPerm,2.136) + 0.0184;
    R11 = math.pow((fn/19.47),6)/(1 + 0.0962*math.pow(fn/19.47,6));
    R12 = 1 / (1 + 0.00245*math.pow(u,2));
    R13 = 0.9408*(math.pow( ereff,R8)) - 0.9603;
    R14 = (0.9408 - R9)*(math.pow(ereff0,R8))-0.9603;
    R15 = 0.707*R10*math.pow(fn/12.3,1.097);
    R16 = 1 + 0.0503*math.pow(sub.relPerm,2)*R11*(1 - math.exp(-math.pow(u/15,6)));
    R17 = R7*(1 - 1.1241*(R12/R16)*math.exp(-0.026*math.pow(fn,1.15656)-R15));
    Zc = Z0*(math.pow((R13/R14),R17)) # characteristic impedance
    return Zc, ereff
    
  def synthMicrostrip(sub,
                      imp: float,
                      eLength: float):
      """ calculate the dimensions of a microstrip line output the width and 
          length of a microstrip line in units of mils based on desired
          impedance and electrical length
      Args:
      sub --> defined by class microstrip_sub, includes, metal thickness, substrate 
              thickness, relative permittivity and frequency of operation
      imp = desired characteristic impedance (e.g., 50 Ohms)
      eLength = deisred electrical length (e.g., 90 degrees)
      """
      eps_0 = 8.854187817e-12
      mu_0 = 4*math.pi*1e-7
      c = 1/math.sqrt(eps_0*mu_0)
    
      # define some constants
      mur = 1 # relative permeability
      cond = sub.sigMetal # conductivity of metal (assumes copper)
      mu = mur * mu_0
    
      lambda_0 = c/sub.freq # wavelength in freespace in m
      lx = 400 # initial length of line in mils
      wmin = 4 # minimum width of conductor in mils
      wmax = 200 # maximum width of conductor in mils
    
      abstol = 1.0e-6
      reltol = 0.1e-6
      maxiters = 50
    
      A = ((sub.relPerm - 1)/(sub.relPerm + 1)) * (0.226 + 0.121/sub.relPerm) \
          + (math.pi/377)*math.sqrt(2*(sub.relPerm+1))*imp
      w_h = 4/(0.5*math.exp(A) - math.exp(-A))
      if w_h > 2:
        B = math.pi*377/(2*imp*math.sqrt(sub.relPerm))
        w_h = (2/math.pi)*(B - 1 - math.log(2*B - 1) + ((sub.relPerm-1)/(2*sub.relPerm))\
              *(math.log(B-1) + 0.293 - 0.517/sub.relPerm));
      
      wx = sub.t_sub*w_h
    
      if wx >= wmax:
        wx = 0.95*wmax
      
      if wx <= wmin:
        wx = wmin
      
      wold = 1.01*wx
      zold, erold = microstrip_calc.calcMicrostrip(sub, wold, lx);
   
      if zold < imp:
        wmax = wold
      else:
        wmin = wold
    
      iters = 0
      done = 0
    
      while done == 0:
        iters = iters + 1;
        z0, er0 = microstrip_calc.calcMicrostrip(sub, wx, lx);
        if z0 < imp:
          wmax = wx
        else:
          wmin = wx
        if abs(z0-imp) < abstol:
          done = 1
        elif abs(wx-wold) < reltol:
          done = 1
        elif iters >= maxiters:
          done = 1
        else:
          dzdw = (z0 - zold)/(wx-wold)
          wold = wx
          zold = z0
          wx = wx - (z0-imp)/dzdw
          if (wx > wmax) or (wx < wmin):
            wx = (wmin + wmax) / 2
      zD, ereff = microstrip_calc.calcMicrostrip(sub, wx, lx);
    
      v = c/math.sqrt(ereff)
      l = 39370.0787*eLength/360*v/sub.freq # 39370.0787 converts m to mil
      return wx, l

def findConnectivity(arr, pixelSize, ports, sides, portPos):
  labeled = measure.label(arr, background=False, connectivity=1)
  if ports == 2:
    c13 = False
    c14 = False
    c23 = False
    c24 = False
    c34 = False
    if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
      labeled[int(portPos[2]/pixelSize)-1,int(portPos[3]/pixelSize)]:
      c12 = True
    else:
      c12 = False
  if ports == 3:
    if sides == 2:
      c14 = False
      c24 = False
      c34 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)]:
        c12 = True
      else:
        c12 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c13 = True
      else:
        c13 = False
      if labeled[int(portPos[2]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c23 = True
      else:
        c23 = False
    if sides == 3:
      c14 = False
      c24 = False
      c34 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[2]/pixelSize)-1,int(portPos[3]/pixelSize)]:
        c12 = True
      else:
        c12 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize),int(portPos[5]/pixelSize)-1]:
        c13 = True
      else:
        c13 = False
      if labeled[int(portPos[2]/pixelSize)-1,int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize),int(portPos[5]/pixelSize)-1]:
        c23 = True
      else:
        c23 = False
  if ports == 4:
    if sides == 2:
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)]:
        c12 = True
      else:
        c12 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c13 = True
      else:
        c13 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[6]/pixelSize)-1,int(portPos[7]/pixelSize)]:
        c14 = True
      else:
        c14 = False
      if labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c23 = True
      else:
        c23 = False
      if labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)] == \
        labeled[int(portPos[6]/pixelSize)-1,int(portPos[7]/pixelSize)]:
        c24 = True
      else:
        c24 = False
      if labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)] == \
        labeled[int(portPos[6]/pixelSize)-1,int(portPos[7]/pixelSize)]:
        c34 = True
      else:
        c34 = False
    if sides == 4:
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)-1]:
        c12 = True
      else:
        c12 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c13 = True
      else:
        c13 = False
      if labeled[int(portPos[0]/pixelSize),int(portPos[1]/pixelSize)] == \
        labeled[int(portPos[6]/pixelSize),int(portPos[7]/pixelSize)-1]:
        c14 = True
      else:
        c14 = False
      if labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)-1] == \
        labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)]:
        c23 = True
      else:
        c23 = False
      if labeled[int(portPos[2]/pixelSize),int(portPos[3]/pixelSize)-1] == \
        labeled[int(portPos[6]/pixelSize),int(portPos[7]/pixelSize)-1]:
        c24 = True
      else:
        c24 = False
      if labeled[int(portPos[4]/pixelSize)-1,int(portPos[5]/pixelSize)] == \
        labeled[int(portPos[6]/pixelSize),int(portPos[7]/pixelSize)-1]:
        c34 = True
      else:
        c34 = False
  return c12, c13, c14, c23, c24, c34

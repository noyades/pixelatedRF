import math
import numpy as np
from skimage import measure

class MicrostripSub:
  def __init__(self, t_metal: float, metal_cond: float, t_sub:float, rel_perm: float, freq: float):
      """ define the microstrip substrate
      Args:
          t_metal = conductor thickness in mils (e.g., 1.4 mil)
          metal_cond = conductivity of the metal in S/m
          t_sub = substrate height between conductor and ground plane in mils (e.g., 19 mils)
          rel_perm = relative permittivity of the substrate material
          freq = desired operating frequency in Hz (e.g., 2.4e9)
      """    
      self.t_metal = t_metal
      self.sig_metal = metal_cond
      self.t_sub = t_sub
      self.rel_perm = rel_perm
      self.freq = freq

class MicrostripCalc:
  @staticmethod
  def ee_handj(sub: MicrostripSub, u: float) -> float:
      """
      Calculate the effective permittivity using Hammerstad and Jensen's equations.

      Args:
          sub (MicrostripSub): The microstrip substrate.
          u (float): The ratio of the strip width to the substrate height.

      Returns:
          float: The effective permittivity.
      """    
      A  = (
           1.0 
           + (1.0 / 49.0) * math.log((u**4 + (u / 52.0)**2) / (u**4 + 0.432))
           + (1.0 / 18.7) * math.log(1.0 + (u / 18.1)**3)
      )
    
      B  = 0.564 * ((sub.rel_perm - 0.9) / (sub.rel_perm + 3.0))**0.053
  
      ee = (
           (sub.rel_perm + 1.0) / 2.0 + ((sub.rel_perm - 1.0) / 2.0)
           * (1.0 + 10.0/u)**(-A * B)
      )
      return ee

  @staticmethod
  def z0_handj(u: float):
      """
        Calculate the characteristic impedance using Hammerstad and Jensen's equations.

        Args:
            u (float): The ratio of the strip width to the substrate height.

        Returns:
            float: The characteristic impedance.
        """
      LIGHTSPEED = 2.99792458e8  # Speed of light in vacuum (m/s)
      FREESPACEZ0 = 4.0 * math.pi *1.0e-7 * LIGHTSPEED # Characteristic impedance of free space (Ohms)
      
      # Calculate the correction factor F
      F = 6.0 + (2.0 * math.pi - 6.0) * math.exp(-1 * math.pow((30.666 / u), 0.7528))
      
      # Calculate the characteristic impedance Z0
      z01 = (FREESPACEZ0 / (2 * math.pi)) * math.log(
          F / u + math.sqrt(1.0 + math.pow((2 / u), 2.0))
      )
      return z01
  
  @staticmethod
  def calc_microstrip(sub: MicrostripSub, width: float, length: float) -> float:
      """
        Calculate the impedance of a microstrip line based on dimensions.

        Args:
            sub (MicrostripSub): The microstrip substrate.
            width (float): Conductor width in mils (e.g., 38 mils).
            length (float): Physical length of conductor (e.g., 1000 mils).

        Returns:
            float: The impedance of the microstrip line.
      """
      u = width / sub.t_sub # ratio of the conductor width to height
      
      if sub.t_metal > 0:
          t = sub.t_metal/sub.t_sub
          # Calculate the correction factor for the conductor thickness
          du1 = (t * math.log(1.0 + 4.0 * math.exp(1) / t / math.pow(1.0 / math.tanh(math.sqrt(6.517 * u)), 2.0))) / math.pi  # from Hammerstad and Jensen: ACCURATE #MODELS FOR MICROSTRIP COMPUTER-AIDED DESIGN
          u1 = u + du1
          # Calculate the correction factor for the relative permittivity
          dur = (du1 / 2) * (1.0 + 1.0 / math.cosh(math.sqrt(sub.rel_perm - 1))) # from Hammerstad and Jensen
          ur = u + dur
      else:
          u1 = u
          ur = u
      # Calculate the effective permittivity    
      Y = MicrostripCalc.ee_handj(sub, ur)
      # Calculate the characteristic impedance
      Z0 = MicrostripCalc.z0_handj(ur) / math.sqrt(Y)
      # Calculate the effective relative permittivity
      ereff0 = Y * math.pow(MicrostripCalc.z0_handj(u1) / MicrostripCalc.z0_handj(ur), 2.0)
      # Calculate the normalized frequency
      fn = 1e-9 * .0254 * sub.freq * sub.t_sub #/1e7 # 1e-9*.0254 converts to GHz*mm

      # Calculate the correction factors P1, P2, P3, and P4
      P1 = 0.27488 + (0.6315 + 0.525 / math.pow(1 + 0.157 * fn, 20)) * u - 0.065683 * math.exp(-8.7513 * u)
      P2 = 0.33622 * (1 - math.exp(-0.03442 * sub.rel_perm))
      P3 = 0.0363 * math.exp(-4.6 * u) * (1 - math.exp(-math.pow((fn / 38.7), 4.97)))
      P4 = 1 + 2.751 * (1 - math.exp(-math.pow((sub.rel_perm / 15.916), 8)))
    
      # Calculate the final correction factor P
      P = P1 * P2 * math.pow((0.1844 + P3 * P4) * fn, 1.5763)
      
      ereff = (sub.rel_perm * P + ereff0) / (1 + P) # equavlent relative dielectric constant
    
      R1 = 0.03891 * (math.pow(sub.rel_perm, 1.4))
      R2 = 0.267 * (math.pow(u, 7.0))
      R3 = 4.766 * math.exp(-3.228 * (math.pow(u, 0.641)))
      R4 = 0.016 + math.pow((0.0514 * sub.rel_perm), 4.524)
      R5 = math.pow(fn / 28.843, 12.0)
      R6 = 22.20 * (math.pow(u, 1.92))
      R7 = 1.206 - 0.3144 * math.exp(-R1) * (1 - math.exp(-R2))
      R8 = 1.0 + 1.275 * (1.0 -  math.exp(-0.004625 * R3 * math.pow(sub.rel_perm, 1.674) * math.pow(fn / 18.365, 2.745)));
      R9 = (5.086 * R4 * R5 / (0.3838 + 0.386 * R4))*(math.exp(-R6) / (1 + 1.2992 * R5)) * math.pow(sub.rel_perm-1, 6) / (1 + 10 * math.pow(sub.rel_perm - 1, 6))
      R10 = 0.00044 * math.pow(sub.rel_perm, 2.136) + 0.0184
      R11 = math.pow((fn / 19.47), 6) / (1 + 0.0962 * math.pow(fn / 19.47, 6))
      R12 = 1 / (1 + 0.00245 * math.pow(u, 2))
      R13 = 0.9408 * (math.pow(ereff, R8)) - 0.9603
      R14 = (0.9408 - R9) * (math.pow(ereff0, R8)) - 0.9603
      R15 = 0.707 * R10 * math.pow(fn / 12.3, 1.097)
      R16 = 1 + 0.0503 * math.pow(sub.rel_perm, 2) * R11 * (1 - math.exp(-math.pow(u / 15, 6)))
      R17 = R7 * (1 - 1.1241 * (R12 / R16) * math.exp(-0.026 * math.pow(fn, 1.15656) - R15))
      Zc = Z0 * (math.pow((R13 / R14), R17)) # characteristic impedance
      return Zc, ereff

  @staticmethod  
  def synth_microstrip(sub: MicrostripSub, imp: float, eLength: float) -> float:
      """ calculate the dimensions of a microstrip line output the width and 
          length of a microstrip line in units of mils based on desired
          impedance and electrical length
      Args:
      sub --> defined by class MicrostripSub, includes, metal thickness, substrate 
              thickness, relative permittivity and frequency of operation
      imp = desired characteristic impedance (e.g., 50 Ohms)
      eLength = deisred electrical length (e.g., 90 degrees)
      """
      eps_0 = 8.854187817e-12
      mu_0 = 4 * math.pi * 1e-7
      c = 1 / math.sqrt(eps_0 * mu_0)
    
      # define some constants
      mur = 1 # relative permeability
      cond = sub.sig_metal # conductivity of metal (assumes copper)
      mu = mur * mu_0
    
      lambda_0 = c / sub.freq # wavelength in freespace in m
      lx = 400 # initial length of line in mils
      wmin = 0.04 # minimum width of conductor in mils
      wmax = 500 # maximum width of conductor in mils
    
      abstol = 1.0e-6
      reltol = 0.1e-6
      maxiters = 50
    
      A = ((sub.rel_perm - 1) / (sub.rel_perm + 1)) * (0.226 + 0.121 / sub.rel_perm) + (math.pi / 377) * math.sqrt(2 * (sub.rel_perm + 1)) * imp
      w_h = 4 / (0.5 * math.exp(A) - math.exp(-A))
      
      if w_h > 2:
        B = math.pi * 377 / (2 * imp * math.sqrt(sub.rel_perm))
        w_h = (2 / math.pi) * (B - 1 - math.log(2 * B - 1) + ((sub.rel_perm - 1) / (2 * sub.rel_perm)) * (math.log(B - 1) + 0.293 - 0.517 / sub.rel_perm))
      
      wx = sub.t_sub * w_h
    
      if wx >= wmax:
        wx = 0.95 * wmax
      
      if wx <= wmin:
        wx = wmin
      
      wold = 1.01 * wx
      zold, erold = MicrostripCalc.calc_microstrip(sub, wold, lx)
   
      if zold < imp:
        wmax = wold
      else:
        wmin = wold
    
      iters = 0
      done = 0
    
      while done == 0:
        iters = iters + 1
        z0, er0 = MicrostripCalc.calc_microstrip(sub, wx, lx)
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
          dzdw = (z0 - zold) / (wx - wold)
          wold = wx
          zold = z0
          wx = wx - (z0 - imp) / dzdw
          if (wx > wmax) or (wx < wmin):
            wx = (wmin + wmax) / 2
      zD, ereff = MicrostripCalc.calc_microstrip(sub, wx, lx)
    
      v = c / math.sqrt(ereff)
      l = 39370.0787 * eLength / 360 * v /sub.freq # 39370.0787 converts m to mil
      return wx, l
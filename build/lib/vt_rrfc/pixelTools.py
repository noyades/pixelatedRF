
import numpy as np
import os
import time

class changePixels:

  def __init__(self, matrix_mask = None, relaxing_time = 0):
    self.__last_array = None
    self.__lastest_array = None
    self.relaxing_time = relaxing_time
    if (type(matrix_mask) != type(None)):
        self.matrix_mask = np.array(matrix_mask, dtype=np.int32)
    else:
        self.matrix_mask = matrix_mask

  def update(self, matrix):
    '''
    Update pixel region according to the new matrix. For the first time it is called, the pixels 
    will be created in the FDTD simulation CAD. In the following update process, it will 
    enable/disable correspoinding pixels.
    Parameters
    ----------
    matrix : numpy.array
      Array (values:0~1) that represent the pixels in the region.
    '''
    if (type(self.matrix_mask) != type(None)):
      enable_positions = np.where(np.transpose(self.matrix_mask) == 1)
      if (len(np.transpose(enable_positions)) != len(matrix)):
        raise Exception("The input matrix can not match the matrix_mask!")
      masked_matrix = self.matrix_mask.copy().astype(np.double)
      for i,position in enumerate(np.transpose(enable_positions)):
        masked_matrix[position[1], position[0]] = matrix[i]
    elif (len(matrix.shape) != 2):
      raise Exception("The input matrix should be two-dimensional when matrix_mask not specified!")
    else:
      masked_matrix = matrix
       
    if (len(masked_matrix.shape) != 2):
      raise Exception("The input matrix should be two-dimensional!")
    if (type(self.__lastest_array) == type(None)):
        self.__lastest_array = np.array(masked_matrix,dtype=np.double)
        self.__last_array = np.array(masked_matrix,dtype=np.double)
    else:
        self.__lastest_array = np.array(masked_matrix,dtype=np.double)
        self.__diff = self.__lastest_array - self.__last_array
        self.__last_array = np.array(masked_matrix,dtype=np.double)
        reconfig_positions = np.where(~np.isclose(np.abs(self.__diff), 0))
        command = ""
        for position in np.transpose(reconfig_positions):
          x_length = self.pixel_x_length * self.__lastest_array[position[0], position[1]]
          y_length = self.pixel_y_length * self.__lastest_array[position[0], position[1]]
          disable_flag = 0
          if (x_length < 0.001):
            x_length = 0
            disable_flag = 1
          if (y_length < 0.001):
            y_length = 0
            disable_flag = 1
          if (np.isclose(x_length, self.pixel_x_length) or x_length > self.pixel_x_length):
            x_length = self.pixel_x_length
          if (np.isclose(y_length, self.pixel_y_length) or y_length > self.pixel_y_length):
            y_length = self.pixel_y_length

          command += 'select("{}");'.format(self.group_name + str(position[0]) + "_" + str(position[1]))
          command += 'set("x span", {:.6f}e-6);'.format(x_length)
          command += 'set("y span", {:.6f}e-6);'.format(y_length)
          if (disable_flag):
            command += 'set("enabled", 0);'
          else:
            command += 'set("enabled", 1);'

        command_list = command.split(";")[:-1]
        block_length = int(len(command_list) / 10000) + 1
        for i in range(0, block_length):
          command_block = ";".join(command_list[i * 10000: (i + 1) * 10000])
          if (command_block != ""):
            command_block += ";"
            time.sleep(self.relaxing_time)
            #self.fdtd_engine.fdtd.eval(command_block)
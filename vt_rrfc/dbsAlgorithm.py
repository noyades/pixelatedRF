# The following are adaptations of algorithms from https://github.com/Hideousmon/SPLayout. Modifications are specific to implementations for RF PCBs and integrated circuits and to operate on flattened pixel arrays generated by pixelatedRF base code. Additionally, modifications were made in DBS to allow more than 1 flip per trial. 

import numpy as np
import math
class pixOps:

  def __init__(self, pixelSize, matrix_mask = None, relaxing_time = 0):
    self.__last_array = None
    self.__lastest_array = None
    self.relaxing_time = relaxing_time
    if (type(matrix_mask) != type(None)):
      self.matrix_mask = np.array(matrix_mask, dtype=np.int32)
    else:
      self.matrix_mask = matrix_mask

  def updatePixels(self, matrix):
    if (type(self.matrix_mask) != type(None)):
      enable_positions = np.where(self.matrix_mask == 1)
      if (len(enable_positions) != len(matrix)):
        raise Exception("The input matrix cannot match the matrix mask!")
      masked_matrix = self.matrix_mask.copy().astype(np.double)
      for i,position in enumerate(np.transpose(enable_positions)):
        masked_matrix[position[1],position[0]] = matrix[i]
    #elif (np.size(matrix) != len():
    #  print(np.shape(matrix)[:,0])
    #  raise Exception("The input matrix should be one-dimensional when matrix_mask not specified!")
    else:
      masked_matrix = matrix
    if (type(self.__lastest_array) == type(None)):
      self.__lastest_array = np.array(masked_matrix,dtype=np.double)
      self._lastest_array = np.array(masked_matrix,dtype=np.double)
      newPixMap = matrix
      return newPixMap
    else:
      self.__lastest_array = np.array(masked_matrix,dtype=np.double)
      self.__diff = self.__lastest_array - self.__last_array
      self.__last_array = np.array(masked_matrix,dtype=np.double)
      reconfig_positions = np.where(abs(self.__diff) != 0)
      newPixMap = np.zeros(len(matrix),dtype=int)
      newPixMap[reconfig_positions] == 1
      return newPixMap

  #def expandPixels(self, matrix):
    

class dbsAlgo:

  def __init__(self,loS,cost_function,rows,symmetry,minmax,simulPositions = 5,max_iteration = 3,callback_function = None,initial_solution = None):
    self.loS = loS
    self.cost_function = cost_function
    self.rows = rows
    self.sym = symmetry
    self.max_iteration = max_iteration
    self.sim_positions = simulPositions

    if (type(initial_solution) != type(None)):
      self.__Sol = initial_solution
    else:
      self.__Sol = np.random.randint(0, 2, size=(self.loS))

    self.cg_curve = np.zeros(max_iteration*self.loS)
    self.cost = np.zeros(1)
    self.__iter = 0
    self.best_solution = np.zeros(loS)
    self.__undisturbed = np.array(range(0,self.loS))
    self.minmax = minmax

    if (callback_function == None):
      def call_back():
        pass
      self.call_back = call_back
    else:
      self.call_back = callback_function

    self.__engine_init()

  def __engine_init(self):
    self.cost = self.cost_function(self.__Sol)
    self.best_solution = self.__Sol.copy()
    self.__iter = 0

  def run(self):
    """
    Run the DBS engine.
    """
    while (self.__iter < self.max_iteration):
      self.__undisturbed = np.array(range(0, self.loS))
      for i in range(0,self.loS,self.sim_positions):
        if self.__undisturbed.size <= 1:
          #self.__iter += 1
          break
        else:
          temp_solution = self.__Sol.copy()
          print('sim=' + str(i) + ' Positions remaining=' + str(self.__undisturbed.size))
          if self.sym == 'x-axis':
          # for structures that require symmetry about the x-axis (e.g., vertical symmetry), this routine will find and flip 
          # a mirror symmetric pixel in that plane.
            cols = int(self.loS / self.rows)
            print('Pixels=' + str(self.loS) + ' Rows=' + str(self.rows) + ' Cols=' + str(cols))
            temp_solution = temp_solution.reshape(self.rows,cols)
            self.__undisturbed = np.array(range(0,int(np.ceil(self.rows/2)*cols))) #only need half the pixels 
            for j in range(0,self.sim_positions):
              perturbate_shuffle = np.random.randint(0,self.__undisturbed.size) 
              perturbate_position = self.__undisturbed[perturbate_shuffle] 
              temp_solution[int(perturbate_position/cols),perturbate_position%cols] = \
                   int((temp_solution[int(perturbate_position/cols),perturbate_position%cols] + 1)%2) 
              temp_solution[self.rows - int(perturbate_position/cols)-1,perturbate_position%cols] = \
                   int((temp_solution[self.rows - int(perturbate_position/cols)-1,perturbate_position%cols] + 1)%2) 
              perturbate_position = self.loS - self.__undisturbed[perturbate_shuffle]
            temp_solution = temp_solution.reshape(self.loS,1)
            new_cost = self.cost_function(temp_solution)
            print('sim=' + str(i) + ' Simult. Positions' + str(j) + ' Positions remaining=' + str(self.__undisturbed.size))
            if self.minmax == 'min':
              if (new_cost <= self.cost):
                self.__Sol = temp_solution
                self.cost = new_cost
                self.best_solution = self.__Sol
            elif self.minmax == 'max':
              if (new_cost >= self.cost):
                self.__Sol = temp_solution
                self.cost = new_cost
                self.best_solution = self.__Sol
          else:
            for j in range(0,self.sim_positions):
              perturbate_shuffle = np.random.randint(0,self.__undisturbed.size)
              perturbate_position = self.__undisturbed[perturbate_shuffle]
              temp_solution[:,perturbate_position] = int((temp_solution[:,perturbate_position] + 1)%2)
              if (i+j != self.loS -1):
                self.__undisturbed = np.delete(self.__undisturbed,perturbate_shuffle)
            new_cost = self.cost_function(temp_solution)
            print('sim=' + str(i) + ' Simult. Positions' + str(j) + ' Positions remaining=' + str(self.__undisturbed.size))
            if self.minmax == 'min':
              if (new_cost <= self.cost):
                self.__Sol = temp_solution
                self.cost = new_cost
                self.best_solution = self.__Sol
            if self.minmax == 'max':
              if (new_cost >= self.cost):
                self.__Sol = temp_solution
                self.cost = new_cost
                self.best_solution = self.__Sol

          print('sim=' + str(int((self.__iter * self.loS + i)/self.sim_positions)))
          self.cg_curve[int((self.__iter * self.loS + i)/self.sim_positions)] = self.cost
          self.call_back()
      
      self.__iter += 1

  def get_remained_size(self):
    """
    Get the size of undisturbed positions.
    Returns
    -------
    out : Int
      Size of undisturbed positions.
    """
    return len(self.__undisturbed)

  def get_remained(self):
    """
    Get the undisturbed positions.
    Returns
    -------
    out : Array
      Undisturbed positions.
    """
    return self.__undisturbed

  def get_iteration_number(self):
    """
    Get the temporal iteration number.
    Returns
    -------
    out : Int
      Iteration number.
    """
    return self.__iter

  def get_cost(self):
    """
    Get the temporal cost.
    Returns
    -------
    out : Float
      cost.
    """
    return self.cost

  def get_best_solution(self):
    """
    Get the temporal best solution.
    Returns
    -------
    out : Array
      Best solution.
    """
    return self.best_solution

class bpsAlgo:
    

    def __init__(self, 
                 noS, 
                 loS, 
                 cost_function, 
                 minmax, 
                 max_iteration = 500, 
                 callback_function=None, 
                 initial_solution = None, 
                 v_max = 6, 
                 inertia_weight = 0.99, 
                 c_1 = 2, 
                 c_2 = 2, 
                 ratio_personal = 0.2, 
                 ratio_global = 0.8):
        """
        Binary Particle Swarm Optimization Algorithm.

        Parameters
        ----------
        noS : Int
            Initial Number of solutions.
        loS : Int
            Length of a single solution.
        cost_function : func
            Cost function for evaluating a single solution, input: Array, size (loS,), output: Float, lower means better .
        minmax : string
            Direction of optimization, to minimize FoM 'min', to maximize 'max'
        max_iteration : Int
            Maximum of iterations (default: 500).
        callback_function : func
            Self-defined callback function that will be called after every iteration (default: None).
        v_max : Float or Int
            Maximum of particle velocity (default: 6).
        inertia_weight : Float
            Intertia weight for particles (default: 0.99).
        c_1 : Float or Int
            Learning rate for self-cognition (default: 2).
        c_2 : Float or Int
            Learning rate for social-cognition (default: 2).
        ratio_personal : Float
            Ratio for self-cognition (default: 0.2).
        ratio_global : Float
            Ratio for social-cognition (default: 0.8).
        """
        self.noS = noS
        self.loS = loS
        self.cost_function = cost_function
        self.max_iteration = max_iteration
        self.callback_function = callback_function
        self.v_max = v_max
        self.inertia_weight = inertia_weight
        self.c_1 = c_1
        self.c_2 = c_2
        self.ratio_personal = ratio_personal
        self.ratio_global = ratio_global
        self.minmax = minmax
        
        #self.__Sol = np.random.randint(0,2,size=(noS,loS)) # Initialize the solutions
        if (type(initial_solution) != type(None)):
          self.__Sol = np.concatenate((initial_solution, np.random.randint(0, 2, size=(noS-1,loS))),0)
        else:
          self.__Sol = np.random.randint(0, 2, size=(noS,loS))
        self.__Best_Sol = self.__Sol.copy()
        self.__v = np.zeros((noS,loS))
        self.cg_curve = np.zeros((max_iteration))
        self.__cost = np.zeros(noS) ## the cost of the population, the lower , the better (1 - FoM)
        self.__engine_flag = 0
        self.__iter = 0
        self.best_solution = np.zeros((1,loS))
        self.max_cost = -math.inf
        self.min_cost = math.inf

        self.engine_init()
        
        if (callback_function == None):
            def call_back():
                pass
            self.call_back = call_back
        else:
            self.call_back = callback_function

    def engine_init(self):
        """
        Initialize the Binary Particle Swarm Optimization, evaluate the first iteration.
        """
        for i in range(0, self.__Sol.shape[0]):
            self.__cost[i] = self.cost_function(self.__Sol[i, :])
        self.min_cost = np.min(self.__cost, axis=0)
        self.max_cost = np.max(self.__cost, axis=0)
        __min_position = np.argmin(self.__cost, axis=0)
        __max_position = np.argmax(self.__cost, axis=0)
        if self.minmax == 'min':
          self.best_solution = self.__Sol[__min_position, :].copy()
        elif self.minmax == 'max':
          self.best_solution = self.__Sol[__max_position, :].copy()
        ## Initialize the iteration
        self.__iter = 0
        self.__engine_flag = 1


    def run(self):
        """
        Run the engine.
        """
        if(self.__engine_flag == 0):
            raise Exception("Engine has not been initialized, run: \"obj.engine_init()\" first")
        while (self.__iter < self.max_iteration):
            if self.minmax == 'min':
                self.cg_curve[self.__iter] = self.min_cost
            if self.minmax == 'max':
                self.cg_curve[self.__iter] = self.max_cost
            self.__iter += 1

            for i in range(0, self.noS):
                self.__v[i,:] = self.inertia_weight*self.__v[i,:] + self.c_1*self.ratio_personal*(self.__Best_Sol[i,:] - self.__Sol[i,:]) + \
                    self.c_2*self.ratio_global*(self.best_solution - self.__Sol[i,:])

                self.__v[i,:] = np.clip(self.__v[i,:], -self.v_max, self.v_max)
                mapped_v =  1/(1+(np.exp((-self.__v[i,:]))))

                self.__Sol[i,:] = np.random.rand(self.loS) <= mapped_v

                ## Calculate the cost
                new_cost = self.cost_function(self.__Sol[i,:])
                if self.minmax == 'min':
                  if (new_cost <= self.__cost[i]) :
                      self.__Best_Sol[i, :] = self.__Sol[i,:].copy()
                      self.__cost[i] = new_cost
                  if new_cost <= self.min_cost:
                      self.best_solution = self.__Sol[i,:].copy()
                      self.min_cost = new_cost
                elif self.minmax == 'max':
                  if new_cost >= self.max_cost:
                      self.__Best_Sol[i, :] = self.__Sol[i,:].copy()
                      self.__cost[i] = new_cost
                  if new_cost >= self.max_cost:
                      self.best_solution = self.__Sol[i,:].copy()
                      self.max_cost = new_cost

            ## Call back function
            self.call_back()

        self.__engine_flag = 0

    def get_iteration_number(self):
        """
        Get the temporal iteration number.
        Returns
        -------
        out : Int
            Iteration number.
        """
        return self.__iter

    def get_min_cost(self):
        """
        Get the temporal minimum of cost.
        Returns
        -------
        out : Float
            Minimum of cost.
        """
        return self.min_cost

    def get_max_cost(self):
        """
        Get the temporal minimum of cost.
        Returns
        -------
        out : Float
            Minimum of cost.
        """
        return self.max_cost

    def get_best_solution(self):
        """
        Get the temporal best solution.
        Returns
        -------
        out : Array
            Best solution.
        """
        return self.best_solution

    def get_total_solutions(self):
        """
        Get the temporal total solutions.
        Returns
        -------
        out : Array
            All the solutions, size: (noS,loS).
        """
        return self.__Sol

    def get_total_cost(self):
        """
        Get the temporal cost for all the solutions.
        Returns
        -------
        out : Array
            cost, size: (noS,1).
        """
        return self.__cost

class bbAlgo:
    """
    Binary Bat Algorithm.
    Parameters
    ----------
    noS : Int
        Number of solutions.
    loS : Int
        Length of a single solution.
    cost_function : func
        Cost function for evaluating a single solution, input: Array, size (loS,), output: Float, lower means better .
    max_iteration : Int
        Maximum of iterations (default: 500).
    callback_function : func
        Self-defined callback function that will be called after every iteration (default: None).
    loudness : Float
        Loudness in Binary Bat Algorithm (default: 0.25).
    pulse_rate : Float
        Pulse rate in Binary Bat Algorithm (default: 0.1).
    """
    def __init__(self, noS , loS, cost_function , minmax, max_iteration = 500,callback_function=None, initial_solution = None,loudness = 0.25, pulse_rate = 0.1):
        self.max_iteration = max_iteration
        self.noS = noS
        self.loS = loS
        self.minmax = minmax
        self.loudness = loudness
        self.pulse_rate = pulse_rate
        self.cost_function = cost_function
        ## some default parameters
        self.__Qmin = 0
        self.__Qmax = 2
        self.__N_iter = 0
        ## initial arrays
        self.__Q = np.zeros((noS,1)) # Frequency
        self.__v = np.zeros((noS,loS)) # Velocities
        #self.__Sol = np.random.randint(0,2,size=(noS,loS)) # Initialize the solutions
        if (type(initial_solution) != type(None)):
          self.__Sol = np.concatenate((initial_solution, np.random.randint(0, 2, size=(noS-1,loS))),0)
        else:
          self.__Sol = np.random.randint(0, 2, size=(noS,loS))
        self.cg_curve = np.zeros((max_iteration))
        self.__cost = np.zeros((noS,1)) ## the cost of the population, the lower , the better (1 - FoM)
        self.__engine_flag = 0
        self.__iter = 0
        self.best_solution = np.zeros((1,loS))
        self.min_cost = math.inf
        self.max_cost = -math.inf
        self.__min_position = math.inf
        self.__max_position = -math.inf
        self.engine_init()
        if (callback_function == None):
            def call_back():
                pass
            self.call_back = call_back
        else:
            self.call_back = callback_function

    def engine_init(self):
        """
        Initialize the Binary Bat Algorithm, evaluate the first iteration.
        """
        for i in range(0, self.noS):
            self.__cost[i] = self.cost_function(self.__Sol[i, :])
        self.min_cost = np.min(self.__cost, axis=0)[0]
        self.max_cost = np.max(self.__cost, axis=0)[0]
        self.__min_position = np.argmin(self.__cost, axis=0)[0]
        self.__max_position = np.argmax(self.__cost, axis=0)[0]
        if self.minmax == 'min':
            self.best_solution = self.__Sol[self.__min_position, :].copy()
        if self.minmax == 'max':
            self.best_solution = self.__Sol[self.__max_position, :].copy()

        ## Initialize the iteration
        self.__iter = 0
        self.__engine_flag = 1

    def run(self):
        """
        Run the engine.
        """
        if(self.__engine_flag == 0):
            raise Exception("Engine has not been initialized, run: \"obj.engine_init()\" first")
        while (self.__iter < self.max_iteration):
            if self.minmax == 'min':
                self.cg_curve[self.__iter] = self.min_cost
            elif self.minmax == 'max':
                self.cg_curve[self.__iter] = self.max_cost
            self.__iter += 1
            for i in range(0,self.noS):
                ## create a temporal solution
                temp_solution = self.__Sol[i,:]
                for j in range(0,self.loS):
                    self.__Q[i] = self.__Qmin + (self.__Qmin - self.__Qmax)*np.random.rand() # Equation 3
                    self.__v[i,j] = self.__v[i,j] + (temp_solution[j] - self.best_solution[j]) * self.__Q[i] # Equation 1

                    V_shaped_transfer_function = abs((2/math.pi)*math.atan((math.pi/2)*self.__v[i,j]))

                    if np.random.rand() < V_shaped_transfer_function :
                        temp_solution[j] = (temp_solution[j] + 1 ) %2

                    if np.random.rand() > self.pulse_rate:
                        temp_solution[j] = self.best_solution[j].copy()

                ## Calculate the cost
                new_cost =  self.cost_function(temp_solution)
                if self.minmax == 'min':
                  if (new_cost <= self.__cost[i]) and (np.random.rand() < self.loudness):
                      self.__Sol[i,:] = temp_solution
                      self.__cost[i] = new_cost
  
                  # Ppdate the current best
                  if new_cost <= self.min_cost:
                      self.best_solution = temp_solution.copy()
                      self.min_cost = new_cost
                if self.minmax == 'max':
                  if (new_cost >= self.__cost[i]) and (np.random.rand() < self.loudness):
                      self.__Sol[i,:] = temp_solution
                      self.__cost[i] = new_cost
  
                  # Ppdate the current best
                  if new_cost >= self.max_cost:
                      self.best_solution = temp_solution.copy()
                      self.max_cost = new_cost

            ## Call back function
            self.call_back()

        self.__engine_flag = 0

    def get_iteration_number(self):
        """
        Get the temporal iteration number.
        Returns
        -------
        out : Int
            Iteration number.
        """
        return self.__iter

    def get_min_cost(self):
        """
        Get the temporal minimum of cost.
        Returns
        -------
        out : Float
            Minimum of cost.
        """
        return self.min_cost

    def get_max_cost(self):
        """
        Get the temporal minimum of cost.
        Returns
        -------
        out : Float
            Minimum of cost.
        """
        return self.max_cost

    def get_best_solution(self):
        """
        Get the temporal best solution.
        Returns
        -------
        out : Array
            Best solution.
        """
        return self.best_solution

    def get_total_solutions(self):
        """
        Get the temporal total solutions.
        Returns
        -------
        out : Array
            All the solutions, size: (noS,loS).
        """
        return self.__Sol

    def get_total_cost(self):
        """
        Get the temporal cost for all the solutions.
        Returns
        -------
        out : Array
            cost, size: (noS,1).
        """
        return self.__cost
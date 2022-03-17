import math
from math import sqrt
from math import e as exp
import seaborn as sns
import statsmodels.api as sm
import random
from scipy import optimize
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter, median_filter

class River:
    def __init__(self):
        self.error = 0
        ##### logit model for probability of amplification
        probfunction = pd.DataFrame([[-2,.333],[-1, .875],[0,1],[1,1], [-10,0], [-3,0]], columns=['initial eDNA', 'probability of amplification'])
        probfunction['copy number'] = probfunction['initial eDNA'].apply(lambda x: 10**x * 3.65*1e5)
        model2 = sm.Logit(probfunction['probability of amplification'].values, probfunction['copy number'].values)
        self.result2 = model2.fit( disp=0)
        self.PofCaptureNet = 0.01

    def change_er(self,er):
        self.error = er    
        
    def init_river_params(self, V,D,u,λ, H):
        self.V = V
        self.u = u
        self.λ = λ
        self.D = D
        self.H = H
        if { 'V' , 'u' , 'λ' , 'D' , 'BV' , 'T' , 'pf' , 'B' } <= self.__dict__.keys():
            self.inf, self.sup = self.Find_detection_range(0.001)
    def init_sampling_strategy(self, pf, boat_V, time):
        self.pf = pf
        self.BV = boat_V
        self.T = time
        if { 'V' , 'u' , 'λ' , 'D' , 'BV' , 'T' , 'pf' , 'B' } <= self.__dict__.keys():
            self.inf, self.sup = self.Find_detection_range(0.001) #0.005
        
    def init_fish(self, dist_bet_fish, biomass):
        self.B = biomass
        self.dist = dist_bet_fish
        if { 'V' , 'u' , 'λ' , 'D' , 'BV' , 'T' , 'pf' , 'B' } <= self.__dict__.keys():
            self.inf, self.sup = self.Find_detection_range(0.001) #0.005
        
    def CtoP(self,c):
        return (self.result2.predict(c)-.5)/.5
        #return (1/(1+np.exp(-(-.83+ .00781*c)))).reshape(-1,1)

        
    def CeDNA_1_Source(self,x):
        constants = (self.B*self.u)/(sqrt(self.V**2 + 4*self.D*self.λ)*self.H)
        if x < 0: 
            result = constants * exp**(  (self.V+ sqrt(self.V**2 + 4*self.D*self.λ))*x / (2*self.D)  )
        else: 
            result = constants * exp**(  (self.V- sqrt(self.V**2 + 4*self.D*self.λ))*x / (2*self.D)  )
        if result < 1: return 0
        return result
    
    def fish_locations_transect(self):
        x = 0
        returnable = []
        while x > self.inf :
            a = -np.random.exponential(self.dist, 1)[0]
            x += a
            returnable += [x]
        returnable = returnable[::-1]
        x = 0        
        while x <  self.BV*self.T + self.sup:
            a = np.random.exponential(self.dist, 1)[0]
            x += a
            returnable += [x]
        return returnable

    def _fish_locations_net(self):
        ret = []
        x = 0
        while x< self.BV*self.T:
            x += np.random.exponential(self.dist, 1)[0]
            ret += [x]
        return ret[:-1]
            
    def average_catch(self, n = 1000):
        lis = np.array([sum([1 if random.random()< self.PofCaptureNet else 0 for x in self._fish_locations_net()]) for x in range(n)])
        return {'mean': lis.mean(), 'std': lis.std(), 'Prob_of_detection': 1 - (np.count_nonzero(lis)/len(lis)), 'list': lis}
    
    #@staticmethod
    def _solved_river_abv(self,x):
        return  -(self.pf/(self.BV*self.H)) *(2*self.B*self.u*self.D) /( 4*self.D*self.λ - self.V*sqrt(self.V**2 + 4*self.D * self.λ)+ self.V**2)* exp**( (self.V - sqrt(self.V**2 + 4*self.D*self.λ))/ (2*self.D) * x )

    def _solved_river_bl(self,x):
        return  (self.pf/(self.BV*self.H)) *(2*self.B*self.u*self.D) /(4*self.D*self.λ + self.V*sqrt(self.V**2 + 4*self.D * self.λ)+ self.V**2) * exp**( (self.V + sqrt(self.V**2 + 4*self.D*self.λ))/ (2*self.D) * x )

    def _sld_intermediary(self,Xi, Xf):
        low, high = sorted([Xi, Xf])
        if low >= 0:
            return abs(self._solved_river_abv(Xf) - self._solved_river_abv(Xi))
        if high <= 0: 
            return abs(self._solved_river_bl(Xf) - self._solved_river_bl(Xi))

        return self._sld_intermediary(low, 0) + self._sld_intermediary(0, high)
    
    def sample_eDNA_transect(self,x0):
        ret = self._sld_intermediary(x0, x0 + self.BV*self.T)  + random.gauss(0, self.error)
        if ret< 0: return 0
        else: return ret
    
    def sample_eDNA_transect_n_sources(self):
        return np.array([self.sample_eDNA_transect( -dis )*(1+random.gauss(0, self.error))  for dis in self.fish_locations_transect()]).sum()
    
    
    def Sample_Multiple_Transects_With_Different_Fish_Distances(self, dist_range =  [0,100], n = 1000):
        store_dist = self.dist
        response = []
        if len(dist_range) == 2:
            distlist = [random.uniform(dist_range[0], dist_range[1]) for i in range(n)]
            
        else: distlist = dist_range
    
        for i in distlist:
            self.dist = i
            response += [self.sample_eDNA_transect_n_sources()]
        
        self.dist = store_dist
        response = self.CtoP(response)
        return {'distances': distlist, 'response': response, 'avg': np.array(response).mean(), 'std': np.array(response).std()}
    
    
    def Catch_Transects_With_Different_Fish_Distances(self, dist_range =  [0,100], n = 1000):
        store_dist = self.dist
        response = []
        distlist = [random.uniform(dist_range[0], dist_range[1]) for i in range(n)]

        if len(dist_range) == 2:
            distlist = [random.uniform(dist_range[0], dist_range[1]) for i in range(n)]
            
        else: distlist = dist_range
        
        for i in distlist:
            self.dist = i
            response += [self.average_catch(n=1)['mean']]
        
        self.dist = store_dist
        response = np.array(response)
        det_dist = np.array([1 if x> 0 else 0 for x in response])
        
        return {'distances': distlist, 'catch': response, 'detection': det_dist,'avg': response.mean(), 'std': response.std(), 'avg_detection': det_dist.mean(), 'std_detection': det_dist.std()}

    def Find_detection_range(self, p):
        max_up = optimize.bisect(lambda d: self.CtoP(self._sld_intermediary(-d, -d+self.BV*self.T))[0] - p, 0, 1e10)
        max_down = optimize.bisect(lambda d: self.CtoP(self._sld_intermediary(-d, -d+self.BV*self.T))[0] - p, -1e10, 0)
        return sorted([max_down, max_up])
    def print_params(self):
        print(' '.join([i+'='+str(self.__dict__[i]) for i in list(self.__dict__.keys())[3:-2]]))

    
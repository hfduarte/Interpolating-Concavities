import numpy as np
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
import shapely.wkt
from shapely.wkt import dumps, loads
import math
import sys

"""
    
    A Moving segment

"""

class MSeg2:
    
    def __init__(self, i_ep_1, f_ep_1, i_ep_2, f_ep_2, starting_time = 0, ending_time = 1, level = 1, final = False):
        self.i_ep_1 = i_ep_1
        self.f_ep_1 = f_ep_1
        self.i_ep_2 = i_ep_2
        self.f_ep_2 = f_ep_2
        
        #self.starting_time = float(starting_time)
        #self.ending_time = float(ending_time)
        #self.final = final
        
        self.starting_time = starting_time
        self.ending_time = ending_time
        self.level = level
        self.final = final

    def at2(self, t, start_time, step):
        if self.starting_time == None:
            if self.ending_time == None:
                self.starting_time = float(start_time + step * self.level)
                self.ending_time = float(self.starting_time + step)
                if self.ending_time > 1:
                    self.ending_time = 1
            else:
                self.starting_time = float(start_time + (step * self.level) + step)
        
        #print(start_time, step)
        #print(self.starting_time, self.ending_time)
        
        if self.starting_time < 0 or self.ending_time > 1:
            return None, None, None, None
        
        #if t < self.starting_time:
        if t <= self.starting_time:
            return None, None, None, None
                
        if t > self.ending_time:
            if self.final:
                return self.i_ep_2[0], self.i_ep_2[1], self.f_ep_2[0], self.f_ep_2[1]
            else:
                return None, None, None, None
            
        i_x = 0
        i_y = 0
        f_x = 0
        f_y = 0
            
        if self.starting_time == 0 and self.ending_time == 1:
            i_x = self.i_ep_1[0] + t * (self.i_ep_2[0] - self.i_ep_1[0])
            i_y = self.i_ep_1[1] + t * (self.i_ep_2[1] - self.i_ep_1[1])

            f_x = self.f_ep_1[0] + t * (self.f_ep_2[0] - self.f_ep_1[0])
            f_y = self.f_ep_1[1] + t * (self.f_ep_2[1] - self.f_ep_1[1])
        else:
            if self.starting_time != 0 and self.ending_time == 1:
                w1_n = (1 - t)
                w1_d = (1 - self.starting_time)
                    
                dt = float(w1_n) / w1_d
                t = float(t - self.starting_time) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])
            elif self.starting_time == 0 and self.ending_time != 1:
                w1_n = (self.ending_time - t)
                w1_d = self.ending_time

                dt = float(w1_n) / w1_d
                t = float(t) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])
            else:
                w1_n = (self.ending_time - t)
                w1_d = self.ending_time - self.starting_time

                dt = float(w1_n) / w1_d
                t = float(t - self.starting_time) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])

        return i_x, i_y, f_x, f_y

    def at(self, t):
        if self.starting_time < 0 or self.ending_time > 1:
            return None, None, None, None
        
        #if t < self.starting_time:
        if t <= self.starting_time:
            return None, None, None, None
                
        if t > self.ending_time:
            if self.final:
                return self.i_ep_2[0], self.i_ep_2[1], self.f_ep_2[0], self.f_ep_2[1]
            else:
                return None, None, None, None
            
        i_x = 0
        i_y = 0
        f_x = 0
        f_y = 0
            
        if self.starting_time == 0 and self.ending_time == 1:
            i_x = self.i_ep_1[0] + t * (self.i_ep_2[0] - self.i_ep_1[0])
            i_y = self.i_ep_1[1] + t * (self.i_ep_2[1] - self.i_ep_1[1])

            f_x = self.f_ep_1[0] + t * (self.f_ep_2[0] - self.f_ep_1[0])
            f_y = self.f_ep_1[1] + t * (self.f_ep_2[1] - self.f_ep_1[1])
        else:
            if self.starting_time != 0 and self.ending_time == 1:
                w1_n = (1 - t)
                w1_d = (1 - self.starting_time)
                    
                dt = float(w1_n) / w1_d
                t = float(t - self.starting_time) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])
            elif self.starting_time == 0 and self.ending_time != 1:
                w1_n = (self.ending_time - t)
                w1_d = self.ending_time

                dt = float(w1_n) / w1_d
                t = float(t) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])
            else:
                w1_n = (self.ending_time - t)
                w1_d = self.ending_time - self.starting_time

                dt = float(w1_n) / w1_d
                t = float(t - self.starting_time) / w1_d

                i_x = dt * self.i_ep_1[0] + t * (self.i_ep_2[0])
                i_y = dt * self.i_ep_1[1] + t * (self.i_ep_2[1])

                f_x = dt * self.f_ep_1[0] + t * (self.f_ep_2[0])
                f_y = dt * self.f_ep_1[1] + t * (self.f_ep_2[1])

        return i_x, i_y, f_x, f_y
    
    def get_val(self, idx):
        if idx == 0:
            return self.i_ep_1
        elif idx == 1:
            return self.f_ep_1
        elif idx == 2:
            return self.i_ep_2
        elif idx == 3:
            return self.f_ep_2
    
    def get_time(self):
        return self.starting_time, self.ending_time
    
    def get_ipoint(self):
        return self.i_ep_1, self.f_ep_1
    
    def to_string(self):
        out = str(self.i_ep_1[0]) + ' ' + str(self.i_ep_1[1]) + ' ' + str(self.f_ep_1[0]) + ' ' + str(self.f_ep_1[1]) + ' '
        out += str(self.i_ep_2[0]) + ' ' + str(self.i_ep_2[1]) + ' ' + str(self.f_ep_2[0]) + ' ' + str(self.f_ep_2[1]) + ' '
        out += str(self.starting_time) + ' ' + str(self.ending_time) + ' ' + str(self.final)
        
        print(out)
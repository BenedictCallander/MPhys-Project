import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import requests
import csv



class galaxy:
    def __init__(self, row, header):
        self.__dict__ = dict(zip(header, row))
    def pos_3D(self):


data = list(csv.reader(open(r'C:\Users\Administrator\Desktop\MPhys Project\Galaxy_class_dev\file.csv')))
instances = [galaxy(i, data[0]) for i in data[1:]]

print(instances[1].posx)
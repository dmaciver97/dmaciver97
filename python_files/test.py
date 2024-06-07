import csv
import random
from tkinter import X
import numpy as np
import matplotlib.pyplot as plt

step_size = 1
q1_dat, q2_dat, q3_dat, q4_dat = [], [], [], []
t_list = []
with open('.\\python_files\\tmp\\data_14-May-202413_58_16.txt', 'r') as input:
    lines = input.readlines()
    for line in lines[1:-1:step_size]:
        print(line)
        items =  line.split(',')
        q1_dat.append(float(items[0]))
        q2_dat.append(float(items[1]))
        q3_dat.append(float(items[2]))
        q4_dat.append(float(items[3]))
        t_list.append(float(items[4]))
print(q1_dat)
Q = [[q1, q2, q3, q4] for q1, q2, q3, q4 in zip(q1_dat, q2_dat, q3_dat, q4_dat)]

print(Q)
test = [q.index(max(q))+1 for q in Q]
print(test)
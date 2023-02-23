import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fitnest import fit

def row_reader(row):
    for num in row:
        if len(num) == 37:
            item1, item2 =[], []
            item1 = num[0:18]
            item2 = num[19:37]
        else:
            item3, item4= [],[]
            item3 = num[0:19]
            item4 = num[20:38]

    new_row =[float(item1), float(item2), float(item3), float(item4)]
    return new_row
        
with open('11pnt6micron3pnt4amp.txt', 'r') as in_file:
    reader = csv.reader(in_file.readlines()[4:])
    
    with open('log.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerows(row_reader(row) for row in reader)

with open('log.csv', 'r') as file:
    reader = csv.reader(file.readlines()[0::2])
    freq, Px, Py =[], [], []
    for row in reader:
        freq.append(float(row[0]))
        Px.append(float(row[1]))
        Py.append(float(row[2]))

params = ['A', 'B']
lorentz = fit(params, freq[4:300], Px[4:300])
lorentz.run()


with open(r'C:\Users\Administrator\programs\python files\results\info\post_summary.csv') as csvfile:
    reader = csv.reader(csvfile)
    headings = next(reader)
    row = next(reader)
    A = float(row[0])
    B = float(row[1])

    fc = (A/B)**0.5
    D = (2*np.pi**2/(24*B))

y_model = [D/(f**2+fc**2) for f in freq]

fig  = plt.figure() 
ax = fig.add_subplot()  
ax.plot(freq, Px)
ax.plot(freq, y_model)
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()


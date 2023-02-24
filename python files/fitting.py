from scipy.optimize import curve_fit


#Simplified Lorentzian fit of power spectrum with constants A and B
def lorentz(x,A,B):
    return (1/(A+B*x**2))

def fit(x_data, y_data):
    index = -1
    error = 1e10
    for i in range(20):
        index += 1
        popt, pcov = curve_fit(lorentz, x_data[index:], y_data[index:], p0 = [10000,10])
        A, B = popt
        y_model = [1/(A+B*x**2) for x in x_data]
        new_error = sum([np.log(((y_m-y)/5)**2) for y_m, y in zip(y_model, y_data)])
        print(new_error)
        if new_error < error:
            error = new_error
            A_best =  A
            B_best = B
            print(index)
    

    print(A_best, B_best)
    return A_best, B_best


if __name__ == '__main__':
    import csv
    import numpy as np
    import matplotlib.pyplot as plt
    with open('log.csv', 'r') as file:
        reader = csv.reader(file.readlines())
        freq, Px, Py =[], [], []
        for row in reader:
            freq.append(float(row[0]))
            Px.append(float(row[1]))
            Py.append(float(row[2]))

    A, B  = fit(freq, Px)
    fc = (A/B)**0.5
    D = (2*np.pi**2/(15*B))
    Px_model = [D/(f**2+fc**2) for f in freq]

    fig  = plt.figure()
    ax1 = fig.add_subplot() 
    ax1.plot(freq[1:], Px[1:])
    ax1.plot(freq[1:], Px_model[1:], label= f'${np.round(D,3)}/(f^2+{np.round(fc,2)}^2)$')
    ax1.legend()
    ax1.set_ylabel('Power $[V^2/Hz]$')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(1,2.0e03)
    plt.show()


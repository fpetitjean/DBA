from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from DBA import performDBA

def main():
    #generating synthetic data
    n_series = 20
    length = 200

    series = list()
    padding_length=30
    indices = range(0, length-padding_length)
    main_profile_gen = np.array(list(map(lambda j: np.sin(2*np.pi*j/len(indices)),indices)))
    for i in range(0,n_series):
        n_pad_left = np.random.randint(0,padding_length)
        #adding zero at the start or at the end to shif the profile
        series_i = np.pad(main_profile_gen,(n_pad_left,padding_length-n_pad_left),mode='constant',constant_values=0)
        #randomize a bit
        series_i = list(map(lambda j:np.random.normal(0,0.03)+j,series_i))
        assert(len(series_i)==length)
        series.append(series_i)
    series = np.array(series)

    #plotting the synthetic data
    for s in series:
        plt.plot(range(0,length), s)
    plt.draw()
    plt.figure()

    #calculating average series with DBA
    average_series = performDBA(series,w=length/10)

    #plotting the average series
    plt.plot(range(0,length), average_series)
    plt.show()

if __name__== "__main__":
    main()

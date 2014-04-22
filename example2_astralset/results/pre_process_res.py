#!/usr/bin/pthon
import numpy as Np
import sys

def get_processed(hpp_fl):
    hlxfl = Np.loadtxt(hpp_fl,dtype= 'str',skiprows =1)
    for i in range(Np.shape(hlxfl)[0]):
        hlxfl[i,0] = hlxfl[i,0][0:8]
    outfl = "processed_%s"%hpp_fl
    Np.savetxt(outfl,hlxfl,fmt = '%s',delimiter = '\t',newline = '\n')
    return hlxfl

if __name__ == "__main__":
    pfl = get_processed(sys.argv[1])
    print pfl


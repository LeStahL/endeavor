import struct
import numpy as np
import matplotlib.pyplot as plt

#for i in range(-15,15, 1):
    #packed = struct.pack('@e', pow(2.,float(i)))
    #unpacked = struct.unpack('@H', packed)[0]
    #packed = struct.pack('@e', -pow(2.,float(i)))
    #unpacked2 = struct.unpack('@H', packed)[0]
    #print( bin(unpacked), bin(unpacked2),pow(2.,float(i)))

data = [19648,22160,20704,22096,21248]

for d in data:
    sign = np.mod(d, 2) == 1
    exponent = int(np.mod(np.floor(d/2), 32))
    if sign:
        print('-', 'e', bin(exponent))
    else:
        print('+', 'e', bin(exponent))

#float sign = mod(d, 2.),
        #exponent = mod(floor(d/2.), 32.),
        #x = floor(d/64.);

    #// Return full float16 number
    #//return mix(1.,-1., sign) * pow(10., exponent-15.) * mix(x,1.+x, clamp(exponent, 0., 1.));
    #return exponent;

import os
from sys import argv

if len(argv) >= 2:
    p = argv[1].rstrip('/')

    if len(argv) >= 3:
        s = [int(i)-1 for i in argv[2:]]
    else:
        s = range(0, 21)

    for i in s:
        cmd = ("cas-offinder "+p+"/targets_%d.txt G "+p+"/outs_%d.txt")%(i+1, i+1)
        print("Running " + cmd + " ...")
        os.system(cmd)

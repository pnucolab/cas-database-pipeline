#!/usr/bin/python
from multiprocessing.pool import ThreadPool
from subprocess import call

import glob
import sys
import os

scriptpath = "/data/django/rgenome_net/manage.py"

def run(cmd):
    print ("Executing {0}...".format("'" + " ".join(cmd) + "'"))
    return cmd, call(cmd)

def batch(t):
    cmds = []
    for fn in glob.glob(sys.argv[1] + "/fixtures/*_{0}_*.json".format(t)):
        fn = os.path.normpath(fn)
        cmds.append( ["pypy", scriptpath, "loaddata", fn] )

    cmds.sort()
    for cmd, rc in ThreadPool(30).imap(run, cmds):
        if rc != 0:
            print("{0} failed with exit status: {1}".format(cmd, rc))

#batch("basic")
#batch("org")
#batch("gene")
#batch("ctc")
#batch("gs")
#batch("tgt")
#batch("tt")
batch("ot")

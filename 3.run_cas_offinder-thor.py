import sys
import os

p = sys.argv[2].rstrip("/")

if sys.argv[1] == "push":
    fp = sys.argv[3].strip("/")
    targets = [p+"/targets_%d.txt"%(i+1) for i in range(21)]
    print ("Removing all files/directories on server...")
    os.system('ssh user036@chundoong0.snu.ac.kr "rm -fr *"')
    print ("Uploading Cas-OFFinder and target files...")
    os.system("scp -r chundoong/runall.py chundoong/cas-offinder-amd " + ' '.join(targets) + " " + fp + " " + "user036@chundoong0.snu.ac.kr:/home/user036/")
    print ("Running Cas-OFFinder...")
    os.system('ssh user036@chundoong0.snu.ac.kr "python runall.py"')
    print ("Done!")
elif sys.argv[1] == "pull":
    print ("Downloading out files...")
    os.system("scp -r user036@chundoong0.snu.ac.kr:/home/user036/outs_*.txt " + p + "/")
else:
    print ("Usage:    3.run_cas_offinder-thor push organism_directory fasta_directory\n"
           "       or 3.run_cas_offinder-thor pull organism_directory")

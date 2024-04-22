import os
import sys
import subprocess
import shutil
import glob

args = sys.argv

para_list = open("parameter_list_full.csv", mode="r")

line_no = 0
for lines in para_list:
    if line_no == int(args[1]):
        line = lines.rstrip("\n")
        list_all = line.split(",")
        rep = list_all[-1]
        other_list = list_all[:-1]
        list = ' '.join(other_list)
    line_no += 1

para_list.close()

ini = 0

if os.path.isfile("regi_time.txt"):
    with open("regi_time.txt") as f:
        for line in f:
            line = line.rstrip("\t")
            tmp = line.split("\t")
            if tmp[0] == args[1]:
                ini += 1

subprocess.run("rm *.out", shell=True)

for i in range(ini, int(rep)):
    if os.path.isfile("regi_state"+args[1]+".txt"):
        subprocess.run("g++ main_restart.cpp genotype.cpp population_restart.cpp -Wall -Wextra -std=c++11 -O3 -o test"+args[1]+".out", shell=True)
        subprocess.run("./test"+args[1]+".out "+"regi_state"+args[1]+".txt", shell=True)
        subprocess.run("rm regi_state*.txt", shell=True)
    else:
        subprocess.run("g++ main.cpp genotype.cpp population.cpp -Wall -Wextra -std=c++11 -O3 -o test"+args[1]+".out", shell=True)
        subprocess.run("./test"+args[1]+".out "+list+" "+str(i), shell=True)
        subprocess.run("rm regi_state*.txt", shell=True)
exit()

import os
import sys
import subprocess
import shutil
import glob

args = sys.argv

para_list = open("parameter_list_raw.csv", mode="r")

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

subprocess.call("g++ *.cpp -Wall -Wextra -std=c++11 -O3 -o test"+args[1]+".out", shell=True)

for i in range(int(rep)):
    subprocess.call("./test"+args[1]+".out "+list, shell=True)

exit()

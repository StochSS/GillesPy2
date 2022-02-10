#!/usr/bin/env python3

import subprocess
import os
import re

versions = [
]

def SortFn(ver):
    x = re.findall("\d+",ver)
    return int(x[0])*10000+int(x[1])*100+int(x[2])

if __name__=="__main__":
    cmd3 = f"git clone https://github.com/StochSS/GillesPy2"
    print(cmd3)
    output3 = subprocess.run(cmd3, shell=True, capture_output=True)

    cmd4 = f"cd GillesPy2; git tag"
    print(cmd4)
    output4 = subprocess.run(cmd4, shell=True, capture_output=True)
    versions = output4.stdout.decode("utf-8").split("\n")

    while("" in versions) :
        versions.remove("")

    versions.sort(key=SortFn)

    for ver in versions:
        print(f"{ver}: ",end='')
        cmd2 = f"cd GillesPy2; git checkout {ver}"
        output2 = subprocess.run(cmd2, shell=True, capture_output=True)
        cmd = f'export PYTHONPATH="{os.getcwd()}/GillesPy2"; python3 ./tyson_oscillator.py'

        output = subprocess.run(cmd, shell=True, capture_output=True)
        out = output.stdout.decode("utf-8").strip()
        try:
            print(float(out))
        except:
            print("error");


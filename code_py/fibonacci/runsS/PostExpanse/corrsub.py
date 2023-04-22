import os

with open("corrsub.txt") as file:
    lines = file.readlines()

for line in lines:
    cc=line.rstrip()    
    if cc == "":
        continue
    if cc[0] == "#":
        continue
    cc = cc.split()
    cc = cc[0]

    #os.system("echo fiboPost.sub " + cc)

    os.system("sbatch fiboPost.sub " + cc)



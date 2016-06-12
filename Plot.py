__date__ = '07.06.2016'
import numpy as np
import matplotlib.pyplot as plt
import sys

tm_file = open("tm_file.txt")
glob_file = open("glob_file.txt")

tm_x = []
tm_y = []
gl_x = []
gl_y = []

tline = tm_file.readline()
while tline:
    s = tline.split()
    x = int(s[0])
    y = int(s[1].replace("\n", ""))
    tm_x.append(x)
    tm_y.append(y)
    tline = tm_file.readline()


gline = glob_file.readline()
while gline:
    s = gline.split()
    x = int(s[0])
    y = int(s[1].replace("\n", ""))
    gl_x.append(x)
    gl_y.append(y)
    gline = glob_file.readline()

tm_file.close()
glob_file.close()

xmax = 40
ymax = 50

x = gl_x
y = gl_y
hist, xedges, yedges = np.histogram2d(x, y, bins=range(ymax))
plt.pcolor(xedges, yedges, hist, cmap=plt.cm.YlOrRd_r)
plt.colorbar()
plt.title("Globular helices")
plt.xlabel("Length")
plt.ylabel("Hydrophobic residues")




plt.show()
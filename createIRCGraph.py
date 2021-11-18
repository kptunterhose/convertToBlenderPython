import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')

save = True


data = []
f = open('irc.csv', 'r')
flines = f.readlines()
for line in flines:
    line = line.split(';')
    data.append((float(line[1].replace(',', '.')), float(line[2].replace(',', '.'))))

xData = []
yData = []
for dat in data:
    yData.append(dat[0])
    xData.append(dat[1])



x = xData
y = yData

fig, ax = plt.subplots()
fig.patch.set_alpha(0.)
fig.set_figwidth(10)
fig.set_figheight(5)
plt.grid(True)
plt.xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], fontsize=15)
plt.yticks(fontsize=15)
ax.set_xlabel('Reaktionskoordinate', fontsize=15)
ax.set_ylabel('Energie', fontsize=15)

#line, = ax.plot(x, y, color='r', linewidth=3)

point, = ax.plot([], [], marker='+', markersize=7000, markeredgewidth=2, color='r', linewidth=10)

plt.plot(x, y, linestyle='-', color='b', linewidth=2) # linestyle: '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
def update(num, x, y, line):
    line.set_data(x[:num], y[:num])
    line.axes.axis([min(xData)*1.1, max(xData)*1.1, min(yData)*1.1, max(yData) - min(yData)*.1])
    return line,

def updateDot(num, x, y, line):
    point.set_data(x[num], y[num])
    point.axes.axis([min(xData)*1.1, max(xData)*1.1, min(yData)*1.1, max(yData) - min(yData)*.1])
    return point,


ani = animation.FuncAnimation(fig, updateDot, len(x), fargs=[x, y, line],
                              interval=50, blit=True)


if save:
    ani.save('test.mp4', codec='png', dpi=300, fps=24, savefig_kwargs={'transparent': True})
    os.system('ffmpeg -i test.mp4 -c:v libvpx -pix_fmt yuva420p -auto-alt-ref 0 -b:v 4M -y test.webm')
plt.show()


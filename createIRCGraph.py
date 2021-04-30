import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')

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
plt.grid()
line, = ax.plot(x, y, color='r', linewidth=3)

def update(num, x, y, line):
    line.set_data(x[:num], y[:num])
    line.axes.axis([min(xData)*1.1, max(xData)*1.1, min(yData)*1.1, max(yData) - min(yData)*.1])
    return line,

ani = animation.FuncAnimation(fig, update, len(x), fargs=[x, y, line],
                              interval=50, blit=True)
ani.save('test.mp4', writer='ffmpeg', codec='h264', dpi=300)
plt.show()

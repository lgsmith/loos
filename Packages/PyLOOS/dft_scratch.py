import numpy as np
import matplotlib.pyplot as plt
rate = 30.0
t = np.arange(0, 10, 1/rate)
x = np.sin(2*np.pi*4*t) + np.sin(2*np.pi*7*t) + np.random.randn(len(t))*0.2
p = np.abs(np.fft.fft(x))**2
frqs = np.linspace(0, rate/2, len(p))
t_fft = np.fft.fftfreq(x.size, 1/rate)
# magic circle below here
E = np.zeros_like(frqs)
for bin in range(len(frqs)):
    w = np.pi * bin/len(frqs)
    k = 2 * np.sin(w)
    y1 = 0.0
    y2 = 0.0
    for sample in x:
        y2 += sample - k * y1
        y1 += k * y2

    E[bin] = y1**2 + y2**2 - k * y1 * y2

plt.plot(t_fft, E)
plt.show()
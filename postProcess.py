import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams.update({'font.size': 15})

def FFT(Fs,data):
    L = len(data)                          # 信号长度
    N = int(np.power(2,np.ceil(np.log2(L))))      # 下一个最近二次幂
    FFT_y = abs(fft(data,N))/L*2                  # N点FFT，但除以实际信号长度 L
    Fre = np.arange(int(N/2))*Fs/N          # 频率坐标
    FFT_y = FFT_y[range(int(N/2))]          # 取一半
    return Fre, FFT_y


nt = 1001
t = 100
Fs = 10
wave = np.zeros(nt)
wave1 = np.zeros(nt)
for i in range(nt):
    p = 'outputEuler0.001/u_' + '%05d'%i + '.out'
    with open(p,encoding = 'utf-8') as f:
        data = np.loadtxt(f)
        #print(data.shape)
        wave[i] = data[200, 5]

for i in range(nt):
    p = 'outputEuler0.001DX4/u_' + '%05d'%i + '.out'
    with open(p,encoding = 'utf-8') as f:
        data = np.loadtxt(f)
        #print(data.shape)
        wave1[i] = data[100, 5]

t = np.linspace(0, t, nt)
fig, ax = plt.subplots()
ax.plot(t, wave, label="dx = 2")
ax.plot(t, wave1, label="dx = 4")
ax.set_xlabel('时间/s') #设置x轴名称 x label
ax.set_ylabel('波幅/m') #设置y轴名称 y label
ax.set_title('波幅H时历图')
ax.legend()

plt.show()

'''
x = np.arange(0, 800, 2)
wave = np.zeros(x.shape)
p = 'outputEuler0.001/h_' + '%05d'%500 + '.out'
with open(p,encoding = 'utf-8') as f:
    data = np.loadtxt(f)
        #print(data.shape)
    wave = data[:, 5]

print(wave.shape)
plt.plot(x, wave)
plt.xlabel('x坐标/s')
plt.ylabel('波幅/m')
plt.title('时间t=50s')
plt.show()
'''
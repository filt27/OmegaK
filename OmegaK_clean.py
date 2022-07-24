import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import interpolate



tsweep = 166e-3  # chirp duration   [s]
bw = 700e6  # bandwidth [Hz]
fc = 24e9   #central frequency [Hz]
step_size = 0.01    # step size [m]



c = 3e8
fs = 44e3
f0 = fc - bw/2
Rs = 1

data = np.loadtxt("path-to-matrix-file.txt")


gamma = bw/tsweep
sweep_samples = len(data[0])        #broj uzoraka
num_of_steps = len(data)
aperture_length = step_size*num_of_steps   # put koji radar prelazi [m]
f = np.linspace(0, fs/2, sweep_samples//2+1)
range = c*f*tsweep/(2*bw)
t = np.linspace(0, tsweep, len(data[0]))

print(" ")
print("** Parameters **")
print("Step size: ", step_size*1e2, " cm")
print("Number of steps: ", num_of_steps)
print("Total aperture: ", aperture_length, " m")
print("Sweep samples: ", sweep_samples)
print("Min frequency: ", f0/1e9, " GHz")
print("Bandwidth: ", bw/1e6, "MHz")


range0 = 0
range1 = c*(fs/2.0)*tsweep/(2*bw)
crange0 = -(num_of_steps - 1)*aperture_length/2.0
crange1 = (num_of_steps - 1)*aperture_length/2.0
raw_extent = (range0, range1, crange0, crange1)

data_show=np.rot90(data)
plt.figure()
ax = sns.heatmap(np.abs(data_show))
plt.title("Captured FMCW signals")


# Hilbert transform #############################################################################################

def hilbert_rvp(x, fs, gamma):
    y = np.fft.fft(x, axis=-1)
    y[:,:y.shape[1]//2+1] = 0 # Zero positive frequencies
    # Residual video phase term compensation
    f = np.linspace(-fs/2, fs/2, y.shape[1])
    y *= np.exp(1j * np.pi * f**2 / gamma)
    return np.fft.ifft(y, axis=-1)

data_hilbert = hilbert_rvp(data, fs, bw/tsweep)


# Hanning window #############################################################################################

H=[]
data_hanning=[]
count = 0
while(count<sweep_samples):
    H.append(0.5+0.5*np.cos(2*np.pi*(count-sweep_samples/2)/sweep_samples))
    count+=1
for i in data_hilbert:
    data_hanning.append(i*H)


plt.figure()
ax = sns.heatmap(np.abs(data_hanning))
plt.title("After Hilbert & RVP")


data_fft1 = np.fft.fft(data_hanning, axis=1)
data_fft1 = np.fft.fftshift(data_fft1, axes=1)



# Zero padding #############################################################################################

zpad = 4096
data_zpad=[]


zeros = np.zeros(len(data_hanning[0]))
count=0
while(count<zpad):
    data_zpad.append(zeros)
    count+=1

data_zpad[int(zpad/2)-int(len(data_hanning)/2):int(zpad/2)+int(len(data_hanning)/2)] = data_hanning


# FFT #############################################################################################


data_fft = np.fft.fft(data_zpad, axis=0)
data_fft = np.fft.fftshift(data_fft, axes=0)

data_fft_show=np.rot90(data_fft)
plt.figure()
ax = sns.heatmap(np.abs(data_fft_show))
plt.title("Data after FFT")


data_compr = np.fft.ifft(data_fft, axis=1)
data_compr = np.fft.fftshift(data_compr, axes=1)


# RFM #############################################################################################

kx = np.linspace(-np.pi/step_size, np.pi/step_size, len(data_fft))
kr = np.linspace(((4*np.pi/c)*(f0-bw/2)), ((4*np.pi/c)*(f0+bw/2)), len(data_fft[0]))

phi_mf = np.zeros((len(data_fft), len(data_fft[0])))
for index_j, j in enumerate(data_fft):
    for index_i, i in enumerate(j):
        phi_mf[index_j, index_i] = Rs*(kr[index_i]**2 - kx[index_j]**2)**0.5

smf = np.e**(1j*phi_mf)

S_mf = data_fft*smf


plt.figure()
ax = sns.heatmap(np.angle(S_mf))
plt.title("Matched Filter")
#plt.show()


data_final_ifft_range = np.fft.ifft(S_mf, axis=1)
data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 1)

data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 0)
data_final_ifft_range = np.fft.ifft(data_final_ifft_range, axis=0)



data_final_ifft_range_show=np.rot90(data_final_ifft_range)
plt.figure()
ax = sns.heatmap(np.abs(data_final_ifft_range_show))
plt.title("Before Stolt")



# Stolt #############################################################################################

kstart = 973
kstop = 1009

Ky_even = np.linspace(kstart,kstop,sweep_samples)

count = 0
ky = []
S_st = []

while(count<len(kx)):
    ky.append((kr**2-kx[count]**2)**0.5)
    func = interpolate.interp1d(ky[count], S_mf[count], fill_value="extrapolate")
    S_st.append(func(Ky_even))
    count+=1



# 2nd Hanning i zeropad

H=[]
count = 1
while(count<=len(S_st[0])):
    H.append(0.5+0.5*np.cos(2*np.pi*(count-len(S_st[0])/2)/len(S_st[0])))
    count+=1

S_st_hanning = []
for i in S_st:
    S_st_hanning.append(i*H)


# IFFT #############################################################################################

zpad2 = 1000
S_st=np.pad(S_st_hanning, zpad2, mode='constant')


data_final_ifft_range = np.fft.ifft(S_st, axis=1)
data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 1)

data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 0)
data_final_ifft_range = np.fft.ifft(data_final_ifft_range, axis=0)




# Display #############################################################################################

fig = plt.figure()

final = np.asmatrix(data_final_ifft_range)
final_radar_image=np.rot90(final)


x_len = final_radar_image.shape[0]
y_len = final_radar_image.shape[1]
resize_factor_x = 0.02
resize_factor_y = 0.05


final_radar_image = final_radar_image[int(x_len/2)-90:int(x_len/2)+30, int(y_len/2)-60:int(y_len/2)+60]

ax = sns.heatmap(np.abs(final_radar_image))

plt.title("SAR image")

ax.set_axis_off()
plt.axis('off')

max_range = (2*np.pi*len(S_st_hanning[0]))/(kstop-kstart)

X = np.arange(-12, 12+1, 1.0)
Y = np.linspace(max_range+1, -max_range, int(max_range/3), dtype = int)

xtick_locations = np.linspace(0, int(final_radar_image.shape[1]), len(X))
ytick_locations = np.linspace(0, int(final_radar_image.shape[0]), len(Y))

ax.set_xticks(xtick_locations)
ax.set_xticklabels(X)

ax.set_yticks(ytick_locations)
ax.set_yticklabels(Y)

plt.show()












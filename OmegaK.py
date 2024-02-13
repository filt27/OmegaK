import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.fftpack import *
from scipy import interpolate
from scipy import signal


filename="Pastic_Glass_Metal_object_Example.txt"
data = np.loadtxt(filename)

c = 3e8         # speed of light [m/s]
fc = 24e9       # central frequency [Hz]
bw = 1300e6     # chirp bandwidth [Hz]
step_size = 0.005 # GBSAR step size [m]
num_of_steps = len(data)        # number of GBSAR steps (in stop-and-go mode)
gbsar_direction = "L"   # GBSAR direction of sensor movement [R-to the right, L-to the left]
tsweep = 155e-3  # chirp sweep time [s]
gamma = bw/tsweep       # chirp rate [Hz/s]
sweep_samples = len(data[0])        # number of frequency samples in chirp
aperture_length = step_size*num_of_steps   # aperture length [m]
f0 = fc - bw/2      # minimum chirp frequency [Hz]
fs = int(sweep_samples/tsweep)  # sampling frequency [Hz]
f = np.linspace(0, fs/2, sweep_samples//2+1)    # sampling frequency linspace
t = np.linspace(0, tsweep, sweep_samples)    # array of times in chirp period

frequency_step = fs/sweep_samples
range_resolution = frequency_step*c/(2*gamma)   # FMCW range resolution [m]



# Print radar parameters
print(" ")
print("** GBSAR parameters **")
print("Step size: ", step_size*1e2, " cm")
print("Number of steps: ", num_of_steps)
print("Total aperture: ", aperture_length, " m")
print("Number of chirp samples: ", sweep_samples)
print("Central frequency: ", fc/1e9, " GHz")
print("Bandwidth: ", bw/1e6, "MHz")


#flip data if needed
if(gbsar_direction=="L"):
        data = np.flip(data, 0)


#plot heatmap of raw data
data_show=np.rot90(data)
plt.figure()
ax = sns.heatmap(np.abs(data_show))
plt.title("Raw data")


# linear mask for low frequencies
num_of_steps_for_mask = 3
mask = np.ones(sweep_samples)
mask[:num_of_steps_for_mask] = np.linspace(0,1,num_of_steps_for_mask)
mask[-num_of_steps_for_mask:] = np.linspace(1,0,num_of_steps_for_mask)
for index, one_step in enumerate(data):
    fmcw_fft = fft(signal.detrend(one_step))
    remove_low_freq = fmcw_fft*mask
    data[index]=ifft(remove_low_freq)


# Hilbert transform, Residual Video Phase compensation and Hanning window
def hilbert_rvp(x, fs, gamma):
    y = np.fft.fft(x, axis=-1)
    y[:,:int(sweep_samples/2)+1] = 0 # Zero positive frequencies
    # Residual video phase term compensation
    f = np.linspace(-fs/2, fs/2, y.shape[1])
    y *= np.exp(1j * np.pi * f**2 / gamma)
    return np.fft.ifft(y, axis=-1)

data_hilbert = hilbert_rvp(data, fs, gamma)

data_hanning=[]

for i in data_hilbert:
    data_hanning.append(i*np.hanning(sweep_samples))



# Zero padding and FFT #############################################################################################

zpad = 256
data_zpad=[]
zeros = np.zeros(sweep_samples)
count=0
while(count<zpad):
    data_zpad.append(zeros)
    count+=1

data_zpad[int(zpad/2)-int(num_of_steps/2):int(zpad/2)+int(num_of_steps/2)] = data_hanning

data_fft = np.fft.fft(data_zpad, axis=0)
data_fft = np.fft.fftshift(data_fft, axes=0)


# Referent Function Multiplication ###################################################################################

Rs = 0.2        # Referent distance between radar and objects [m]
kx = np.linspace(-np.pi/step_size, np.pi/step_size, len(data_fft))
kr = np.linspace(((4*np.pi/c)*(fc-bw/2)), ((4*np.pi/c)*(fc+bw/2)), sweep_samples)


phi_mf = np.zeros((len(data_fft), sweep_samples))
for index_j, j in enumerate(data_fft):
    for index_i, i in enumerate(j):
        phi_mf[index_j, index_i] = -Rs*kr[index_i] + Rs*(kr[index_i]**2 - kx[index_j]**2)**0.5

smf = np.e**(1j*phi_mf)

S_mf = data_fft*smf
S_mf = data_fft


# Stolt transform #############################################################################################

Ky_even = kr
count = 0
ky = []
S_st = []

while(count<len(kx)):
    ky.append((kr**2-kx[count]**2)**0.5)        #-Rs*kr[index_i] + Rs*(kr[index_i]**2 - kx[index_j]**2)**0.5
    func = interpolate.interp1d(ky[count], S_mf[count], fill_value="extrapolate")
    S_st.append(func(Ky_even))
    count += 1


S_st_hanning = []
for i in S_st:
    S_st_hanning.append(i * np.hanning(len(S_st[0])))



# Zero padding and 2D IFFT #############################################################################################

zpad2 = 256
S_st = np.pad(S_st_hanning, zpad2, mode='constant')


data_final_ifft_range = np.fft.ifft(S_st, axis=1)
data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 1)

data_final_ifft_range = np.fft.fftshift(data_final_ifft_range, axes = 0)
data_final_ifft_range = np.fft.ifft(data_final_ifft_range, axis=0)


# Visualization #############################################################################################


final_radar_image = np.asmatrix(data_final_ifft_range)
final_radar_image = np.rot90(final_radar_image)

plt.figure()
ax = sns.heatmap(np.abs(final_radar_image), cbar = False)
plt.show()

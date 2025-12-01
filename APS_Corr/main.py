import os
import wave
import shutil
import struct
import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import sys

# Notes from Dr. Xiao meeting
# Check math on eigenvector decomposition, rank of 100x100 proximity matrix
# Max/min of eigenvalues/vectors
# Run MDS multiple times, compare answers, see if they're unique

# Todo from Dr. Hanieh
# Use match filters (I.E., transfer-functioned value using 5second wav)
# Use 5 second white noise
# More accurate accoustic propogation (
# Consider only nearest neighbors (Nearest 8? Maybe nearest 26?)


def get_pos(n):
    #0 is bottom left, 1 is directly to the right, last is topright most...
    i, j = n%width, int(n/height)
    return poss[i][j][0], poss[i][j][1]

def get_timedelay(x, y):
    global samplingrate, speedofsound
    # Distance between points x and y, so
    posx = get_pos(x)
    posy = get_pos(y)
    distance = math.dist(posx, posy) # euclidian distance in meters
    timedelay = distance / speedofsound # distance / speedofsound = timedelay
    numsamples = int(timedelay * samplingrate) # timedelay is in seconds, samplingrate is in Hz, round down
    return numsamples # We want the final timedelay in numsamples

def cross_correlate(fullsound, signal):
    lenfull, lensignal = len(fullsound), len(signal)
    # fullsound = fullsound.astype(np.int64)
    signal = signal.astype(np.int64)
    corrs = []
    full = []
    for offset in range(0, lenfull-lensignal, 1):
        # sum = np.dot(signal, fullsound[offset:offset+lensignal])
        # corrs.append(sum)
        full.append(fullsound[offset:offset+lensignal])
    full = np.asarray(full)
    corrs = np.dot(full, signal)
    # print(corrs)
    return corrs

def square_matrix(arr):
    retM = np.zeros((len(arr), len(arr)))
    for index in range(len(arr)):
        retM[index][index] = np.real(arr[index])
    return retM

def mds(dists):
    numPoints = len(dists)
    C = np.eye(len(dists)) - (np.ones((len(dists), len(dists)))/len(dists))  # centering matrix C
    Dsquared = np.zeros((len(dists), len(dists))) # piece-wise squared distance matrix
    for i in range(len(dists)):
        for j in range(len(dists)):
            Dsquared[i][j] = (dists[i][j]*dists[i][j])
    B = -0.5*np.linalg.multi_dot((C,Dsquared,C))  # B = -0.5*C*D^(2)*C

    ret = np.linalg.eig(B)
    # Get m largest eigenvalues
    theoreticalEigenValues = np.argsort(ret.eigenvalues)[-2:]  # Get largest eigenvalues/vectors
    # print(theoreticalEigenValues)
    eigenValues = ret.eigenvalues[theoreticalEigenValues]
    # And corresponding eigenvectors
    eigenVectors = ret.eigenvectors[:, theoreticalEigenValues]
    eigenVectors = eigenVectors / np.linalg.norm(eigenVectors, axis=0, ord=2)
    d2c = (np.dot(eigenVectors, np.sqrt(square_matrix(eigenValues))))  # Recalculate new points
    # d2c = np.asmatrix(d2c)

    newdists = np.zeros((numPoints, numPoints))
    for i in range(numPoints):
        for j in range(numPoints):
            newdists[i][j] = (np.linalg.norm(d2c[i] - d2c[j]))

    rmse = 0  # Calculate difference between initial data and answer to geometry
    for i in range(numPoints):
        for j in range(numPoints):
            rmse += ((dists[i][j]-newdists[i][j])**2)/(numPoints**2)

    rmse = math.sqrt(rmse)
    d2c = np.real(d2c)
    # print(d2c)
    fig4, ax4 = plt.subplots(1, 1)
    ax4.scatter(d2c[:, 0], d2c[:, 1])
    for i in range(numbernodes):
        ax4.annotate(i, (d2c[:, 0][i], d2c[:, 1][i]))
    fig4.savefig("fig4")
    return rmse

def animate_cross_correlate(fullsound, signal, corrs):
    lenfull, lensignal = len(fullsound), len(signal)
    print("Signal: ", signal)
    animfig, (animax1, animax2) = plt.subplots(2,1)
    lines = [animax1.plot([], [])[0] for _ in range(2)] + [animax2.plot([], [])[0]] #lines to animate

    def init():
        animax1.set_xlim(0, lenfull)
        animax1.set_ylim(-10000, 10000)
        animax2.set_xlim(0, lenfull)
        animax2.set_ylim(-40000000000, 40000000000)
        xdata = [i for i in range(lenfull)]
        ydata = fullsound
        lines[1].set_data(xdata, ydata)
        return lines
    def update(frame):
        xdata = [(frame + i ) % lenfull for i in range(lensignal)]
        ydata = signal
        xdata2 = [i for i in range(int(frame))]
        ydata2 = corrs[0:int(frame)]
        lines[0].set_data(xdata, ydata)
        lines[2].set_data(xdata2, ydata2)
        return lines
    ani = FuncAnimation(animfig, update, frames=np.linspace(0, lenfull-lensignal, num=100), init_func = init, blit=True)
    writergif = animation.PillowWriter(fps=30)
    ani.save("animation.gif", writer=writergif)
    exit()


def cross_correlate2(fullsound, signal):
    lenfull, lensignal = len(fullsound), len(signal)
    # fullsound = fullsound.astype(np.int64)
    signal = signal.astype(np.int64)
    corrs = []
    full = []
    for offset in range(0, lenfull-lensignal, 1):
        sum = np.dot(signal, fullsound[offset:offset+lensignal])
        corrs.append(sum)
        # full.append(fullsound[offset:offset+lensignal])
    # full = np.asarray(full)
    # corrs = np.dot(full, signal)
    # print(corrs)
    return corrs

def readdata(wav):
    global numchans, samplingrate, numsamps, sampwidth
    s = wav.readframes(wav.getnframes())
    unpstr = "<{0}h".format(wav.getnframes())
    x = list(struct.unpack(unpstr, s))
    return np.asarray(x)

def readdata32(wav):
    global numchans, samplingrate, numsamps, sampwidth
    s = wav.readframes(wav.getnframes())
    unpstr = "<{0}i".format(wav.getnframes())
    x = list(struct.unpack(unpstr, s))
    return np.asarray(x)


# For simulation purposes, they will be arranged in a grid with 3m spacing, with node 0 at bottom left (0, 0) and the last node at the topright most
width = 5
height = 5
numbernodes = width * height

# Other simulation variables
speedofsound = 343  # Units of m/s. For this case, a distance of 3 meters means a time delay of .1
poss = []
for i in range(height):
    tmppos = []
    for j in range(width):
        tmppos.append([3 * i + np.random.normal(0, 0.25), 3 * j + np.random.normal(0, 0.25)])
    poss.append(tmppos)

# Ensures a clear trial1 directory exists
if os.path.exists("./trial1/"):
    shutil.rmtree("./trial1/")
    os.mkdir("./trial1/")
else:
    os.mkdir("./trial1/")

for i in range(numbernodes):
    os.mkdir(f"./trial1/{i}/")

pos = []
for node in range(numbernodes):
    pos.append(get_pos(node))
pos = np.asarray(pos)
fig, ax = plt.subplots(1, 1)
ax.scatter(x=pos[:,0], y=pos[:,1])
for i in range(numbernodes):
    ax.annotate(i, (pos[:,0][i], pos[:,1][i]))
fig.savefig("fig")

chirp = wave.Wave_write("./chirp.wav")
chirp.setparams((1, 2, 4800, 480, 'NONE', 'not compressed'))
# numchans, sampwidth, samplingrate, numsamps, comptype, compname = chirp.getparams()
numchans, sampwidth, samplingrate, numsamps, comptype, compname = (1, 2, 4800, 4800, 'NONE', 'not compressed')
t = np.linspace(0, 1, numsamps, endpoint=False)
audio_data = ((5000 * (np.sin(2 * np.pi * 440 * t) + 0.5*np.sin(2 * np.pi * 720 * t) + 0.33*np.sin(2 * np.pi * 1050 * t)))*np.power(np.sin(np.pi*t),1.5)).astype(np.int16)
# print(audio_data)
# audio_data = np.concatenate((audio_data, audio_data, audio_data, audio_data, audio_data))
# numsamps = numsamps*5
chirp.writeframes(audio_data.tobytes())
chirp.close()
onechan = audio_data

# chirp = wave.Wave_read("./Normalized_5s_White_Noise_SR_48000_AmpCoeff_4.wav")
# numchans, sampwidth, samplingrate, numsamps, comptype, compname = chirp.getparams()
# print(chirp.getparams())
# onechan = readdata32(chirp)


# Go through and roughly simulate the time delays for each trial, and store in trial1
for speaker in range(numbernodes):
    for listener in range(numbernodes):
        filepath = f"./trial1/{speaker}/{listener}.wav"
        sampledelay = get_timedelay(speaker, listener)
        distance = math.dist(get_pos(speaker), get_pos(listener))
        tmpsamps = numsamps + samplingrate
        tmpdata = np.zeros(tmpsamps)
        if distance > 0:
            tmpdata[sampledelay:sampledelay+numsamps] = onechan/(pow(distance, 1))
        else:
            tmpdata[sampledelay:sampledelay+numsamps] = onechan
        tmpfile = wave.Wave_write(filepath)
        tmpfile.setparams((numchans, sampwidth, samplingrate, tmpsamps, 'NONE', 'not compressed'))
        if sampwidth == 2:
            tmpdata = tmpdata.astype(np.int16)
            noise = np.random.normal(loc=0, scale=100, size=tmpsamps).astype(np.int16)
        elif sampwidth == 4:
            tmpdata = tmpdata.astype(np.int32)
            noise = np.random.normal(loc=0, scale=100, size=tmpsamps).astype(np.int32)
        tmpdata2 = tmpdata + noise
        if speaker==0 and listener==8:
            # plt.plot(tmpdata)
            # plt.plot(noise)
            fig2, ax2 = plt.subplots(1, 1)
            ax2.plot(tmpdata2)
            ax2.plot(tmpdata)
            fig2.savefig("fig2")
        tmpfile.writeframes(tmpdata2.tobytes())
        tmpfile.close()

timedelays = np.zeros((numbernodes, numbernodes)) # The recalculated time of flights using crosscorrelation
timedelays2 = np.zeros((numbernodes, numbernodes)) # The ground truth
# onechantmp = onechan[0:12000]
onechantmp = onechan
# Go through and crosscorrelate everything
for speaker in range(numbernodes):
    for listener in range(numbernodes):
        filepath = f"./trial1/{speaker}/{listener}.wav"
        tmpwav = wave.open(filepath)

        # corr = signal.correlate(tmpdata, onechantmp, mode="full", method="direct") / (numsamps*10) # attempt 1
        # lags = signal.correlation_lags(len(tmpdata), len(onechantmp))
        # corr = signal.correlate(onechantmp, tmpdata, mode="full", method="direct") / (numsamps*2) # attempt 2
        # lags = signal.correlation_lags(len(onechantmp), len(tmpdata))
        # lag = lags[np.argmax(corr)]

        if sampwidth == 2:
            tmpdata = readdata(tmpwav)
            corrs = cross_correlate(tmpdata, onechantmp) # attempt 3
            if speaker == 0 and listener == 0:
                animate_cross_correlate(tmpdata, onechantmp, corrs)
        elif sampwidth == 4:
            tmpdata = readdata32(tmpwav)
            corrs = cross_correlate2(tmpdata, onechantmp) # attempt 3
        lag = np.argmax(corrs)

        timedelays[speaker][listener] = lag # put the time of arrival into timedelays
        timedelays2[speaker][listener] = get_timedelay(speaker, listener) # store ground truth

        # if speaker==0 and listener==8:
        #     fig, (ax_orig, ax_noise, ax_corr) = plt.subplots(3, 1, figsize=(4.8, 4.8))
        #     ax_orig.plot(tmpdata)
        #     ax_orig.set_title('Original signal')
        #     ax_orig.set_xlabel('Sample Number')
        #     ax_noise.plot(onechantmp)
        #     ax_noise.set_title('Signal with noise')
        #     ax_noise.set_xlabel('Sample Number')
        #     ax_corr.plot(lags, corr)
        #     ax_corr.set_title('Cross-correlated signal')
        #     ax_corr.set_xlabel('Lag')
        #     ax_orig.margins(0, 0.1)
        #     ax_noise.margins(0, 0.1)
        #     ax_corr.margins(0, 0.1)
        #     fig.tight_layout()
        if speaker==0 and listener==8:
            fig3, ax3 = plt.subplots(1,1)
            ax3.plot(corrs)
            fig3.savefig("fig3")
        tmpwav.close()

timedelays = np.asarray(timedelays)
timedelays /= samplingrate # now in units of seconds
timedelays *= speedofsound # now in units of meters
for speaker in range(numbernodes):
    timedelays[speaker] -= timedelays[speaker][speaker] # subtract the values on the diagonal from each row to get time of flight
    timedelays[speaker] = abs(timedelays[speaker])

timedelays2 = np.asarray(timedelays2)
timedelays2 /= samplingrate # now in units of seconds
timedelays2 *= speedofsound # now in units of meters
for speaker in range(numbernodes):
    timedelays2[speaker] -= timedelays2[speaker][speaker] # subtract the values on the diagonal from each row to get time of flight
    timedelays2[speaker] = abs(timedelays2[speaker])

print("Correlations RMSE: ", np.sqrt(np.mean((timedelays2-timedelays)**2)))
# print(timedelays)
rmse = mds(timedelays)
print("MDS RMSE: ", rmse)

# plt.show()
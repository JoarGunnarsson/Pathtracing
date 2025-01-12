import numpy as np
import scipy.signal
from numpy.fft import fft2, ifft2, fftshift


def gaussfft(inpic, sigma2):
	pfft = np.fft.fft2(inpic)
	[h, w] = np.shape(inpic)
	[x, y] = np.meshgrid(np.linspace(0, 1-1/w, w),np.linspace(0, 1-1/h, h))
	ffft = np.exp(sigma2 * (np.cos(2*np.pi*x) + np.cos(2*np.pi*y) - 2))
	pixels = np.real(np.fft.ifft2(ffft * pfft))
	return pixels

def medfilt(Image, wheight, wwidth = -1):
	if wwidth == -1:
		wwidth = wheight
	result = scipy.ndimage.median_filter(Image, size=(wheight, wwidth))
	return result
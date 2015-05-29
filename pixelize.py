import numpy as np
import matplotlib.pyplot as py
import math
from mpl_toolkits.mplot3d import Axes3D

pix = 0.08
slit_width = 1.0
R = 0.7*0.75

pixel = pix*pix
width_pix = slit_width/pixel

perct_last_w_pix = width_pix - math.floor(width_pix)
width_npix = math.ceil(width_pix)

length_pix = R*2./pixel
perct_last_l_pix = length_pix - math.floor(length_pix)
length_npix = math.ceil(width_pix)

pixels = np.ones((length_npix,width_npix))
pixels[:,width_npix-1] *= perct_last_l_pix
pixels[0,:] *= perct_last_w_pix


l, w = np.meshgrid(np.arange(0,length_npix,1),np.arange(0,width_npix,1))
py.pcolor(l,w,pixels)
py.grid()
py.show()




import numpy as np
def bp(img,pix,u1,u2,npad=1024):
    n=img.shape; P=np.zeros((npad,npad)); P[:n[0],:n[1]]=img
    F=np.fft.fft2(P); u=np.fft.fftfreq(npad,d=pix/206265.)/1e3; U=np.hypot(*np.meshgrid(u,u))
    t=np.clip((U-u1)/max(u2-u1,1e-9),0,1); F*=0.5*(1-np.cos(np.pi*t))
    return np.real(np.fft.ifft2(F))[:n[0],:n[1]]

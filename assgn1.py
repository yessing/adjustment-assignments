#
import numpy as np

'''
def pfpb(xi,yi,a0=1.818,b0=3.273):
    return 2*(b0-yi+xi*a0)/(a0*a0+1)

def pfpa(xi,yi,a0=1.818,b0=3.273):
    item1=2*(a0*xi-yi+b0)*xi*(a0*a0+1)
    item2=2*((a0*xi-yi+b0)**2)*a0
    return (item1-item2)/((a0*a0+1)**2)

def Dptln(xi,yi,a0,b0):
    return (a0*xi - yi + b0)**2/(a0*a0+1)




for i in range(1):
    mb=np.matrix([[pfpa(1.1,4.2,a_,b_),pfpb(1.1,4.2,a_,b_)],\
                  [pfpa(1.8,6.8,a_,b_),pfpb(1.8,6.8,a_,b_)],\
                  [pfpa(2.6,8.0,a_,b_),pfpb(2.6,8.0,a_,b_)]])

    sig0 = 1
    dx=np.matrix([[1,0,-0.5],\
                  [0,2,0],\
                  [-0.5,0,1]])

    dy=np.matrix([[1,0,0],\
                  [0,1,0],\
                  [0,0,2]])

    Lma = np.matrix([[Dptln(1.1,4.2,a_,b_)],\
                    [Dptln(1.8,6.8,a_,b_)],\
                    [Dptln(2.8,8.0,a_,b_)]])
    Lmi = mb*np.matrix([[a_],[b_]])

    px = sig0 * dx**(-1)
    py = sig0 * dy**(-1)

    PxB= px * mb
    pv=sig0*(dx+dy)**(-1)

    #xhat = (PxB.T*py*PxB)**(-1) * PxB.T * py*px*Lmi
    xhat = (mb.T*pv*mb)**(-1)*mb.T*pv*Lmi
    a_=a_+xhat[0,0]
    b_=b_+xhat[1,0]
    #print(xhat)
    print(a_,b_)


def pfpx(x0,y0,a,b):
    return 2*(a*x0-y0+b)*a/(a*a+1)

def pfpy(x0,y0,a,b):
    return -2*(a*x0-y0+b)/(a*a+1)
'''
def pfpb(x0,y0,a,b):return 1
def pfpa(x0,y0,a,b):return x0
def pfpx(x0,y0,a,b):return a
def pfpy(x0,y0,a,b):return -1
def f0(x0,y0,a,b):return a*x0-y0+b

a_ = 1.818
b_ = 3.273

for i in range(10):
    ma=np.matrix(np.zeros(shape=(3,6)))
    ma[0,0]=pfpx(1.1,4.2,a_,b_)
    ma[1,1]=pfpx(1.8,6.8,a_,b_)
    ma[2,2]=pfpx(2.8,8.0,a_,b_)

    ma[0,3]=pfpy(1.1,4.2,a_,b_)
    ma[1,4]=pfpy(1.8,6.8,a_,b_)
    ma[2,5]=pfpy(2.8,8.0,a_,b_)

    mc=np.matrix([1.5,1])

    mb=np.matrix([[pfpa(1.1,4.2,a_,b_),pfpb(1.1,4.2,a_,b_)],\
                  [pfpa(1.8,6.8,a_,b_),pfpb(1.8,6.8,a_,b_)],\
                  [pfpa(2.6,8.0,a_,b_),pfpb(2.6,8.0,a_,b_)]])

    mw = np.matrix([[f0(1.1,4.2,a_,b_)],\
                     [f0(1.8,6.8,a_,b_)],\
                     [f0(2.8,8.0,a_,b_)]])

    mwx = np.matrix([f0(1.5,6,a_,b_)])

    sig0 = 1
    dxdy=np.matrix(np.zeros(shape=(6,6)))
    dxdy[0,0]=1
    dxdy[0,2]=-0.5
    dxdy[1,1]=2
    dxdy[2,0]=-0.5
    dxdy[2,2]=1
    dxdy[3,3]=1
    dxdy[4,4]=1
    dxdy[5,5]=2
    pv=sig0*(dxdy)**(-1)

    naa = ma * pv**(-1) * ma.T
    nbb = mb.T * naa**(-1) * mb
    ncc = mc * nbb**(-1) * mc.T
    we = mb.T * naa**(-1) * mw

    xhat = -(nbb**(-1) - nbb**(-1)*mc.T*ncc**(-1)*mc*nbb**(-1))\
             *we-nbb**(-1)*mc.T*ncc**(-1)*mwx
    vmaj = -pv**(-1)*ma.T * naa**(-1)*(mw+mb*xhat)

    
    a_=a_+xhat[0,0]
    b_=b_+xhat[1,0]
    print(a_,b_)

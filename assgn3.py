'''
A network consists of distance and angle.
Helmert variance component estimation is tested.
Note：
DMS were converted to decimal degree by hand

x = [xc,yc,xb,yb] t = 4
v = [beta1 ... beta8, edge1 ... edge5]
'''

import numpy as np
RHOSEC = 206264.81
TO_M = 100

def getAz(xa,ya,xb,yb):
    dx = xb - xa
    dy = yb - ya
    az = 0
    azbase = np.arctan(dy/dx)
    if dx>0 and dy>0:
        az = azbase
    elif dx<0 and dy>0:
        az = np.pi + azbase
    elif dx<0 and dy<0:
        az = np.pi + azbase
    else:
        az = 2 * np.pi + azbase
    return az
def getAng(pta,pt1,pt2):
    az1 = getAz(pta[0], pta[1], pt1[0], pt1[1])
    az2 = getAz(pta[0], pta[1], pt2[0], pt2[1])
    ang = az2 - az1
    if ang < 0:
        ang += 2 * np.pi
    return ang * RHOSEC

def getDist(xa,ya,xb,yb,n = 2):
    dx = xb - xa
    dy = yb - ya
    return (dx*dx +dy*dy) **(1/n)
def getDistp(pta,ptb,n = 2):
    return getDist(pta[0], pta[1], ptb[0], ptb[1])

def pAzpx(xa,ya,xb,yb,nkn=[1,1]):
    '''
    linearize Az = f(x)
    Parameters
    ----------
    xa,ya : from point cord
    xb,tb : from point cord
    nkn : whether point is known 0 or unknown 1
    if two points are input the suppose to be a len = 2 list
    '''
    dx = xb - xa
    dy = yb - ya
    sqdis = getDist(xa, ya, xb, yb,1)
    px1_ =  dy/sqdis
    py1_ = -dx/sqdis
    px1 = nkn[0] * px1_
    py1 = nkn[0] * py1_
    px2 = nkn[1] * (-px1_)
    py2 = nkn[1] * (-py1_)
    return [px1 * RHOSEC,py1 * RHOSEC,px2*RHOSEC,py2*RHOSEC]

# reverse-clockwise 1 to 2
def pBetapx(xa,ya,x1,y1,x2,y2,nkn=[1,1,1]):
    b1 = pAzpx(xa, ya, x1, y1,[nkn[0],nkn[1]])
    b2 = pAzpx(xa, ya, x2, y2,[nkn[0],nkn[2]])
    return [b2[0]-b1[0],b2[1]-b1[1],-b1[2],-b1[3],b2[2],b2[3]]

def subpBeta(bij6,indexes):
    r = [0.0,0.0,0.0,0.0,0.0,0.0]
    for it in range(6):
        if indexes[it] > 0:
            r[indexes[it]-1] = bij6[it]
    if len(r) == 2:
        r = r + [0.0,0.0]
    return r[0:4]

def pEdgepx(x1,y1,x2,y2,nkn = [1,1]):
    dx = x2 - x1
    dy = y2 - y1
    dis = getDist(x1, y1, x2, y2)
    r = [-dx/dis*nkn[0],-dy/dis*nkn[0],dx/dis*nkn[1],dy/dis*nkn[1]]
    return r

def powerEdge(dist,sig):
    return (sig/(0.001 + 1e-6 * dist))**2

def helmertEst(n1:int,n2:int,nt:int,mb:np.matrix,mp:np.matrix,xhat,lmin):
    mb1 = mb[0:n1,0:nt]
    mp1 = mp[0:n1,0:n1]
    mb2 = mb[n1:n1+n2,0:nt]
    mp2 = mp[n1:n1+n2,n1:n1+n2]
    mn1 = mb1.T * mp1 * mb1
    mn2 = mb2.T * mp2 * mb2
    mn = mn1 + mn2
    mni = mn **(-1)

    lmin1 = lmin[0:n1,0]
    lmin2 = lmin[n1:n1+n2]
    v1 = mb1 * xhat - lmin1
    v2 = mb2 * xhat - lmin2
    mtheta = np.matrix([[(v1.T * mp1 *v1)[0,0]],[(v2.T * mp2 *v2)[0,0]]])
    # msmaj = np.matrix([\
    #     [n1 - 2*np.trace(mni*mn1) + np.trace(mni*mn1*mni*mn1) , np.trace(mni*mn1*mni*mn2)],\
    #     [np.trace(mni*mn1*mni*mn2) , n2 - 2*np.trace(mni*mn2) + np.trace(mni*mn2*mni*mn2)]])
    '''
    elementa not at diag my be negative, causing result sigma0 to be negative when the level of
    two types of observations differs too much, e.g. asec and meter. asec may work well
    with millimeter
    '''
    msmaj = np.matrix([\
        [n1 - 2*np.trace(mni*mn1) + np.trace(mni*mn1*mni*mn1) , 0],\
        [0 , n2 - 2*np.trace(mni*mn2) + np.trace(mni*mn2*mni*mn2)]])
    r = msmaj**(-1)*mtheta
    # print(r)
    return r

if __name__ == "__main__":
    sig0beta = 1    #unit asec
    sig0 = 0.001        #unit cm


    ptA = 564231.7913, 499937.7026
    ptB = 565858.0833, 499247.9980
    ptC = 564029.6555, 499613.4205 
    ptD = 565549.0958, 498813.5240
    
    for itvar in range(20):   # iteration of variance components conv. at ~13
        for itnl in range(5):  # iteration for non linear
        # by adding +- to make xc yc xx yd align
            mb = np.matrix([\
                subpBeta( pBetapx(ptB[0], ptB[1], ptA[0], ptA[1], ptC[0], ptC[1],nkn=[0,0,1]) , [0,0,0,0,1,2]), \
                subpBeta( pBetapx(ptB[0], ptB[1], ptC[0], ptC[1], ptD[0], ptD[1],nkn=[0,1,1]) , [0,0,1,2,3,4]), \
                subpBeta( pBetapx(ptA[0], ptA[1], ptC[0], ptC[1], ptD[0], ptD[1],nkn=[0,1,1]) , [0,0,1,2,3,4]), \
                subpBeta( pBetapx(ptA[0], ptA[1], ptD[0], ptD[1], ptB[0], ptB[1],nkn=[0,1,0]) , [0,0,3,4,0,0]), \
                subpBeta( pBetapx(ptC[0], ptC[1], ptD[0], ptD[1], ptB[0], ptB[1],nkn=[1,1,0]) , [1,2,3,4,0,0]), \
                subpBeta( pBetapx(ptC[0], ptC[1], ptB[0], ptB[1], ptA[0], ptA[1],nkn=[1,0,0]) , [1,2,0,0,0,0]), \
                subpBeta( pBetapx(ptD[0], ptD[1], ptB[0], ptB[1], ptA[0], ptA[1],nkn=[1,0,0]) , [3,4,0,0,0,0]), \
                subpBeta( pBetapx(ptD[0], ptD[1], ptA[0], ptA[1], ptC[0], ptC[1],nkn=[1,0,1]) , [3,4,0,0,1,2]), \
                pEdgepx(ptB[0], ptB[1],ptD[0], ptD[1], [0,1]),\
                pEdgepx(ptC[0], ptC[1],ptD[0], ptD[1], [1,1]),\
                pEdgepx(ptC[0], ptC[1],ptA[0], ptA[1], [1,0]),\
                pEdgepx(ptA[0], ptA[1],ptD[0], ptD[1], [0,1]),\
                pEdgepx(ptC[0], ptC[1],ptB[0], ptB[1], [1,0])\
                ])

            mlmaj = np.matrix([\
                42006.43,\
                237167.42,\
                293253.78,\
                62984.68,\
                59261.29,\
                249705.88,\
                305793.21,\
                45779.12,\
                533.1374,\
                1717.1362,\
                382.1227,\
                1731.7819,\
                1864.6199,\
                ]).T

            the1 = sig0beta
            the2 = sig0
            the12 = the1 * the1
            mp = np.matrix(np.diag([the12,the12,the12,the12,the12,the12,the12,the12]+[\
                powerEdge(533.1374,the2),powerEdge(1717.1362,the2),powerEdge(382.1227,the2),\
                powerEdge(1731.7819,the2),powerEdge(1864.6199,the2)]))
            mx0 = np.matrix([ptC[0],ptC[1],ptD[0],ptD[1]]).T

            mw = np.matrix([\
                getAng(ptB, ptA, ptC),\
                getAng(ptB, ptC, ptD),\
                getAng(ptA, ptC, ptD),\
                getAng(ptA, ptD, ptB),\
                getAng(ptC, ptD, ptB),\
                getAng(ptC, ptB, ptA),\
                getAng(ptD, ptB, ptA),\
                getAng(ptD, ptA, ptC),\
                getDistp(ptB,ptD),\
                getDistp(ptC,ptD),\
                getDistp(ptC,ptA),\
                getDistp(ptA,ptD),\
                getDistp(ptC,ptB),\
                ]).T
            
            mn = mb.T * mp * mb
            dx = mn **(-1) * mb.T * mp * (mlmaj - mw)
            mx = mx0 + dx
            ptC = mx[0,0],mx[1,0]
            ptD = mx[2,0],mx[3,0]

            if np.max(np.fabs(dx)) < 1e-7:  #iteration check
                print(f"iteration for non-linear ends at{itnl+1}")
                print(mx)
                break
        
        thetas = helmertEst(8, 5, 4, mb, mp, dx, mlmaj - mw)
        sig0 = thetas[1,0] ** (1/2) / thetas[0,0] ** (1/2)
        print(f"sig in {itvar} times iteration for variation equals {sig0}")

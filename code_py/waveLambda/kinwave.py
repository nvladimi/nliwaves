# this is a comment

import numpy as np

eps=1e-12

#-- define functions --

def kfun(p, q, r):
    
    k = p + q - r
    return(k)

def Qfun(k,p,q):

    #print("Qfun arguments:", k, p, q)
    
    pp = sum(p*p)
    qq = sum(q*q)
    kk = sum(k*k)
    kp = sum(p*k)

    ss = kk * pp

    if ss < eps:
        return(0)
        
    Q = (qq/ss)**(1/8) * (np.sqrt(ss) + kp)

    return(Q)

def Ufun(k,p,q):
    
    U = - Qfun(-k,p,q) - Qfun(-k,q,p) +  Qfun(p,q,-k)
    return(U)

def U1fun(k,p,q):
    
    U1 = Qfun(k,p,q) + Qfun(k,q,p) +  Qfun(p,q,k)
    return(U1)

def Afun(k,p,q):

    U = Ufun(k,p,q)

    #print("Afun arguments, U:", k, p, q, U)

    if np.abs(U) < eps:
        return(0)
        
    denom = (sum(k*k))**(1/4) - (sum(p*p))**(1/4) - (sum(q*q))**(1/4)
    A = - U / denom
    return(A)

def A1fun(k,p,q):

    U1 = Ufun(k,p,q)

    #print("A1fun arguments, U1:", k, p, q, U1)

    if np.abs(U1) < eps:
        return(0)
    
    denom = (sum(k*k))**(1/4) + (sum(p*p))**(1/4) + (sum(q*q))**(1/4)
    A1 = - U1 / denom
    return(A1)


def Zfun(k0,k1,k2,k3):
    
    Z =  U1fun( k0,    k1, -k0-k1 ) * A1fun( k2,    k3, -k2-k3 ) \
        + Ufun( k0+k1, k0,  k1    ) *  Afun( k2+k3, k2,  k3    ) \
        - Ufun( k0,    k2,  k0-k2 ) *  Afun( k3,    k1,  k3-k1 ) \
        - Ufun( k2,    k0,  k2-k0 ) *  Afun( k1,    k3,  k1-k3 ) \
        - Ufun( k0,    k3,  k0-k3 ) *  Afun( k2,    k1,  k2-k1 ) \
        - Ufun( k3,    k0,  k3-k0 ) *  Afun( k1,    k2,  k1-k2 )
    return(Z)
        
def Vfun(k0,k1,k2,k3):
    
    kk0 = sum(k0**2)
    kk1 = sum(k1**2)

    if kk0<eps:
        return(0)
    if kk1<eps:
        return(0)
    
    kk2 = sum(k2**2)
    kk3 = sum(k3**2)
    
    k02 = np.sqrt( sum((k0+k2)**2) )
    k12 = np.sqrt( sum((k1+k2)**2) )
    k03 = np.sqrt( sum((k0+k3)**2) )
    k13 = np.sqrt( sum((k1+k3)**2) )

    v = 2 * kk0 * np.sqrt(kk1) + 2 * kk1 * np.sqrt(kk0) \
        - np.sqrt(kk0*kk1) * (k02 + k12 + k03 +k13)
    V = (kk2 * kk3 / kk0 / kk1)**(1/8) * v
    return(V)

def Wfun(p,q,k,r):
    
    W =   Vfun(-p,-q,k,r) + Vfun(k,r,-p,-q) - Vfun(-p,k,-q,r) \
        - Vfun(-q,k,-p,r) - Vfun(-p,r,-q,k) - Vfun(-q,r,-p,k)
    return(W)

def Lfun(p,q,k,r):

    #print("Lfun arguments:", p, q, k, r)

    L = Wfun(p,q,k,r) + \
        0.25*( Zfun(p,q,k,r) + Zfun(q,p,k,r) +  Zfun(k,r,p,q) + Zfun(r,k,p,q) )
    return(L)

#------------------------------------------------------------------------------

def Nfun(r,q,p,p1):
    
    #print(r, q, p, p1)
    
    k  = p + q - r
    q1 = p + q - p1

    pp = sum(p*p)
    qq = sum(q*q)
    rr = sum(r*r)
    kk = sum(k*k)

    denom = pp**(0.25) + qq**(0.25) - rr**(0.25) - kk**(0.25)
    if np.abs(denom) < eps:
        return(0)

    L1 = Lfun(p,q,k,r)
    if np.abs(L1) < eps:
        return(0)
    
    L2 = Lfun(k,r,p1,q1)
    if np.abs(L2) < eps:
        return(0)

    L3 = Lfun(p,q,p1,q1)
    if np.abs(L3) < eps:
        return(0)

    #print("L1, L2, L3 = ", L1, L2, L3)

    c = ( rr**(-0.125) + kk**(-0.125) ) / denom

    #print("rr, kk, denom, c = ", rr, kk, denom, c)
    
    N = c * L1 * L2 / L3

    #print("N = ", N)
    
    return(N)

#------------------------------------------------------------------------------

def NLfun(r,q,p,p1):
    
    #print(r, q, p, p1)
    
    k  = p + q - r
    q1 = p + q - p1

    L1 = Lfun(p,q,k,r)
    if np.abs(L1) < eps:
        return(0)
    
    L2 = Lfun(k,r,p1,q1)
    if np.abs(L2) < eps:
        return(0)

    L3 = Lfun(p,q,p1,q1)
    if np.abs(L3) < eps:
        return(0)
    
    NL = L1 * L2 / L3
    
    return(NL)

#------------------------------------------------------------------------------

def NFfun(r,q,p,p1):
    
    #print(r, q, p, p1)
    
    k  = p + q - r

    pp = sum(p*p)
    qq = sum(q*q)
    rr = sum(r*r)
    kk = sum(k*k)
    
    if np.abs(rr) < eps:
        return(0)
    
    if np.abs(kk) < eps:
        return(0)
     
    denom = pp**(0.25) + qq**(0.25) - rr**(0.25) - kk**(0.25)
    
    if np.abs(denom) < eps:
        return(0)

    NF = ( rr**(-0.125) + kk**(-0.125) ) / denom
   
    return(NF)

#------------------------------------------------------------------------------

def denomFfun(r,q,p,p1):
    
    #print(r, q, p, p1)
    
    k  = p + q - r

    pp = sum(p*p)
    qq = sum(q*q)
    rr = sum(r*r)
    kk = sum(k*k)
        
    denom = pp**(0.25) + qq**(0.25) - rr**(0.25) - kk**(0.25)
    
    return(denom)


#------------------------------------------------------------------------------

def L4Dfun(m, mr, fname):

    L4D = np.zeros((2*m+1,2*m+1,2*m+1,2*m+1))
    
    r  = np.array((mr, 0))
    
    for q0 in np.arange(-m,m+1):
        for q1 in np.arange(-m,m+1):
            
            q = np.array((q0, q1))

            for p0 in np.arange(-m,m+1):
                for p1 in np.arange(-m,m+1):

                    p = np.array((p0, p1))
                    k  = p + q - r

                    L4D[q0+m,q1+m,p0+m,p1+m] =  Lfun(p,q,k,r)

    L4D = L4D/mr**3

    with open(fname, 'wb') as f:
        np.save(f, np.array([m, mr]))
        np.save(f, L4D)

#------------------------------------------------------------------------------

def L2Dfun(m, dq, pp, fname):
    
    r  = np.array((1, 0))

    n = len(pp)
    
    L2D = np.zeros((n,2*m+1,2*m+1))
    
    for i in np.arange(n):
        p = pp[i]
    
        for q0 in np.arange(-m,m+1):
            for q1 in np.arange(-m,m+1):
            
                q = np.array((q0, q1))*dq
                k  = p + q - r

                L2D[i,q0+m,q1+m] =  Lfun(p,q,k,r)

    with open(fname, 'wb') as f:
        np.save(f, np.array([m, dq]))
        np.save(f, pp)
        np.save(f, L2D)
        

#------------------------------------------------------------------------------

def Nsum(p,p1, q, mr, dr, arrout = False):

    NF = np.zeros((2*mr+1, 2*mr+1))
    NL = np.zeros((2*mr+1, 2*mr+1))
    Nr = np.zeros((2*mr+1, 2*mr+1))

    for i0 in np.arange(-mr,mr+1):
        for i1 in np.arange(-mr,mr+1):
            
            r = np.array((i0,i1))*dr
            #Nr[mr+i0,mr+i1] = Nfun(r,q,p,p1)
            NL[mr+i0,mr+i1] = NLfun(r,q,p,p1)
            NF[mr+i0,mr+i1] = NFfun(r,q,p,p1)
    
    S = np.sum(NL*NF)*dr*dr

    if arrout:
        return(S, NF, NL)
    else:
        return(S)


#------------------------------------------------------------------------------

def N1sum(p,p1, q, mr, dr, arrout = False):

    NF = np.zeros((2*mr+1, 2*mr+1))

    for i0 in np.arange(-mr,mr+1):
        for i1 in np.arange(-mr,mr+1):
            
            r = np.array((i0,i1))*dr
            NF[mr+i0,mr+i1] = NFfun(r,q,p,p1)
    
    S = np.sum(NF)*dr*dr

    if arrout:
        return(S, NF)
    else:
        return(S)
    
#------------------------------------------------------------------------------

def integratePV(p, p1, q, mr, dr, arrout):

    debug = True
    drmin = dr/32
    
    #-- compute denominator field --
    
    D = np.zeros((2*mr+1, 2*mr+1))

    for i0 in np.arange(-mr,mr+1):
        for i1 in np.arange(-mr,mr+1):
            
            r = np.array((i0,i1))*dr
            D[mr+i0,mr+i1] = denomFfun(r,q,p,p1)

    #-- find and classify cells where denominator changes sign --

    type1 =  set({})
    type2 =  set({})
    type3 =  set({})
    type4 =  set({})
    type5 =  set({})
    type6 =  set({})

    signD = np.sign(D)

    change0 = signD[:-1, :] * signD[1:, :]
    change1 = signD[:, :-1] * signD[:, 1:]

    ind0 = np.argwhere(change0 < 0)
    ind1 = np.argwhere(change1 < 0)

    for i in ind0:
        for j in ind1:
            if (i == j).all():
                type6.add(tuple(i))
            if (i == j-[1,0]).all():
                type5.add(tuple(i))
            if (i-[0,1] == j).all():
                type3.add(tuple(j))
            if (i+[1,-1] == j).all():
                type4.add(tuple(i-[0,1]))
        for j in ind0:
            if (i+[0,1] == j).all():
                type1.add(tuple(i))
    for i in ind1:
        for j in ind1:
            if (i+[1,0] == j).all():
                type2.add(tuple(i))
                                
    if 1==2:
        
        #for i in ind0:
        #    print(i, D[i[0],i[1]],  D[i[0]+1, i[1] ])
        #print(type(ind1))
        #for i in ind1:
        #    print(i, D[i[0],i[1]],  D[i[0], i[1]+1 ])
        #
        #print("number of crossings = ", len(ind0) + len(ind1))

        def print_corners(i):
            print("[{:3.0f}, {:3.0f}]  {:3.0f} {:3.0f} {:3.0f} {:3.0f}".format(i[0], i[1],
                  signD[i]*signD[i],              signD[i[0],i[1]+1]*signD[i],
                  signD[i[0]+1, i[1]+1]*signD[i], signD[i[0]+1,i[1]]*signD[i] ) )

        print("number type1 cells = ", len(type1) )
        for i in type1:
            print_corners(i)

        print("number type2 cells = ", len(type2) )
        for i in type2:
            print_corners(i)

        print("number type3 cells = ", len(type3) )
        for i in type3:
            print_corners(i)

        print("number type4 cells = ", len(type4) )
        for i in type4:
            print_corners(i)

        print("number type5 cells = ", len(type5) )
        for i in type5:
            print_corners(i)

        print("number type6 cells = ", len(type6) )
        for i in type6:
            print_corners(i)

        print("number of crossings = ", len(ind0) + len(ind1))
 
        print("number of cells = ",
              len(type1) + len(type2) + len(type3) +len(type4) + len(type5) + len(type6) ) 


    #-- find adjustments to sum in each cells with interface --

    def evalFun(x,y):
        r = np.array((x,y))
        # v = NFfun(r,q,p,p1)
        v = 1/denomFfun(r,q,p,p1)
        # v = 5*x - 4*y + 3
        # v=1
        #v = evalDenom(x,y)
        return(v)
    
    def evalDenom(x,y):
        r = np.array((x,y))
        v = denomFfun(r,q,p,p1)
        return(v)

    def findrootX(x0,x1,y):
        dx = x1 - x0
        f0 = evalDenom(x0,y)
        f1 = evalDenom(x1,y)
        #print("findrootX:  ", f0,f1)
        while dx > drmin:
            xm = (x0+x1)/2
            fm = evalDenom(xm,y)
            if (f0*fm) < 0:
                x1 = xm
                f1 = fm
            else:
                x0 = xm
                f0 = fm
            dx = dx/2
        #print("findrootX: ", fm)
        return(xm)

    def findrootY(y0,y1,x):
        dy = y1 - y0
        f0 = evalDenom(x,y0)
        f1 = evalDenom(x,y1)
        #print("findrootY:  ", f0,f1)
        while dy > drmin:
            ym = (y0+y1)/2
            fm = evalDenom(x,ym)
            if (f0*fm) < 0:
                y1 = ym
                f1 = fm
            else:
                y0 = ym
                f0 = fm
            dy = dy/2
        #print("findrootY: ", fm)
        return(ym)

    
    def cellcoords(cell,dr):
        return((cell[0]-mr)*dr, (cell[0]-mr+1)*dr, (cell[1]-mr)*dr, (cell[1]-mr+1)*dr)

    def celleval(x0,x1,y0,y1):
        a = evalFun(x0,y0)
        b = evalFun(x0,y1)
        c = evalFun(x1,y1)
        d = evalFun(x1,y0)
        return(a,b,c,d)

    def offcut(x0,y0,x1,y1):
        dx = np.abs(x0-x1)
        dy = np.abs(y0-y1)
        l  = np.sqrt(dx*dx + dy*dy)
        dx = dy/l * drmin*2
        dy = dx/l * drmin*2
        return(dx, dy)
    
    def volume12():
        v1 = a * dr * (A  + F1  + D)
        v2 = d * dr * (D  + F1 + G1)
        v3 = b * dr * (G2 + F2  + B)
        v4 = c * dr * (B  + C  + G2)
        dV = (v1 + v2 + v3 + v4)/6 - dr*dr*(A+B+C+D)/4
        return(dV)

    def volume36():
        v1 = a * c  * (A  + F1  + G1)
        v2 = b * dr * (B  + C   + F2)
        v3 = d * dr * (D  + C   + G2)
        v4 = (a*dr + c*dr - a*c) * (C + F2  + G2)
        dV = (v1 + v2 + v3 + v4)/6 - dr*dr*(A+B+C+D)/4
        return(dV)
             
    for cell in type1:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        X0 = findrootX(x0, x1, y0) 
        X1 = findrootX(x0, x1, y1)
        dx, dy = offcut(X0,y0, X1,y1)
        if X0<X1:
            b = x1 - X0
            a = X0 - x0
            c = x1 - X1
            d = X1 - x0
            A, D, C, B = celleval(x0,x1,y0,y1)
            G1 = evalFun(X1-dx, y1+dy)
            G2 = evalFun(X1+dx, y1-dy)
            F1 = evalFun(X0-dx, y0+dy)
            F2 = evalFun(X0+dx, y0-dy)
            dV = volume12()
            #print(  np.sign(np.array((A, D, F1, G1,  C, B, F2, G2))*np.sign(A)).astype("int") )
            #print(D, G1, G2, C)
        else:
            a = x1 - X0
            b = X0 - x0
            d = x1 - X1
            c = X1 - x0
            B, C, D, A = celleval(x0,x1,y0,y1)
            G1 = evalFun(X1+dx, y1+dy)
            G2 = evalFun(X1-dx, y1-dy)
            F1 = evalFun(X0+dx, y0+dy)
            F2 = evalFun(X0-dx, y0-dy)
            dV = volume12()
            #print(  np.sign(np.array((A, D, F1, G1,  C, B, F2, G2))*np.sign(A)).astype("int") )
        #print(dV, dr*dr*(A+B+C+D)/4 )
 
            
    for cell in type2:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        Y0 = findrootY(y0, y1, x0) 
        Y1 = findrootY(y0, y1, x1)
        dx, dy = offcut(x0,Y0, x1,Y1)
        if Y0<Y1:
            b = y1 - Y0
            a = Y0 - y0
            c = y1 - Y1
            d = Y1 - y0
            A, B, C, D = celleval(x0,x1,y0,y1)
            F1 = evalFun(x0+dx, Y0-dy)
            F2 = evalFun(x0-dx, Y0+dy)
            G1 = evalFun(x1+dx, Y1-dy)
            G2 = evalFun(x1-dx, Y1+dy)
            dV = volume12()
            print(  np.sign(np.array((A, D, F1, G1,  C, B, F2, G2))*np.sign(A)).astype("int") )

        else:
            c = y1 - Y0
            d = Y0 - y0
            b = y1 - Y1
            a = Y1 - y0
            D, C, B, A = celleval(x0,x1,y0,y1)
            G1 = evalFun(x0-dx, Y0-dy)
            G2 = evalFun(x0+dx, Y0+dy)
            F1 = evalFun(x1-dx, Y1-dy)
            F2 = evalFun(x1+dx, Y1+dy)
            dV = volume12()
            #print(  np.sign(np.array((A, D, F1, G1,  C, B, F2, G2))*np.sign(A)).astype("int") )
        #print(dV, dr*dr*(A+B+C+D)/4 )


            
    for cell in type3:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        Y0 = findrootY(y0, y1, x0) 
        X1 = findrootX(x0, x1, y1)
        dx, dy = offcut(x0,Y0, X1,y1)
        b = x1 - X1
        a = X1 - x0       
        c = y1 - Y0
        d = Y0 - y0
        D, A, B, C = celleval(x0,x1,y0,y1)
        F1 = evalFun(X1-dx, y1+dy)
        F2 = evalFun(X1+dx, y1-dy)
        G1 = evalFun(x0-dx, Y0+dy)
        G2 = evalFun(x0+dx, Y0-dy)
        dV = volume36()
        #print(dV, dr*dr*(A+B+C+D)/4 )

    for cell in type4:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        X1 = findrootX(x0, x1, y1)
        Y1 = findrootY(y0, y1, x1)
        dx, dy = offcut(X1,y1, x1,Y1)
        c = x1 - X1
        d = X1 - x0       
        a = y1 - Y1
        b = Y1 - y0
        C, D, A, B = celleval(x0,x1,y0,y1)
        F1 = evalFun(x1+dx, Y1+dy)
        F2 = evalFun(x1-dx, Y1-dy)
        G1 = evalFun(X1+dx, y1+dy)
        G2 = evalFun(X1-dx, y1-dy)
        dV = volume36()
        #print(dV, dr*dr*(A+B+C+D)/4 )

    for cell in type5:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        X0 = findrootX(x0, x1, y0)
        Y1 = findrootY(y0, y1, x1) 
        dx, dy = offcut(X0,y0, x1,Y1)
        a = x1 - X0
        b = X0 - x0       
        d = y1 - Y1
        c = Y1 - y0
        B, C, D, A = celleval(x0,x1,y0,y1)
        F1 = evalFun(X0+dx, y0-dy)
        F2 = evalFun(X0-dx, y0+dy)
        G1 = evalFun(x1+dx, Y1-dy)
        G2 = evalFun(x1-dx, Y1+dy)
        dV = volume36()
        #print(dV, dr*dr*(A+B+C+D)/4 )
        
    for cell in type6:
        x0, x1, y0, y1 = cellcoords(cell,dr)
        X0 = findrootX(x0, x1, y0)
        Y0 = findrootY(y0, y1, x0)
        dx, dy = offcut(x0,Y0, X0,y0)
        d = x1 - X0
        c = X0 - x0       
        b = y1 - Y0
        a = Y0 - y0
        A, B, C, D = celleval(x0,x1,y0,y1)
        F1 = evalFun(x0-dx, Y0-dy)
        F2 = evalFun(x0+dx, Y0+dy)
        G1 = evalFun(X0-dx, y0-dy)
        G2 = evalFun(X0+dx, y0+dy)
        dV = volume36()
        #print(dV, dr*dr*(A+B+C+D)/4 )

        
    #S = np.sum(NF)*dr*dr
    S = 0
    return(S, signD)



#------------------------------------------------------------------------------

def N2Darr(p,p1, mq, dq, mr, dr):

    Nq = np.zeros((2*mq+1, 2*mq+1))

    for i0 in np.arange(-mq,mq+1):
        for i1 in np.arange(-mq,mq+1):

            q = np.array((i0,i1))*dq

            Nq[mq+i0, mq+i1] = Nsum(p,p1,q, mr, dr)
    
    return(Nq)

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

def write_N2D_pi2(nq, nr):
    
    p  = [1, 0]
    p1 = [0, 1]

    fname = "sumN_pi2_q" + str(nq).zfill(2) + "r" + str(nr).zfill(2)+ ".npy"

    dq = 1/nq
    dr = 1/nr
    mq = nq*nq
    mr = nr*nr
    
    Nlambda = N2Darr(np.array(p), np.array(p1), mq, dq, mr, dr)

    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, (mq, dq, mr, dr))
        np.save(f, Nlambda)

        
def write_N2D_pi4(nq, nr):
    
    p  = [1, 0]
    p1 = [1, 1]/np.sqrt(2)

    fname = "sumN_pi4_q" + str(nq).zfill(2) + "r" + str(nr).zfill(2)+ ".npy"

    dq = 1/nq
    dr = 1/nr
    mq = nq*nq
    mr = nr*nr
    
    Nlambda = N2Darr(np.array(p), np.array(p1), mq, dq, mr, dr)

    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, (mq, dq, mr, dr))
        np.save(f, Nlambda)



#------------------------------------------------------------------------------

def N2Dfld_pi2(q, mr, dr, fname):

    p  = [1, 0]
    p1 = [1, 1]

    S, NF, NL = Nsum(np.array(p), np.array(p1), np.array(q), mr, dr, arrout=True)
    
    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, q)
        np.save(f, (mr, dr, S))
        np.save(f, NF)
        np.save(f, NL)

        
def N2Dfld_pi4(q, mr, dr, fname):

    p  = [1, 0]
    p1 = [1, 1]/np.sqrt(2)

    S, NF, NL = Nsum(np.array(p), np.array(p1), np.array(q), mr, dr, arrout=True)
    
    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, q)
        np.save(f, (mr, dr, S))
        np.save(f, NF)
        np.save(f, NL)


        
def N1fld(q, mr, dr, fname):

    p  = [1, 0]
    p1 = [0, 1]

    S, NF = N1sum(np.array(p), np.array(p1), np.array(q), mr, dr, arrout=True)
    
    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, q)
        np.save(f, (mr, dr, S))
        np.save(f, NF)

#------------------------------------------------------------------------------
        
def testNfld(q, mr, dr, fname):

    p  = [1, 0]
    p1 = [0, 1]

    S, NF = integratePV(np.array(p), np.array(p1), np.array(q), mr, dr, arrout=True)
    
    with open(fname, 'wb') as f:
        np.save(f, p)
        np.save(f, p1)
        np.save(f, q)
        np.save(f, (mr, dr, S))
        np.save(f, NF)

#------------------------------------------------------------------------------

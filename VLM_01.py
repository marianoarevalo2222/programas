import numpy as np
k=0.5
c=1.0
xyn=[[[0.25*c/2,-2],[0.25*c/2,-1],[0.25*c/2,0],[0.25*c/2,1],[0.25*c/2,2]],[[0.25*c/2+c/2,-2],[0.25*c/2+c/2,-1],[0.25*c/2+c/2,0],[0.25*c/2+c/2,1],[0.25*c/2+c/2,2]]]
xy=[[0.75*c/2,-3*k],[0.75*c/2,-k],[0.75*c/2,k],[0.75*c/2,3*k],[0.75*c/2+c/2,-3*k],[0.75*c/2+c/2,-k],[0.75*c/2+c/2,k],[0.75*c/2+c/2,3*k]]

VAB=[]
VAinf=[]
VBinf=[]
r0=[]
r1=[]
r2=[]
kr=0
k1=0

for i in range(0,len(xyn)) :
    for j in range(0,len(xyn[i])-1) :
        r0.append([xyn[i][j+1][0]-xyn[i][j][0],xyn[i][j+1][1]-xyn[i][j][1],0])        
        #print(r0[kr])
        kr=kr+1

for ip in range(0,len(xy)) :
    for i in range(0,len(xyn)) :
        for j in range(0,len(xyn[i])-1) :
            r1.append([xy[ip][0]-xyn[i][j][0],xy[ip][1]-xyn[i][j][1],0])
            r2.append([xy[ip][0]-xyn[i][j+1][0],xy[ip][1]-xyn[i][j+1][1],0])
            #print(r2[k1])
            k1=k1+1
r00=[]
for jp in range(0,k1) :
    for krr in range(0,kr) :
        r00.append(r0[krr])
        
for kp in range(0,k1) :
    crossv=(np.cross(r1[kp],r2[kp]))
    abscross=np.dot(crossv,crossv)
    r11=np.dot(r1[kp],r1[kp])
    r22=np.dot(r2[kp],r2[kp])
    escv=np.dot(r00[kp],(r1[kp]/r11-r2[kp]/r22))
    VAB.append(crossv*escv/abscross/4/3.14)
    VAinf.append((1.0+r1[kp][0]/r11**0.5)/(-r1[kp][1])/4/3.14)
    VBinf.append((1.0+r2[kp][0]/r22**0.5)/(-r2[kp][1])/4/3.14)
    #print(VBinf[kp],kp)

Acoef=np.zeros((len(xy),len(xy)))
bc=[]
vinf=100
alpha=0.1
for i in range(0,len(xy)) :
    	bc.append(-vinf*alpha)

for i in range(0,len(xy)) :
    j=0
    for kp in range(0+i*(len(xy)),len(xy)+i*(len(xy))) :
        Acoef[i][j]=VAB[kp][2]+VAinf[kp]+VBinf[kp]
        #print(i,j,Acoef[i][j],kp)
        j=j+1

#print(Acoef)

InvAcoef=np.linalg.inv(Acoef)
sol=np.dot(InvAcoef,bc)
print(sol)

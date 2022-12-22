import cmath
import math
from statistics import mean

import numpy as np

alfaV  = 10
alfaI  = 1*alfaV
gamma  = 10
beta   = 0.001
lam0   = 0.5
cf     = beta

# Reference Values IEEE14
angVR   = np.array([0,-0.1349,-0.3295,-0.2655,-0.2294,-0.379,-0.3429,-0.3429,-0.3836,-0.3899,-0.3872,-0.3998,-0.4002,-0.4172])
vR      = np.array([1.06,1.045,1.01,1.0095,1.018,1.07,1.0305,1.09,1.0067,1.0071,1.0331,1.0456,1.0357,0.993])

Pij     = np.array([0.5776,0.1141,0.027,0.2583,0.1147,0.0629,0.0645,0.122,-0.0885,0.3843,2.4038,-0.9935,-0.3253,1.114,0.8641,0.7806,0.6439,0.3843,0.2155,0])
Qij     = np.array([-0.0181,0.048,0.0222,0.1525,0.1441,0.113,-0.0284,-0.0041,-0.0781,0.2303,-0.3782,0.1355,0.1248,0.029,-0.0604,-0.0165,-0.1619,-0.0862,0.0159,0.3684])
Pji     = np.array([-0.5602,-0.1124,-0.0268,-0.2531,-0.1119,-0.0616,-0.0644,-0.1201,0.0909,-0.3843,-2.3021,1.0401,0.3336,-1.0541,-0.8544,-0.7482,-0.6439,-0.3843,-0.2152,0])
Qji     = np.array([0.035,-0.0446,-0.022,-0.1422,-0.1382,-0.11,0.0288,0.0081,0.083,-0.2095,0.6304,0.0146,-0.1389,0.1649,0.0777,0.0754,0.2691,0.118,0.0096,-0.3483])

         
#----------------------------------------------------------------------------------------------------

# Model Parameters
bus_f   = np.array([2,6,12,6,6,11,9,9,14,7,1,3,3,1,5,2,5,4,4,8]) 
bus_to  = np.array([5,12,13,13,11,10,10,14,13,9,2,2,4,5,4,4,6,7,9,7])
rLine   = np.array([0.057,0.1229,0.2209,0.0662,0.095,0.0821,0.0318,0.1271,0.1709,0,0.0194,0.047,0.067,0.054,0.0134,0.0581,0,0,0.005,0])
xLine   = np.array([0.1739,0.2558,0.1999,0.1303,0.1989,0.1921,0.0845,0.2704,0.348,0.11,0.0592,0.198,0.171,0.223,0.0421,0.1763,0.252,0.2091,0.5562,0.1762])
bLine   = np.array([0.034,0,0,0,0,0,0,0,0,0,0.0528,0.0438,0.0346,0.0492,0.0128,0.0374,0,0,0,0])
lineStat= np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])

tr_f    = np.array([])
tr_to   = np.array([])
rTr     = np.array([])
xTr     = np.array([])
bTr     = np.array([])
#a       = np.array([])

bus_f   = np.append(bus_f,  tr_f,   axis=None)
bus_to  = np.append(bus_to, tr_to,  axis=None)

rLine   = np.append(rLine, rTr, axis=None)
xLine   = np.append(xLine, xTr, axis=None)
bLine   = np.append(bLine, bTr, axis=None)

nBusN    = 14
numLines = 20 

#----------------------------------------------------------------------------------------------------

vRe    =   np.zeros((nBusN))
vIm    =   np.zeros((nBusN))
xV     =   np.zeros((nBusN))
xAngV  =   np.zeros((nBusN))

IRMag  =   np.zeros((2*numLines))
angIR  =   np.zeros((2*numLines))
xI     =   np.zeros((2*numLines))
xAngI  =   np.zeros((2*numLines))

ILines  =   np.zeros((nBusN,nBusN),dtype = 'complex_')
yLines  =   np.zeros((nBusN,nBusN),dtype = 'complex_')
bLines  =   np.zeros((nBusN,nBusN),dtype = 'complex_')
#aTr     =   np.zeros((max(max(tr_f),max(tr_to)),max(max(tr_f),max(tr_to))))


for i in range(len(vR)):
    vRe[i] = cmath.rect(vR[i],angVR[i]).real
    vIm[i] = cmath.rect(vR[i],angVR[i]).imag


for i in range(numLines):
    fromBus = int(bus_f[i])-1
    toBus   = int(bus_to[i])-1

    ILines[fromBus][toBus] = np.conj( (Pij[i] + 1j*Qij[i]) / (vRe[fromBus] + 1j*vIm[fromBus]) )
    ILines[toBus][fromBus] = np.conj( (Pji[i] + 1j*Qji[i]) / (vRe[toBus]   + 1j*vIm[toBus]) )

    yLines[fromBus][toBus] = 1/(rLine[i]*lineStat[i] + 1j*xLine[i]*lineStat[i])
    yLines[toBus][fromBus] = yLines[fromBus][toBus]

    bLines[fromBus][toBus] = 1j*bLine[i]
    bLines[toBus][fromBus] = bLines[fromBus][toBus]


#for i in range(len(a)):
#    fromBusTr = tr_f[i]-1
#    toBusTr   = tr_to[i]-1
#    aTr[fromBusTr][toBusTr] = a[i] 
#    aTr[toBusTr][fromBusTr] = a[i]


IRe = np.zeros((np.count_nonzero(ILines)))
IIm = np.zeros((np.count_nonzero(ILines)))
ii  = 0
for i in range(nBusN):
    for j in range(nBusN):
        if yLines[i][j] != 0: 
            IRe[ii] = ILines[i][j].real
            IIm[ii] = ILines[i][j].imag
            ii=ii+1           


for i in range(len(IRe)):
    IRMag[i] = cmath.polar(IRe[i] + 1j*IIm[i])[0]
    angIR[i] = cmath.polar(IRe[i] + 1j*IIm[i])[1]

#---------------------------------------------------------------------------------------------------
#PARAMETROS AREAS

#Definicion de areas

area1 = np.array([1,2,3,4,5])
area2 = np.array([5,6,9,10,11,12,13,14])
area3 = np.array([4,7,8,9])
front = np.array([ [0,5,4] , [5,0,9] , [4,9,0] ])

nBus  = np.array([len(area1),len(area2),len(area3)])
    
refBus = [1,6,7]

yLines1  = np.zeros((nBus[0],nBus[0]) , dtype='complex_')
bLines1  = np.zeros((nBus[0],nBus[0]) , dtype='complex_')
for i in range(len(area1)):
    for j in range(len(area1)):
        yLines1[i][j] = yLines[area1[i]-1][area1[j]-1]
        bLines1[i][j] = bLines[area1[i]-1][area1[j]-1]

yLines2  = np.zeros((nBus[1],nBus[1]) , dtype='complex_')
bLines2  = np.zeros((nBus[1],nBus[1]) , dtype='complex_')
for i in range(len(area2)):
    for j in range(len(area2)):
        yLines2[i][j] = yLines[area2[i]-1][area2[j]-1]
        bLines2[i][j] = bLines[area2[i]-1][area2[j]-1]

yLines3  = np.zeros((nBus[2],nBus[2]) , dtype='complex_')
bLines3  = np.zeros((nBus[2],nBus[2]) , dtype='complex_')
for i in range(len(area3)):
    for j in range(len(area3)):
        yLines3[i][j] = yLines[area3[i]-1][area3[j]-1]
        bLines3[i][j] = bLines[area3[i]-1][area3[j]-1]

nLines =  [int(np.count_nonzero(yLines1)/2) , int(np.count_nonzero(yLines2)/2) , int(np.count_nonzero(yLines3)/2) ]
     
corrientesArea1 = np.array([1,2,3,4,5,6,7,8,9,10,11,14,15,16])
corrientesArea2 = np.array([17,18,19,20,21,28,29,30,31,32,33,34,35,36,37,38,39,40])
corrientesArea3 = np.array([12,13,22,23,24,25,26,27])

#---------------------------------------------------------------------------------------------------
# MEASUREMENTS

sigMagV = 0.02
sigAngV = 0.005
sigMagI = 0.02
sigAngI = 0.005

erV     = np.random.normal(0, sigMagV, vR.size)
erAngV  = np.random.normal(0, sigAngV, vR.size)
erI     = np.random.normal(0, sigMagI, IRMag.size)
erAngI  = np.random.normal(0, sigAngI, IRMag.size)

vM      = vR + erV
angVM   = angVR + erAngV
vMRe    =   np.zeros((nBusN))
vMIm    =   np.zeros((nBusN))
for i in range(len(vM)):
    vMRe[i] = cmath.rect(vM[i],angVM[i]).real
    vMIm[i] = cmath.rect(vM[i],angVM[i]).imag

IM = IRMag + erI
angIM = angIR + erAngI
IMRe    =   np.zeros((2*numLines))
IMIm    =   np.zeros((2*numLines))
for i in range(len(IM)):
    IMRe[i] = cmath.rect(IM[i],angIM[i]).real
    IMIm[i] = cmath.rect(IM[i],angIM[i]).imag

Z = np.concatenate([vMRe, vMIm, IMRe, IMIm] , dtype=object)


#Estimador area 1 ----------------------------------------------

varMagV1 = np.ones(nBus[0])*sigMagV**2
varAngV1 = np.ones(nBus[0])*sigAngV**2
varMagI1 = np.ones(2*nLines[0])*sigMagI**2
varAngI1 = np.ones(2*nLines[0])*sigAngI**2

varianza1   = np.concatenate([varMagV1,varAngV1,varMagI1,varAngI1] , dtype=object)
Wpol1       = np.diag(varianza1)

A1      = np.zeros((2*nLines[0],nBus[0]) , dtype=object)
ySerie1 = np.zeros((2*nLines[0],2*nLines[0]) , dtype=object)

k=0
for i in range(nBus[0]):
    for j in range(nBus[0]):
        if  abs(yLines1[i][j]) != 0:
            rowA = np.zeros(nBus[0])
            rowA[i] = 1
            rowA[j] = -1
            A1[k][:] = rowA
            ySerie1[k][k] = yLines1[i][j]
            k=k+1

VMatrix1 = np.eye(nBus[0])

Yio1 = np.zeros(2*nLines[0], dtype = 'complex_')
ii  = 0
for i in range(nBus[0]):
    for j in range(nBus[0]):
        if abs(yLines1[i][j]) != 0: 
            Yio1[ii] = (bLines1[i][j])/2
            ii=ii+1    

medMap1 = np.zeros(2*nLines[0])
yShunt1  = np.zeros((2*nLines[0],nBus[0]), dtype = 'complex_')

for i in range(2*nLines[0]):
    medMap1[i] = np.where(A1[i][:]==1)[0][0]
    yShunt1[i][int(medMap1[i])] = Yio1[i]

M1 = np.vstack(ySerie1.dot(A1) + yShunt1)

M1real = np.zeros(M1.shape)
M1imag = np.zeros(M1.shape)
for i in range(M1.shape[0]):
    for j in range(M1.shape[1]):
        M1real[i][j] = (M1[i][j]).real
        M1imag[i][j] = (M1[i][j]).imag

# Submatrices de columnas de B
b1ColsMats1 = np.concatenate([VMatrix1                ,  M1real , np.zeros(VMatrix1.shape), M1imag ])
b1ColsMats2 = np.concatenate([np.zeros(VMatrix1.shape), -M1imag , VMatrix1                , M1real ])
B1 =  np.hstack((b1ColsMats1,b1ColsMats2))

W1 = Wpol1

B12 =  np.vstack([ B1[4][:] ,
                  [B1[nBus[0] + 11][:] + B1[nBus[0] + 12][:] + B1[nBus[0] + 13][:]][0],
                  [B1[2*nBus[0] + 2*nLines[0] - 1][:]][0],
                  [B1[2*nBus[0] + 2*nLines[0] + 11][:] + B1[2*nBus[0] + 2*nLines[0] + 12][:] + B1[2*nBus[0] + 2*nLines[0] + 13][:]][0],
                ])

B13 =  np.vstack([ B1[3][:] ,
                  [B1[nBus[0] + 8][:] + B1[nBus[0] + 9][:] + B1[nBus[0] + 10][:]][0],
                  [B1[nBus[0] + 2*nLines[0] + 3][:]][0],
                  [B1[2*nBus[0] + 2*nLines[0] + 8][:] + B1[2*nBus[0] + 2*nLines[0] + 9][:] + B1[2*nBus[0] + 2*nLines[0] + 10][:]][0],
                ])

W12 =np.diag(np.ones(4)*np.sqrt(2/beta))
W13 =np.diag(np.ones(4)*np.sqrt(2/beta))

BVirt1 = np.vstack([ B1,
                    B12,    
                    B13 
                  ])

WVirt1 = np.vstack([ np.hstack([W1                                   , np.zeros((W1.shape[0],W12.shape[1]))      , np.zeros((W1.shape[0],W13.shape[1]))  ]),
                     np.hstack([np.zeros((W12.shape[0],W1.shape[1])) , W12                                       , np.zeros((W12.shape[0],W13.shape[1])) ]),
                     np.hstack([np.zeros((W13.shape[0],W1.shape[1])) ,np.zeros((W13.shape[0],W12.shape[1]))      , W13                                   ])])



WVirt1 = WVirt1.astype('float')
WVirt1_inv = np.linalg.inv(WVirt1)

H1Virt = (np.linalg.inv( (BVirt1.T).dot( WVirt1_inv.dot(BVirt1) ) )).dot( (BVirt1.T).dot(WVirt1_inv) )




#Estimador area 2 ----------------------------------------------

varMagV2 = np.ones(nBus[1])*sigMagV**2
varAngV2 = np.ones(nBus[1])*sigAngV**2
varMagI2 = np.ones(2*nLines[1])*sigMagI**2
varAngI2 = np.ones(2*nLines[1])*sigAngI**2
varianza2= np.concatenate([varMagV2,varAngV2,varMagI2,varAngI2] , dtype=object)
Wpol2    = np.diag(varianza2)

A2     = np.zeros((2*nLines[1],nBus[1]) , dtype=object)
ySerie2 = np.zeros((2*nLines[1],2*nLines[1]) , dtype=object)
k=0
for i in range(nBus[1]):
    for j in range(nBus[1]):
        if  abs(yLines2[i][j]) != 0:
            rowA = np.zeros(nBus[1])
            rowA[i] = 1
            rowA[j] = -1
            A2[k][:] = rowA
            ySerie2[k][k] = yLines2[i][j]
            k=k+1

VMatrix2 = np.eye(nBus[1])

Yio2 = np.zeros(2*nLines[1], dtype = 'complex_')
ii  = 0
for i in range(nBus[1]):
    for j in range(nBus[1]):
        if abs(yLines2[i][j]) != 0: 
            Yio2[ii] = (bLines2[i][j])/2
            ii=ii+1    

medMap2 = np.zeros(2*nLines[1])
yShunt2  = np.zeros((2*nLines[1],nBus[1]), dtype = 'complex_')

for i in range(2*nLines[1]):
    medMap2[i] = np.where(A2[i][:]==1)[0][0]
    yShunt2[i][int(medMap2[i])] = Yio2[i]

M2 = np.vstack(ySerie2.dot(A2) + yShunt2)

M2real = np.zeros(M2.shape)
M2imag = np.zeros(M2.shape)
for i in range(M2.shape[0]):
    for j in range(M2.shape[1]):
        M2real[i][j] = (M2[i][j]).real
        M2imag[i][j] = (M2[i][j]).imag

# Submatrices de columnas de B
b2ColsMats1 = np.concatenate([VMatrix2                ,  M2real , np.zeros(VMatrix2.shape), M2imag ])
b2ColsMats2 = np.concatenate([np.zeros(VMatrix2.shape), -M2imag , VMatrix2                , M2real ])
B2 =  np.hstack((b2ColsMats1,b2ColsMats2))

W2 = Wpol2

B21 =  np.vstack([ B2[0][:] ,
                  [B2[nBus[1]][:]][0] ,
                  [B2[nBus[1] + 2*nLines[1]][:]][0],
                  [B2[2*nBus[1] + 2*nLines[1]][:]][0],
                ])

B23 =  np.vstack([ B2[2][:] ,
                  [B2[nBus[1] + 5][:] + B2[nBus[1] + 6][:]][0],
                  [B2[nBus[1] + 2*nLines[1] + 2][:]][0],
                  [B2[2*nBus[1] + 2*nLines[1] + 5][:] + B2[2*nBus[1] + 2*nLines[1] + 6][:]][0],
                ])

W21 =np.diag(np.ones(4)*np.sqrt(2/beta))
W23 =np.diag(np.ones(4)*np.sqrt(2/beta))

BVirt2 = np.vstack([ B2,
                     B21,    
                     B23 
                  ])

WVirt2 = np.vstack([ np.hstack([W2                                   , np.zeros((W2.shape[0],W21.shape[1]))      , np.zeros((W2.shape[0],W23.shape[1]))  ]),
                     np.hstack([np.zeros((W21.shape[0],W2.shape[1])) , W21                                       , np.zeros((W12.shape[0],W13.shape[1])) ]),
                     np.hstack([np.zeros((W23.shape[0],W2.shape[1])) , np.zeros((W23.shape[0],W21.shape[1]))     , W23                                   ])])

WVirt2 = WVirt2.astype('float')
WVirt2_inv = np.linalg.inv(WVirt2)

H2Virt = (np.linalg.inv( (BVirt2.T).dot( WVirt2_inv.dot(BVirt2) ) )).dot( (BVirt2.T).dot(WVirt2_inv))



#Estimador area 3 ----------------------------------------------

varMagV3 = np.ones(nBus[2])*sigMagV**2
varAngV3 = np.ones(nBus[2])*sigAngV**2
varMagI3 = np.ones(2*nLines[2])*sigMagI**2
varAngI3 = np.ones(2*nLines[2])*sigAngI**2
varianza3= np.concatenate([varMagV3,varAngV3,varMagI3,varAngI3] , dtype=object)
Wpol3    = np.diag(varianza3)

A3     = np.zeros((2*nLines[2],nBus[2]) , dtype=object)
ySerie3 = np.zeros((2*nLines[2],2*nLines[2]) , dtype=object)
k=0
for i in range(nBus[2]):
    for j in range(nBus[2]):
        if  abs(yLines3[i][j]) != 0:
            rowA = np.zeros(nBus[2])
            rowA[i] = 1
            rowA[j] = -1
            A3[k][:] = rowA
            ySerie3[k][k] = yLines3[i][j]
            k=k+1

VMatrix3 = np.eye(nBus[2])

Yio3 = np.zeros(2*nLines[2], dtype = 'complex_')
ii  = 0
for i in range(nBus[2]):
    for j in range(nBus[2]):
        if abs(yLines3[i][j]) != 0: 
            Yio3[ii] = (bLines3[i][j])/2
            ii=ii+1    

medMap3 = np.zeros(2*nLines[2])
yShunt3  = np.zeros((2*nLines[2],nBus[2]), dtype = 'complex_')

for i in range(2*nLines[2]):
    medMap3[i] = np.where(A3[i][:]==1)[0][0]
    yShunt3[i][int(medMap3[i])] = Yio3[i]

M3 = np.vstack(ySerie3.dot(A3) + yShunt3)

M3real = np.zeros(M3.shape)
M3imag = np.zeros(M3.shape)
for i in range(M3.shape[0]):
    for j in range(M3.shape[1]):
        M3real[i][j] = (M3[i][j]).real
        M3imag[i][j] = (M3[i][j]).imag

# Submatrices de columnas de B
b3ColsMats1 = np.concatenate([VMatrix3                ,  M3real , np.zeros(VMatrix3.shape), M3imag ])
b3ColsMats2 = np.concatenate([np.zeros(VMatrix3.shape), -M3imag , VMatrix3                , M3real ])
B3 =  np.hstack((b3ColsMats1,b3ColsMats2))

W3 = Wpol3

B31 =  np.vstack([ B3[0][:] ,
                  [B3[nBus[2]][:] + B3[nBus[2] + 1 ][:]][0] ,
                  [B3[nBus[2] + 2*nLines[2]][:]][0],
                  [B3[2*nBus[2] + 2*nLines[2]][:] + B3[2*nBus[2] + 2*nLines[2] +1 ][:]][0],
                ])

B32 =  np.vstack([ B3[3][:] ,
                  [B3[nBus[2] + 6][:] + B3[nBus[2] + 7][:]][0],
                  [B3[nBus[2] + 2*nLines[2] + 3][:]][0],
                  [B3[2*nBus[2] + 2*nLines[2] + 6][:] + B3[2*nBus[2] + 2*nLines[2] + 7][:]][0],
                ])

W31 =np.diag(np.ones(4)*np.sqrt(2/beta))
W32 =np.diag(np.ones(4)*np.sqrt(2/beta))

BVirt3 = np.vstack([ B3,
                     B31,    
                     B32 
                  ])

WVirt3 = np.vstack([ np.hstack([W3                                   , np.zeros((W3.shape[0],W31.shape[1]))      , np.zeros((W3.shape[0],W32.shape[1]))  ]),
                     np.hstack([np.zeros((W31.shape[0],W3.shape[1])) , W31                                       , np.zeros((W31.shape[0],W32.shape[1])) ]),
                     np.hstack([np.zeros((W32.shape[0],W3.shape[1])) , np.zeros((W32.shape[0],W31.shape[1]))     , W32                                   ])])

WVirt3 = WVirt3.astype('float')
WVirt3_inv = np.linalg.inv(WVirt3)

H3Virt = (np.linalg.inv( (BVirt3.T).dot( WVirt3_inv.dot(BVirt3) ) )).dot( (BVirt3.T).dot(WVirt3_inv))


# % MEDICIONES POR AREA

#Area 1--------------------------
vM1 = vM[(area1-1)]
angVM1 = angVM[(area1-1)]
vMRe1    =   np.zeros((len(vM1)))
vMIm1    =   np.zeros((len(vM1)))
for i in range(len(vM1)):
    vMRe1[i] = cmath.rect(vM1[i],angVM1[i]).real
    vMIm1[i] = cmath.rect(vM1[i],angVM1[i]).imag

IM1 = IM[(corrientesArea1-1)]
angIM1 = angIM[(corrientesArea1-1)]
IMRe1    =   np.zeros((len(IM1)))
IMIm1    =   np.zeros((len(IM1)))
for i in range(len(IM1)):
    IMRe1[i] = cmath.rect(IM1[i],angIM1[i]).real
    IMIm1[i] = cmath.rect(IM1[i],angIM1[i]).imag

Z1 = np.concatenate([vMRe1, IMRe1, vMIm1, IMIm1])

#Area 2--------------------------

vM2     = vM[(area2-1)]
angVM2  = angVM[(area2-1)]
vMRe2   = np.zeros((len(vM2)))
vMIm2   = np.zeros((len(vM2)))
for i in range(len(vM2)):
    vMRe2[i] = cmath.rect(vM2[i],angVM2[i]).real
    vMIm2[i] = cmath.rect(vM2[i],angVM2[i]).imag

IM2     = IM[(corrientesArea2-1)]
angIM2  = angIM[(corrientesArea2-1)]
IMRe2   = np.zeros((len(IM2)))
IMIm2   = np.zeros((len(IM2)))
for i in range(len(IM2)):
    IMRe2[i] = cmath.rect(IM2[i],angIM2[i]).real
    IMIm2[i] = cmath.rect(IM2[i],angIM2[i]).imag

Z2 = np.concatenate([vMRe2, IMRe2, vMIm2, IMIm2])

#Area 3--------------------------

vM3     = vM[(area3-1)]
angvM3  = angVM[(area3-1)]
vMRe3   = np.zeros((len(vM3)))
vMIm3   = np.zeros((len(vM3)))
for i in range(len(vM3)):
    vMRe3[i] = cmath.rect(vM3[i],angvM3[i]).real
    vMIm3[i] = cmath.rect(vM3[i],angvM3[i]).imag

IM3     = IM[(corrientesArea3-1)]
angIM3  = angIM[(corrientesArea3-1)]
IMRe3   = np.zeros((len(IM3)))
IMIm3   = np.zeros((len(IM3)))
for i in range(len(IM3)):
    IMRe3[i] = cmath.rect(IM3[i],angIM3[i]).real
    IMIm3[i] = cmath.rect(IM3[i],angIM3[i]).imag

Z3 = np.concatenate([vMRe3, IMRe3, vMIm3, IMIm3])

# Corrientes fronteras

Ibus5_1 = (M1[11][:] + M1[12][:] + M1[13][:]).dot(vMRe1 + 1j*vMIm1)
Ibus4_1 = (M1[8][:]  + M1[9][:]  + M1[10][:]).dot(vMRe1 + 1j*vMIm1)

Ibus5_2 =            (M2[0][:]).dot(vMRe2 + 1j*vMIm2)
Ibus9_2 = (M2[5][:] + M2[6][:]).dot(vMRe2 + 1j*vMIm2)
   
Ibus4_3 = (M3[0][:] + M3[1][:]).dot(vMRe3 + 1j*vMIm3)
Ibus9_3 = (M3[6][:] + M3[7][:]).dot(vMRe3 + 1j*vMIm3)

lam1 = np.array([0,1,0,1,0,1,0,1])*lam0
lam2 = np.array([0,1,0,1,0,1,0,1])*lam0
lam3 = np.array([0,1,0,1,0,1,0,1])*lam0

iter=0
diferk = 1

xV1     =   np.zeros((nBus[0]))
xAngV1  =   np.zeros((nBus[0]))
xI1     =   np.zeros((2*nLines[0]))
xAngI1  =   np.zeros((2*nLines[0]))

xV2     =   np.zeros((nBus[1]))
xAngV2  =   np.zeros((nBus[1]))
xI2     =   np.zeros((2*nLines[1]))
xAngI2  =   np.zeros((2*nLines[1]))

xV3     =   np.zeros((nBus[2]))
xAngV3  =   np.zeros((nBus[2]))
xI3     =   np.zeros((2*nLines[2]))
xAngI3  =   np.zeros((2*nLines[2]))

dif =[]

while diferk > 1e-6 :

    lam1k = lam1
    lam2k = lam2
    lam3k = lam3

    Zvirt1_2 = [    (beta*vMRe1[4]          - gamma*vMRe1[4]        + gamma*vMRe2[0]        - lam1[0])/beta,
                    cf*(beta*Ibus5_1.real   - gamma*Ibus5_1.real    + gamma*Ibus5_2.real    - lam1[1])/beta,
                    (beta*vMIm1[4]          - gamma*vMIm1[4]        + gamma*vMIm2[0]        - lam1[2])/beta,
                    cf*(beta*Ibus5_1.imag   - gamma*Ibus5_1.imag    + gamma*Ibus5_2.imag    - lam1[3])/beta
                ]

    Zvirt1_3 = [    (beta*vMRe1[3]          - gamma*vMRe1[3]        + gamma*vMRe3[0]        - lam1[4])/beta,
                    cf*(beta*Ibus4_1.real   - gamma*Ibus4_1.real    + gamma*Ibus4_3.real    - lam1[5])/beta,
                    (beta*vMIm1[3]          - gamma*vMIm1[3]        + gamma*vMIm3[0]        - lam1[6])/beta,
                    cf*(beta*Ibus4_1.imag   - gamma*Ibus4_1.imag    + gamma*Ibus4_3.imag    - lam1[7])/beta
                ]

    Z1Virt = np.concatenate([Z1, Zvirt1_2, Zvirt1_3], dtype=object)

    
    xEst1 = H1Virt.dot(Z1Virt)

    xVRe1 = xEst1[0       : nBus[0]]
    xVIm1 = xEst1[nBus[0] :    :   ]

    
    for i in range(nBus[0]):
        xV1[i]    = cmath.polar(xVRe1[i] + 1j*xVIm1[i])[0]
        xAngV1[i] = cmath.polar(xVRe1[i] + 1j*xVIm1[i])[1]

    
    xIRe1 = (np.hstack((M1real,-M1imag))).dot(np.concatenate([xVRe1,xVIm1]))
    xIIm1 = (np.hstack((M1imag, M1real))).dot(np.concatenate([xVRe1,xVIm1]))

    for i in range(2*nLines[0]):
        xI1[i]    = cmath.polar(xIRe1[i] + 1j*xIIm1[i])[0]
        xAngI1[i] = cmath.polar(xIRe1[i] + 1j*xIIm1[i])[1]

    #-----------------------------------------

    Zvirt2_1 = [    (beta*vMRe2[0]          - gamma*vMRe2[0]        + gamma*vMRe1[4]        - lam2[0])/beta,
                    cf*(beta*Ibus5_2.real   - gamma*Ibus5_2.real    + gamma*Ibus5_1.real    - lam2[1])/beta,
                    (beta*vMIm2[0]          - gamma*vMIm2[0]        + gamma*vMIm1[4]        - lam2[2])/beta,
                    cf*(beta*Ibus5_2.imag   - gamma*Ibus5_2.imag    + gamma*Ibus5_1.imag    - lam2[3])/beta
                ]

    Zvirt2_3 = [    (beta*vMRe2[2]          - gamma*vMRe2[2]        + gamma*vMRe3[3]        - lam2[4])/beta,
                    cf*(beta*Ibus9_2.real   - gamma*Ibus9_2.real    + gamma*Ibus9_3.real    - lam2[5])/beta,
                    (beta*vMIm2[2]          - gamma*vMIm2[2]        + gamma*vMIm3[3]        - lam2[6])/beta,
                    cf*(beta*Ibus9_2.imag   - gamma*Ibus9_2.imag    + gamma*Ibus9_3.imag    - lam2[7])/beta
                ]       

    Z2Virt = np.concatenate([Z2, Zvirt2_1, Zvirt2_3], dtype=object)

    
    xEst2 = H2Virt.dot(Z2Virt)

    xVRe2 = xEst2[0       : nBus[1]]
    xVIm2 = xEst2[nBus[1] :    :   ]

    
    for i in range(nBus[1]):
        xV2[i]    = cmath.polar(xVRe2[i] + 1j*xVIm2[i])[0]
        xAngV2[i] = cmath.polar(xVRe2[i] + 1j*xVIm2[i])[1]

    
    xIRe2 = (np.hstack((M2real,-M2imag))).dot(np.concatenate([xVRe2,xVIm2]))
    xIIm2 = (np.hstack((M2imag, M2real))).dot(np.concatenate([xVRe2,xVIm2]))

    for i in range(2*nLines[1]):
        xI2[i]    = cmath.polar(xIRe2[i] + 1j*xIIm2[i])[0]
        xAngI2[i] = cmath.polar(xIRe2[i] + 1j*xIIm2[i])[1]


    #----------------------------------

    Zvirt3_1 = [    (beta*vMRe3[0]          - gamma*vMRe3[0]        + gamma*vMRe1[3]        - lam3[0])/beta,
                    cf*(beta*Ibus4_3.real   - gamma*Ibus4_3.real    + gamma*Ibus4_1.real    - lam3[1])/beta,
                    (beta*vMIm3[0]          - gamma*vMIm3[0]        + gamma*vMIm1[3]        - lam3[2])/beta,
                    cf*(beta*Ibus4_3.imag   - gamma*Ibus4_3.imag    + gamma*Ibus4_1.imag    - lam3[3])/beta
                ]

    Zvirt3_2 = [    (beta*vMRe3[3]          - gamma*vMRe3[3]        + gamma*vMRe2[2]        - lam3[4])/beta,
                    cf*(beta*Ibus9_3.real   - gamma*Ibus9_3.real    + gamma*Ibus9_2.real    - lam3[5])/beta,
                    (beta*vMIm3[3]          - gamma*vMIm3[3]        + gamma*vMIm2[2]        - lam3[6])/beta,
                    cf*(beta*Ibus9_3.imag   - gamma*Ibus9_3.real    + gamma*Ibus9_2.imag    - lam3[7])/beta
                ]

    Z3Virt = np.concatenate([Z3, Zvirt3_1, Zvirt3_2], dtype=object)

    xEst3 = H3Virt.dot(Z3Virt)

    xVRe3 = xEst3[0       : nBus[2]]
    xVIm3 = xEst3[nBus[2] :    :   ]

    for i in range(nBus[2]):
        xV3[i]    = cmath.polar(xVRe3[i] + 1j*xVIm3[i])[0]
        xAngV3[i] = cmath.polar(xVRe3[i] + 1j*xVIm3[i])[1]

    
    xIRe3 = (np.hstack((M3real,-M3imag))).dot(np.concatenate([xVRe3,xVIm3]))
    xIIm3 = (np.hstack((M3imag, M3real))).dot(np.concatenate([xVRe3,xVIm3]))

    for i in range(2*nLines[2]):
        xI3[i]    = cmath.polar(xIRe3[i] + 1j*xIIm3[i])[0]
        xAngI3[i] = cmath.polar(xIRe3[i] + 1j*xIIm3[i])[1]

    Ibus5_1 = (M1[11][:] + M1[12][:] + M1[13][:]).dot(xVRe1 + 1j*xVIm1)
    Ibus5_2 = (M2[0][:]).dot(xVRe2 + 1j*xVIm2)

    Ibus4_1 = (M1[8][:]  + M1[9][:]  + M1[10][:]).dot(xVRe1 + 1j*xVIm1)
    Ibus4_3 = (M3[0][:] + M3[1][:]).dot(xVRe3 + 1j*xVIm3)
    
    Ibus9_2 = (M2[5][:] + M2[6][:]).dot(xVRe2 + 1j*xVIm2)
    Ibus9_3 = (M3[6][:] + M3[7][:]).dot(xVRe3 + 1j*xVIm3)

    xFront12 = [ xVRe1[4]       -   xVRe2[0]    ,
                Ibus5_1.real    +   Ibus5_2.real,
                xVIm1[4]        -   xVIm2[0]    ,
                Ibus5_1.imag    +   Ibus5_2.imag
               ]

    xFront13 = [ xVRe1[3]       -   xVRe3[0]    ,
                 Ibus4_1.real   +   Ibus4_3.real,
                 xVIm1[3]       -   xVIm3[0]    ,
                 Ibus4_1.imag   +   Ibus4_3.imag
                ]

    xFront23 = [ xVRe2[2]       -   xVRe3[3]    ,
                 Ibus9_2.real   +   Ibus9_3.real,
                 xVIm2[2]       -   xVIm3[3]    ,
                 Ibus9_2.imag   +   Ibus9_3.imag,
                ]

    zDif = np.concatenate([xFront12, xFront13, xFront23], dtype=object)

    difer = np.sqrt(sum(zDif**2)/len(zDif))
    

    iter=iter+1
    dif.append(difer)

    if abs(difer-diferk) < 1e-5:
        print('conv')
        break
                
    if iter > 2 and (difer-diferk) > 0:
        print('no conv')
        break
        

    diferk=difer;
    
    lam1 = lam1k + [alfaV*( xVRe1[4]        - xVRe2[0]      ),
                    alfaI*( Ibus5_1.real    + Ibus5_2.real  ),
                    alfaV*( xVIm1[4]        - xVIm2[0]      ),
                    alfaI*( Ibus5_1.imag    + Ibus5_2.imag  ),
                    alfaV*( xVRe1[3]        - xVRe3[0]      ),
                    alfaI*( Ibus4_1.real    + Ibus4_3.real  ),
                    alfaV*( xVIm1[3]        - xVIm3[0]      ),
                    alfaI*( Ibus4_1.imag    + Ibus4_3.imag  )]        

    lam2 = lam2k + [alfaV*(  xVRe2[0]       - xVRe1[4]      ),
                    alfaI*(  Ibus5_2.real   + Ibus5_1.real  ),
                    alfaV*(  xVIm2[0]       - xVIm1[4]      ),
                    alfaI*(  Ibus5_2.imag   + Ibus5_1.imag  ),
                    alfaV*(  xVRe2[2]       - xVRe3[3]      ),
                    alfaI*(  Ibus9_2.real   + Ibus9_3.real  ),
                    alfaV*(  xVIm2[2]       - xVIm3[2]      ),
                    alfaI*(  Ibus9_2.imag   + Ibus9_3.imag  )]

    lam3 = lam3k + [alfaV*(  xVRe3[0]       - xVRe1[3]      ),
                    alfaI*(  Ibus4_3.real   + Ibus4_1.real  ),
                    alfaV*(  xVIm3[0]       - xVIm1[3]      ),
                    alfaI*(  Ibus4_3.imag   + Ibus4_1.imag  ),
                    alfaV*(  xVRe3[3]       - xVRe2[2]      ),
                    alfaI*(  Ibus9_3.real   + Ibus9_2.real  ),
                    alfaV*(  xVIm3[3]       - xVIm2[2]      ),
                    alfaI*(  Ibus9_3.imag   + Ibus9_2.imag  )]


zDifV = [ [xV1[4] - xV2[0] , xAngV1[4] - xAngV2[0]],
          [xV1[3] - xV3[0] , xAngV1[3] - xAngV3[0]],
          [xV2[2] - xV3[3] , xAngV2[2] - xAngV3[3]]
        ]

print('zDifV')
print(zDifV)

MAEV = [ [mean(abs(xV1-vR[area1-1])) , mean(abs(xAngV1-angVR[area1-1]))],
         [mean(abs(xV2-vR[area2-1])) , mean(abs(xAngV2-angVR[area2-1]))],
         [mean(abs(xV3-vR[area3-1])) , mean(abs(xAngV3-angVR[area3-1]))]
        ]

print('MAEV')
print(MAEV)

#Computational Thermodynamics (CTD)
#Note: Function div(k*grad(T))=0 is defined based on influx positive principle
#User specify the mesh refinment in x and y direction
import numpy as np
xRefine = [0, 0.1, 0.3, 0.6]
yRefine = [0, 0.1, 0.2, 0.3]
TIni  = [1000]
kIni  = [[[0.1,0.6],[0,0.3],100],1e-3] # [[[xzone],[yzone],value1],default]
TIniBFX = [[[-1e-6,0.05],[-1e-6,0.31],1,[320]],
           [[0.59,0.61],[-1e-6,0.31],2,[0]]] #X direction face boundary conditions
TIniBFY = [[[-1e-6,0.05],[-1e-6,0.05],2,[0]],
           [[0.09,0.61],[-1e-6,0.05],2,[100]],
           [[-1e-6,0.61],[0.29,0.31],3,[300,20]]]
           #[[xzone],[yzone],type,[values]]
           # default to 0. internal face 1. Dirichlet specify T 2. Von Neumman specifcy fluxIn q (Watt/m^2)([0] for adiabatic) 3. Mixed specify far field convection Tinf hinf(convection coefficient)



#Field Evaluation
nVX = len(xRefine) #number of vertices in X
nVY = len(yRefine) #number of vertices in Y
nCel = (nVX-1)*(nVY-1) #number of cells
##Get a length of boundary field definition
lenBF = 0
for i in range(0,len(TIniBFX)):
    if lenBF < len(TIniBFX[i][3]):
        lenBF = len(TIniBFX[i][3])
for i in range(0,len(TIniBFY)):
    if lenBF < len(TIniBFY[i][3]):
        lenBF = len(TIniBFY[i][3])
lenBF += 1#save for type space
#Intialize the mesh points
yRefine,xRefine = np.meshgrid(yRefine,xRefine) #First index is X dir and the Second is Y dir [X index][Y index]

#Create mapping from cell number to X Y cell indices
##Forward mapping
def mapCel(i,nVX=nVX,nCel=nCel):
    if (i>=nCel) | (i<0):
        print("Error: Maping Indexing Out of Range")
        return None
    indX = int(np.mod(i+1,nVX-1)-1)
    indY = int(np.floor((i+1)/(nVX-1)))
    if indX == -1: #Last cell of the row
        indX = int(nVX-2)
        indY = int(indY-1)
    return [indX,indY]
##Reverse mapping
def reMapCel(indX,indY,nVX=nVX,nVY=nVY):
    if (indX >= (nVX-1)) | (indY >= (nVY-1)) | (indX < 0) | (indY < 0):
        print("Error: Reverse Mapping Indexing Out of Range")
        return None
    return int(indX+indY*(nVX-1))

#Calculate surface Area, centroid, distance, volume 
##X sweep
areaY = []#Area facing Y directions
centCellX = []#Cell centroid X locations
disCellX = []#Nodal distance in X direction (Stored distance to the right of the cell center)
for indX in range(0,nVX-1):
    areaY += [xRefine[indX+1][0] - xRefine[indX][0]] #[indX from 0 to nVX-2 for S and N]
    centCellX += [(xRefine[indX+1][0]+xRefine[indX][0])*0.5] #[indX from 0 to nVX-2 for x coordinate]
    if indX == 0: #For the left most one 
        disCellX += [centCellX[indX]-xRefine[0][0]] #Cell centroid X subtract out the x coord of left bottom corner
    else:
        disCellX += [centCellX[indX]-centCellX[indX-1]] #[indX from 0 to nVX-2: indX for W dist and indX+1 for E dist]
disCellX += [xRefine[-1][-1] - centCellX[indX]] #One more distance for the right most cell
##Y sweep
areaX = []#Area facing X direction
centCellY = []#Cell centroid Y locations
disCellY = []#Nodal distance in Y direction (Stored distance to the bottom of the cell center)
for indY in range(0,nVY-1):
    areaX += [yRefine[0][indY+1] - yRefine[0][indY]] #[indY from 0 to nVY-2 for W and E]
    centCellY += [(yRefine[0][indY+1] + yRefine[0][indY])*0.5] #[indY from 0 to nVY-2 for y coordinate]
    if indY == 0: #for the bottom most one
        disCellY += [centCellY[indY]-yRefine[0][0]] #Cell centroid Y subtract out the y coord of left bottom corner
    else:
        disCellY += [centCellY[indY]-centCellY[indY-1]] #[indY from 0 to nVY-2: indY for S dist and indX+1 for N dist]
disCellY += [yRefine[-1][-1] - centCellY[indY]] #One more distance for the  up most cell
##Cell sweep
voluCell = [] #1D array through all the cells 
for i in range(0,nCel):
    voluCell += [areaX[mapCel(i)[1]]*areaY[mapCel(i)[0]]] 

#Calculate interpolation factor
##X sweep
intCoefKX = []#Interpolation factor for k in X direction
for indX in range(0,nVX-1):
    if indX == 0: #For the left most one 
        intCoefKX += [0] #0 for the left most face, it takes the k from the left most cell
    else:
        intCoefKX += [voluCell[reMapCel(indX-1,0)]/(voluCell[reMapCel(indX,0)]+voluCell[reMapCel(indX-1,0)])] #[indX from 0 to nVX-2: indX for interW and indX+1 for interE]
intCoefKX += [1] #1 for the right most cace, it takes the k from the right most cell 
##Y sweep
intCoefKY = []#Interpolation factor for k in Y direction
for indY in range(0,nVY-1):
    if indY == 0: #For the bottom most one 
        intCoefKY += [0] #0 for the bottom most face, it takes the k from the bottom most cell
    else:
        intCoefKY += [voluCell[reMapCel(0,indY-1)]/(voluCell[reMapCel(0,indY)]+voluCell[reMapCel(0,indY-1)])] #[indY from 0 to nVY-2: indY for interS and indY+1 for interN]
intCoefKY += [1] #1 for the up most cace, it takes the k from the up most cell

#Calculate boundary conditions
#Set Temperature
if len(TIni) == 1:#Uniform field
    T = np.array([TIni[0]]*nCel)
else:
    T = np.array([TIni[-1]]*nCel)#Last entry is the default value
    for i in range(0,nCel):
        xCoor = centCellX[mapCel(i)[0]]
        yCoor = centCellY[mapCel(i)[1]]
        for j in range(0,len(TIni)-1):#Exclude the last one which is the default value
            if (xCoor>=TIni[j][0][0]) & (xCoor<=TIni[j][0][1]) & (yCoor>=TIni[j][1][0]) & (yCoor<=TIni[j][1][1]):
                T[i] = TIni[j][2]
#Set k conductivity (Side note hard coded here, need improvement)
if len(kIni) == 1:#Uniform field
    k = [kIni[0]]*nCel
else:
    k = np.array([kIni[-1]]*nCel)#Last entry is the default value
    for i in range(0,nCel):
        xCoor = centCellX[mapCel(i)[0]]
        yCoor = centCellY[mapCel(i)[1]]
        for j in range(0,len(kIni)-1):#Exclude the last one which is the default value
            if (xCoor>=kIni[j][0][0]) & (xCoor<=kIni[j][0][1]) & (yCoor>=kIni[j][1][0]) & (yCoor<=kIni[j][1][1]):
                k[i] = kIni[j][2]

#Interpolate coefficient to faces (Standard way face value will be stored in this code)
##X sweep
kX = np.zeros([nVX,nVY-1])#k in X direction face
for indX in range(0,nVX):
    for indY in range(0,nVY-1): #For x direction, we have nVX face in X direction but only nVY-1 rows in Y direction
        gCur = intCoefKX[indX]
        if gCur == 0: #left most
            kX[indX,indY] = k[reMapCel(indX,indY)]
        elif gCur == 1: #right most
            kX[indX,indY] = k[reMapCel(indX-1,indY)] #indX-1 go to the cell to the left
        else:
            kl = k[reMapCel(indX-1,indY)] #intCoefX is based on left
            kr = k[reMapCel(indX,indY)]
            kX[indX,indY] = (kl*kr)/(kl*(1-gCur)+kr*gCur)#Hamornic Inter to conserve flux
##Y sweep
kY = np.zeros([nVX-1,nVY])#k in Y direction face
for indX in range(0,nVX-1):
    for indY in range(0,nVY):
        gCur = intCoefKY[indY]
        if gCur == 0: #bottom most
            kY[indX,indY] = k[reMapCel(indX,indY)]
        elif gCur == 1: #up most
            kY[indX,indY] = k[reMapCel(indX,indY-1)] #go to the cell one level lower
        else:
            kb = k[reMapCel(indX,indY-1)] #intCoefY is based on the bottom
            ku = k[reMapCel(indX,indY)]
            kY[indX,indY] = (kb*ku)/(kb*(1-gCur)+ku*gCur) #Again not efficient here

#Setup boundary conditions to each face
##X sweep
TBFX = np.zeros([nVX,nVY-1,lenBF]) #Temperature infomation for X facing faces. Third element is for mix boundary cond
for indX in range(0,nVX):
    for indY in range(0,nVY-1): #For x direction, we have nVX face in X direction but only nVY-1 rows in Y direction
        xCoor = xRefine[indX][indY]
        yCoor = yRefine[indX][indY]
        for j in range(0, len(TIniBFX)):
            if (xCoor>=TIniBFX[j][0][0]) & (xCoor<=TIniBFX[j][0][1]) & (yCoor>=TIniBFX[j][1][0]) & (yCoor<=TIniBFX[j][1][1]):
                TBFX[indX][indY][0] = TIniBFX[j][2]
                TBFX[indX][indY][1:1+len(TIniBFX[j][3])] = np.array(TIniBFX[j][3])
##Y sweep For each face: First element is type following elements are the corresponding values T or q or Tinf hinf
TBFY = np.zeros([nVX-1,nVY,lenBF])#Temperature infomation for Y facing faces. Third element is for mix boundary cond
for indX in range(0,nVX-1):
    for indY in range(0,nVY):
        xCoor = xRefine[indX][indY]
        yCoor = yRefine[indX][indY]
        for j in range(0, len(TIniBFY)):
            if (xCoor>=TIniBFY[j][0][0]) & (xCoor<=TIniBFY[j][0][1]) & (yCoor>=TIniBFY[j][1][0]) & (yCoor<=TIniBFY[j][1][1]):
                TBFY[indX][indY][0] = TIniBFY[j][2]
                TBFY[indX][indY][1:1+len(TIniBFY[j][3])] = np.array(TIniBFY[j][3])


#Construct the matrix Ax = b
A = np.zeros([nCel,nCel])
b = np.zeros(nCel)
for i in range(0,nCel):
    indX = mapCel(i)[0]
    indY = mapCel(i)[1]
    indXNei = [indX,indX+1,indX,indX-1]
    indYNei = [indY-1,indY,indY+1,indY]
    indFac  = ['Y','X','Y','X']#South East North West
    indFacX = [indX,indX+1,indX,indX]
    indFacY = [indY,indY,indY+1,indY]
    for j in range(0,4):
        if indFac[j] == 'Y':
            TBFCur = TBFY[indFacX[j]][indFacY[j]]
            kCur = kY[indFacX[j]][indFacY[j]]
            areaCur = areaY[indFacX[j]]
            disCellCur = disCellY[indFacY[j]]
        else:
            TBFCur = TBFX[indFacX[j]][indFacY[j]]
            kCur = kX[indFacX[j]][indFacY[j]]
            areaCur = areaX[indFacY[j]]
            disCellCur = disCellX[indFacX[j]]
        KAD = kCur*areaCur/disCellCur
        if TBFCur[0] == 0: #Internal
            A[i][reMapCel(indXNei[j],indYNei[j])] += KAD
            A[i][i] += -KAD
        elif TBFCur[0] == 1: #Dirichlet
            b[i] += -TBFCur[1]*KAD
            A[i][i] += -KAD
        elif TBFCur[0] == 2: #Neumman
            b[i] += -TBFCur[1]*areaCur  
        else: #Mixing
            KD = kCur/disCellCur
            multiMix = (areaCur*TBFCur[2]*KD)/(TBFCur[2]+KD)
            b[i] += -multiMix*TBFCur[1]
            A[i][i] += -multiMix


##Normalize matrix
for i in range(0,nCel):
    b[i] = b[i]/A[i][i]
    A[i] = A[i]/A[i][i]

xExa = np.linalg.solve(A, b)

#Steepest Descent
TStep = []
numOL = 20000
for i in range(0,numOL):
    r = b - np.matmul(A,T)
    alpha = np.dot(r,r)/np.dot(r,np.matmul(A,r))
    T = T+alpha*r
    TStep += [T*1]


"""
##Print out the equations line by line to check with textbook solution
multiPrint = [0.005,240.002,340,0.006,440.002,640,0.006998,243.9624,345.9406]
roundTo = 60
for i in range(0,nCel):
    strPrint = 'T'+str(i+1)+' = '
    for j in range(0,len(A[i])):
        if (np.around(A[i][j],roundTo) != 0) & (j != i) :
            strPrint += str(-np.around(A[i][j],roundTo)*multiPrint[i])+' *T'+str(j+1)
            strPrint += ' + '
    strPrint += str(np.around(b[i],roundTo)*multiPrint[i])
    strPrint += '\n'
    print(strPrint)
"""


##Print out the results
for i in range(0,len(TStep),1):
    print(TStep[i])
    print('\n')

print('Exact\n')
print(xExa)
print('Found\n')
print(TStep[i])


print("end")

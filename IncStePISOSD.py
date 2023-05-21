#Incompressible Stead State Structural CFD with PISO Step

#Add Package
import numpy as np

############################################################################ User Input Section
#Domain Parameters
##User specify the mesh refinment in x and y direction
xRefine = [0., 0.3, 0.8, 2.] #Vertices along X direction
yRefine = [0., 0.2, 0.5, 1.] #Vertices along Y direction

#Boundary Parameters
##         [[xZoneLowBound, xZoneUpBound],[yZoneLowBound, yZoneUpBound],Boundary Type,[Boundary Values]]
##         Boundary Type 0 --> Internal Face (don't need to set this, default automatically)
##                       1 --> Dirichlet Face     Boundary Values --> The value fixed to (positive in +x and +y direction regardless inlet or outlet)
##                       2 --> Neumman Face       Boundary Values --> The gradient fixed to (positive in +x and +y direction regardless inlet or outlet)
##         Faces not covered will be automatically internal faces
UIniBFX = [[[0,0],[0,1],1,[2.]],
           [[2,2],[0,1],2,[0.]]] #X direction face boundary conditions
VIniBFX = [[[0,0],[0,1],1,[0.]],
           [[2,2],[0,1],2,[0.]]] #X direction face boundary conditions
PIniBFX = [[[0,0],[0,1],2,[0.]],
           [[2,2],[0,1],1,[78364.]]] #X direction face boundary conditions
UIniBFY = [[[0,2],[0,0],1,[0.]],
           [[0,2],[1,1],1,[0.]]]
VIniBFY = [[[0,2],[0,0],1,[0.]],
           [[0,2],[1,1],1,[0.]]]
PIniBFY = [[[0,2],[0,0],2,[0.]],
           [[0,2],[1,1],2,[0.]]]

#Initial Parameters
##      [ [[xZone1LowBound, xZone1UpBound],[yZone1LowBound, yZone1UpBound], valueForZone1]  ,  [[xZone2LowBound, xZone2UpBound],[yZone2LowBound, yZone2UpBound], valueForZone2]  , valueForAllOtherZone]
##      valueForAllOtherZone is always required as the last element, all other zones definition are optional
UIni  = [1.8]
VIni  = [0]
PIni  = [78364.]
nuIni  = [1.5e-5] 

#Numerical Parameters
numOut = 100    #Outer loop iteration
numIn  = 100    #Inner loop iteration
intSav = 10     #Save interval (iteration)
numSD  = 5000    #Iteration in steepest descent
facRelOut = 0.1 #Relaxation factor for momentum prediction step [0:No progress] [1:Full Progress]
facRelIn = 0.1  #Relaxation factor for pressure correction step [0:No progress] [1:Full Progress]
###########################################################################


#Field Evaluation
nVX = len(xRefine) #number of vertices in X
nVY = len(yRefine) #number of vertices in Y
nCel = (nVX-1)*(nVY-1) #number of cells
##Get a length of boundary field definition
lenBF = len(UIniBFX[0][3]) + 1
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
def setInternal(UIni,nCel=nCel,centCellX=centCellX,centCellY=centCellY):
    if len(UIni) == 1:#Uniform field
        U = np.array([UIni[0]]*nCel)
    else:
        U = np.array([UIni[-1]]*nCel)#Last entry is the default value
        for i in range(0,nCel):
            xCoor = centCellX[mapCel(i)[0]]
            yCoor = centCellY[mapCel(i)[1]]
            for j in range(0,len(UIni)-1):#Exclude the last one which is the default value
                if (xCoor>=UIni[j][0][0]) & (xCoor<=UIni[j][0][1]) & (yCoor>=UIni[j][1][0]) & (yCoor<=UIni[j][1][1]):
                    U[i] = UIni[j][2]
    return U
U  = setInternal(UIni)
V  = setInternal(VIni)
P  = setInternal(PIni)
nu = setInternal(nuIni)

#Interpolate coefficient to faces 
##X sweep
nuX = np.zeros([nVX,nVY-1])#k in X direction face
for indX in range(0,nVX):
    for indY in range(0,nVY-1): #For x direction, we have nVX face in X direction but only nVY-1 rows in Y direction
        gCur = intCoefKX[indX]
        if gCur == 0: #left most
            nuX[indX,indY] = nu[reMapCel(indX,indY)]
        elif gCur == 1: #right most
            nuX[indX,indY] = nu[reMapCel(indX-1,indY)] #indX-1 go to the cell to the left
        else:
            nul = nu[reMapCel(indX-1,indY)] #intCoefX is based on left
            nur = nu[reMapCel(indX,indY)]
            nuX[indX,indY] = nur*gCur+nul*(1-gCur)#Hamornic Inter to conserve flux
##Y sweep
nuY = np.zeros([nVX-1,nVY])#k in Y direction face
for indX in range(0,nVX-1):
    for indY in range(0,nVY):
        gCur = intCoefKY[indY]
        if gCur == 0: #bottom most
            nuY[indX,indY] = nu[reMapCel(indX,indY)]
        elif gCur == 1: #up most
            nuY[indX,indY] = nu[reMapCel(indX,indY-1)] #go to the cell one level lower
        else:
            nub = nu[reMapCel(indX,indY-1)] #intCoefY is based on the bottom
            nuu = nu[reMapCel(indX,indY)]
            nuY[indX,indY] = nuu*gCur+nub*(1-gCur) #Again not efficient here

#Setup boundary conditions to each face
def setBound(UIniBFX,UIniBFY,nVX=nVX,nVY=nVY,lenBF=lenBF,xRefine=xRefine,yRefine=yRefine):
    ##X sweep
    UBFX = np.zeros([nVX,nVY-1,lenBF]) #Temperature infomation for X facing faces. Third element is for mix boundary cond
    for indX in range(0,nVX):
        for indY in range(0,nVY-1): #For x direction, we have nVX face in X direction but only nVY-1 rows in Y direction
            xCoor = xRefine[indX][indY]
            yCoor = yRefine[indX][indY]
            for j in range(0, len(UIniBFX)):
                if (xCoor>=UIniBFX[j][0][0]) & (xCoor<=UIniBFX[j][0][1]) & (yCoor>=UIniBFX[j][1][0]) & (yCoor<=UIniBFX[j][1][1]):
                    UBFX[indX][indY][0] = UIniBFX[j][2]
                    UBFX[indX][indY][1:1+len(UIniBFX[j][3])] = np.array(UIniBFX[j][3])

    ##Y sweep For each face: First element is type following elements are the corresponding values T or q or Tinf hinf
    UBFY = np.zeros([nVX-1,nVY,lenBF])#Temperature infomation for Y facing faces. Third element is for mix boundary cond
    for indX in range(0,nVX-1):
        for indY in range(0,nVY):
            xCoor = xRefine[indX][indY]
            yCoor = yRefine[indX][indY]
            for j in range(0, len(UIniBFY)):
                if (xCoor>=UIniBFY[j][0][0]) & (xCoor<=UIniBFY[j][0][1]) & (yCoor>=UIniBFY[j][1][0]) & (yCoor<=UIniBFY[j][1][1]):
                    UBFY[indX][indY][0] = UIniBFY[j][2]
                    UBFY[indX][indY][1:1+len(UIniBFY[j][3])] = np.array(UIniBFY[j][3])
    return UBFX,UBFY
UBFX,UBFY = setBound(UIniBFX,UIniBFY)
VBFX,VBFY = setBound(VIniBFX,VIniBFY)
PBFX,PBFY = setBound(PIniBFX,PIniBFY)

#Construct the matrix Ax = b
def absd(velo):
    if velo>0:
        return 1.
    elif velo == 0.:
        return 0.5
    else:
        return 0.

##Assemble diffusion matrix
ADiff = np.zeros([nCel*2,nCel*2]) #This is for first [U,V] LHS
bDiff = np.zeros(nCel*2) #This is for first [U,V] LHS  ADiff[U,V] + bDiff = sum(nu*grad(U)*S)
for i in range(0,nCel):#Go through each of the cell
    indX = mapCel(i)[0]
    indY = mapCel(i)[1]

    indXNei = [indX,indX+1,indX,indX-1]
    indYNei = [indY-1,indY,indY+1,indY]
    indFac  = ['s','e','n','w']#South East North West
    signFac = [-1 , 1 , 1 ,-1 ]
    indFacX = [indX,indX+1,indX,indX]
    indFacY = [indY,indY,indY+1,indY]
    for j in range(0,4):
        if (indFac[j] == 's') | (indFac[j] == 'n'):
            UBFCur = UBFY[indFacX[j]][indFacY[j]]
            VBFCur = VBFY[indFacX[j]][indFacY[j]]
            nuCur  =  nuY[indFacX[j]][indFacY[j]]
            areaCur = areaY[indFacX[j]]
            disCellCur = disCellY[indFacY[j]]
        else:
            UBFCur = UBFX[indFacX[j]][indFacY[j]]
            VBFCur = VBFX[indFacX[j]][indFacY[j]]
            nuCur   = nuX[indFacX[j]][indFacY[j]]
            areaCur = areaX[indFacY[j]]
            disCellCur = disCellX[indFacX[j]]
        
        ## Construct matrix
        ### X Direction
        if UBFCur[0] == 1: #Dirichlet
            NAD = nuCur*areaCur/disCellCur
            bDiff[i] += NAD*UBFCur[1]
            ADiff[i][i] += -NAD
        elif UBFCur[0] == 2: #Neumman
            bDiff[i] += signFac[j]*nuCur*areaCur*UBFCur[1]
        else: #Internal
            NAD = nuCur*areaCur/disCellCur
            ADiff[i][reMapCel(indXNei[j],indYNei[j])] += NAD
            ADiff[i][i] += -NAD
        ### Y Direction
        if VBFCur[0] == 1: #Dirichlet
            NAD = nuCur*areaCur/disCellCur
            bDiff[i+nCel] += NAD*VBFCur[1]
            ADiff[i+nCel][i+nCel] += -NAD
        elif VBFCur[0] == 2: #Neumman
            bDiff[i+nCel] += signFac[j]*nuCur*areaCur*VBFCur[1]
        else: #Internal
            NAD = nuCur*areaCur/disCellCur
            ADiff[i+nCel][reMapCel(indXNei[j],indYNei[j])+nCel] += NAD
            ADiff[i+nCel][i+nCel] += -NAD

##Assemble pressure matrix
APres = np.zeros([nCel*2,nCel])
bPres = np.zeros(nCel*2) #APres[P] + bPres = -grad(P)*Volume
for i in range(0,nCel):#Go through each of the cell
    indX = mapCel(i)[0]
    indY = mapCel(i)[1]

    indFac  = ['s','e','n','w']#South East North West
    indFacX = [indX,indX+1,indX,indX]
    indFacY = [indY,indY,indY+1,indY]
    dXCur = areaY[indX]
    dYCur = areaX[indY]
    volCur = voluCell[i]

    for j in range(0,4):
        if (indFac[j] == 's') | (indFac[j] == 'n'):
            PBFCur = PBFY[indFacX[j]][indFacY[j]]
            disCellCur = disCellY[indFacY[j]]
        else:
            PBFCur = PBFX[indFacX[j]][indFacY[j]]
            disCellCur = disCellX[indFacX[j]]


        if indFac[j] == 's':
            multiCur = volCur/dYCur
            if PBFCur[0] == 1: #Dirichlet
                bPres[i+nCel] += multiCur*PBFCur[1]
            elif PBFCur[0] == 2: #Neumman
                bPres[i+nCel] += -PBFCur[1]*disCellCur*multiCur
                APres[i+nCel][i] += multiCur
            else: #Internal
                interCur = intCoefKY[indY]
                APres[i+nCel][i] += interCur*multiCur
                APres[i+nCel][reMapCel(indX,indY-1)] += (1-interCur)*multiCur
        elif indFac[j] == 'e':
            multiCur = -volCur/dXCur
            if PBFCur[0] == 1: #Dirichlet
                bPres[i] += multiCur*PBFCur[1]
            elif PBFCur[0] == 2: #Neumman
                bPres[i] += PBFCur[1]*disCellCur*multiCur
                APres[i][i] += multiCur
            else: #Internal
                interCur = intCoefKX[indX+1]
                APres[i][i] += (1-interCur)*multiCur
                APres[i][reMapCel(indX+1,indY)] += interCur*multiCur
        elif indFac[j] == 'n':
            multiCur = -volCur/dYCur
            if PBFCur[0] == 1: #Dirichlet
                bPres[i+nCel] += multiCur*PBFCur[1]
            elif PBFCur[0] == 2: #Neumman
                bPres[i+nCel] += PBFCur[1]*disCellCur*multiCur
                APres[i+nCel][i] += multiCur
            else: #Internal
                interCur = intCoefKY[indY+1]
                APres[i+nCel][i] += (1-interCur)*multiCur
                APres[i+nCel][reMapCel(indX,indY+1)] += interCur*multiCur
        else:
            multiCur = volCur/dXCur
            if PBFCur[0] == 1: #Dirichlet
                bPres[i] += PBFCur[1]*multiCur
            elif PBFCur[0] == 2: #Neumman
                bPres[i] += -PBFCur[1]*disCellCur*multiCur
                APres[i][i] += multiCur
            else: #Internal
                interCur = intCoefKX[indX]
                APres[i][i] += interCur*multiCur
                APres[i][reMapCel(indX-1,indY)] += (1-interCur)*multiCur

##Assemble advection matrix & Assemble sumPhi matrix
def calAdvMat(U,V,UBFX=UBFX,UBFY=UBFY,VBFX=VBFX,VBFY=VBFY,areaX=areaX,areaY=areaY,disCellX=disCellX,disCellY=disCellY,nCel=nCel,mapCel=mapCel,reMapCel=reMapCel,absd=absd):
    
    AAdv = np.zeros([nCel*2,nCel*2]) #This is for first [U,V] LHS
    bAdv = np.zeros(nCel*2) #This is for first [U,V] LHS  AAdv[U,V] + bAdv = sum(U*phi)
    bSumPhi = np.zeros(nCel) #This sum the explicit LHS   ASumPhi*[U,V] + bSumPhi
    ASumPhi = np.zeros([nCel,nCel*2]) #This sum the flux across all cells LHS
    for i in range(0,nCel):#Go through each of the cell
        indX = mapCel(i)[0]
        indY = mapCel(i)[1]

        indFac  = ['s','e','n','w']#South East North West
        indFacX = [indX,indX+1,indX,indX]
        indFacY = [indY,indY,indY+1,indY]
        for j in range(0,4):
            if (indFac[j] == 's') | (indFac[j] == 'n'):
                UBFCur = UBFY[indFacX[j]][indFacY[j]]
                VBFCur = VBFY[indFacX[j]][indFacY[j]]
                areaCur = areaY[indFacX[j]]
                disCellCur = disCellY[indFacY[j]]
            else:
                UBFCur = UBFX[indFacX[j]][indFacY[j]]
                VBFCur = VBFX[indFacX[j]][indFacY[j]]
                areaCur = areaX[indFacY[j]]
                disCellCur = disCellX[indFacX[j]]
            
            if indFac[j] == 's':
                if VBFCur[0] == 1: #Dirichlet
                    velAdv = VBFCur[1]
                    bSumPhi[i] += -VBFCur[1]*areaCur
                elif VBFCur[0] == 2: #Neumman
                    velAdv = -VBFCur[1]*disCellCur + V[reMapCel(indX,indY)]
                    bSumPhi[i] += -areaCur*(-VBFCur[1]*disCellCur)
                else: #Internal
                    velNei = V[reMapCel(indX,indY-1)]
                    velCen = V[reMapCel(indX,indY)]
                    velAdv =  absd(velNei)*velNei + absd(-velNei)*velCen
                phiAdv = -velAdv*areaCur
            elif indFac[j] == 'e':
                if UBFCur[0] == 1: #Dirichlet
                    velAdv = UBFCur[1]
                    bSumPhi[i] += areaCur*UBFCur[1]
                elif UBFCur[0] == 2: #Neumman
                    velAdv = UBFCur[1]*disCellCur + U[reMapCel(indX,indY)]
                    bSumPhi[i] += areaCur*UBFCur[1]*disCellCur
                else: #Internal
                    velNei = U[reMapCel(indX+1,indY)]
                    velCen = U[reMapCel(indX,indY)]
                    velAdv = absd(velCen)*velCen + absd(-velCen)*velNei
                phiAdv = velAdv*areaCur
            elif indFac[j] == 'n':
                if VBFCur[0] == 1: #Dirichlet
                    velAdv = VBFCur[1]
                    bSumPhi[i] += areaCur*VBFCur[1]
                elif VBFCur[0] == 2: #Neumman
                    velAdv = VBFCur[1]*disCellCur + V[reMapCel(indX,indY)]
                    bSumPhi[i] += areaCur*VBFCur[1]*disCellCur
                else: #Internal
                    velNei = V[reMapCel(indX,indY+1)]
                    velCen = V[reMapCel(indX,indY)]
                    velAdv = absd(velCen)*velCen + absd(-velCen)*velNei
                phiAdv = velAdv*areaCur
            else: #West face
                if UBFCur[0] == 1: #Dirichlet
                    velAdv = UBFCur[1]
                    bSumPhi[i] += -areaCur*UBFCur[1]
                elif UBFCur[0] == 2: #Neumman
                    velAdv = -UBFCur[1]*disCellCur + U[reMapCel(indX,indY)]
                    bSumPhi[i] += -areaCur*(-UBFCur[1]*disCellCur)
                else: #Internal
                    velNei = U[reMapCel(indX-1,indY)]
                    velCen = U[reMapCel(indX,indY)]
                    velAdv = absd(velNei)*velNei + absd(-velNei)*velCen
                phiAdv = -velAdv*areaCur
            
            ## Construct Matrix
            if indFac[j] == 's':
                if VBFCur[0] == 1: #Dirichlet
                    bAdv[i+nCel] += VBFCur[1]*phiAdv
                elif VBFCur[0] == 2: #Neumman
                    bAdv[i+nCel] += -VBFCur[1]*disCellCur*phiAdv
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv
                else: #Internal
                    AAdv[i+nCel][reMapCel(indX,indY-1)+nCel] +=  absd(V[reMapCel(indX,indY-1)])*phiAdv
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += absd(-V[reMapCel(indX,indY-1)])*phiAdv
                if UBFCur[0] == 1: #Dirichlet
                    bAdv[i] += UBFCur[1]*phiAdv
                elif UBFCur[0] == 2: #Neumman
                    bAdv[i] += -UBFCur[1]*disCellCur*phiAdv
                    AAdv[i][reMapCel(indX,indY)] += phiAdv
                else: #Internal
                    AAdv[i][reMapCel(indX,indY-1)] += absd(V[reMapCel(indX,indY-1)])*phiAdv
                    AAdv[i][reMapCel(indX,indY)] += absd(-V[reMapCel(indX,indY-1)])*phiAdv
            elif indFac[j] == 'e':
                if VBFCur[0] == 1: #Dirichlet
                    bAdv[i+nCel] += phiAdv*VBFCur[1]
                elif VBFCur[0] == 2: #Neumman
                    bAdv[i+nCel] += phiAdv*VBFCur[1]*disCellCur
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv
                else: #Internal
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv*absd(U[reMapCel(indX,indY)])
                    AAdv[i+nCel][reMapCel(indX+1,indY)+nCel] += phiAdv*absd(-U[reMapCel(indX,indY)])
                if UBFCur[0] == 1: #Dirichlet
                    bAdv[i] += phiAdv*UBFCur[1]
                elif UBFCur[0] == 2: #Neumman
                    bAdv[i] += phiAdv*UBFCur[1]*disCellCur
                    AAdv[i][reMapCel(indX,indY)] += phiAdv
                else: #Internal
                    AAdv[i][reMapCel(indX,indY)] += phiAdv*absd(U[reMapCel(indX,indY)])
                    AAdv[i][reMapCel(indX+1,indY)] += phiAdv*absd(-U[reMapCel(indX,indY)])
            elif indFac[j] == 'n':
                if VBFCur[0] == 1: #Dirichlet
                    bAdv[i+nCel] += phiAdv*VBFCur[1]
                elif VBFCur[0] == 2: #Neumman
                    bAdv[i+nCel] += phiAdv*VBFCur[1]*disCellCur
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv
                else: #Internal
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv*absd(V[reMapCel(indX,indY)])
                    AAdv[i+nCel][reMapCel(indX,indY+1)+nCel] += phiAdv*absd(-V[reMapCel(indX,indY)])
                if UBFCur[0] == 1: #Dirichlet
                    bAdv[i] += phiAdv*UBFCur[1]
                elif UBFCur[0] == 2: #Neumman
                    bAdv[i] += phiAdv*UBFCur[1]*disCellCur
                    AAdv[i][reMapCel(indX,indY)] += phiAdv
                else: #Internal
                    AAdv[i][reMapCel(indX,indY)] += phiAdv*absd(V[reMapCel(indX,indY)])
                    AAdv[i][reMapCel(indX,indY+1)] += phiAdv*absd(-V[reMapCel(indX,indY)])
            else: #West face
                if VBFCur[0] == 1: #Dirichlet
                    bAdv[i+nCel] += phiAdv*VBFCur[1]
                elif VBFCur[0] == 2: #Neumman
                    bAdv[i+nCel] += phiAdv*(-VBFCur[1]*disCellCur)
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv
                else: #Internal
                    AAdv[i+nCel][reMapCel(indX,indY)+nCel] += phiAdv*absd(-U[reMapCel(indX-1,indY)])
                    AAdv[i+nCel][reMapCel(indX-1,indY)+nCel] += phiAdv*absd(U[reMapCel(indX-1,indY)])
                if UBFCur[0] == 1: #Dirichlet
                    bAdv[i] += phiAdv*UBFCur[1]
                elif UBFCur[0] == 2: #Neumman
                    bAdv[i] += phiAdv*(-UBFCur[1]*disCellCur)
                    AAdv[i][reMapCel(indX,indY)] += phiAdv
                else: #Internal
                    AAdv[i][reMapCel(indX,indY)] += phiAdv*absd(-U[reMapCel(indX-1,indY)])
                    AAdv[i][reMapCel(indX-1,indY)] += phiAdv*absd(U[reMapCel(indX-1,indY)])

    for i in range(0,nCel):#Go through each of the cell
        indX = mapCel(i)[0]
        indY = mapCel(i)[1]

        indFac  = ['s','e','n','w']#South East North West
        indFacX = [indX,indX+1,indX,indX]
        indFacY = [indY,indY,indY+1,indY]
        for j in range(0,4):
            if (indFac[j] == 's') | (indFac[j] == 'n'):
                UBFCur = UBFY[indFacX[j]][indFacY[j]]
                VBFCur = VBFY[indFacX[j]][indFacY[j]]
                areaCur = areaY[indFacX[j]]
            else:
                UBFCur = UBFX[indFacX[j]][indFacY[j]]
                VBFCur = VBFX[indFacX[j]][indFacY[j]]
                areaCur = areaX[indFacY[j]]
            
            if indFac[j] == 's':
                if VBFCur[0] == 2: #Neumman
                    ASumPhi[i][reMapCel(indX,indY)+nCel] += -areaCur
                elif VBFCur[0] == 0: #Internal
                    velNei = V[reMapCel(indX,indY-1)]
                    ASumPhi[i][reMapCel(indX,indY)+nCel] += -areaCur*absd(-velNei)
                    ASumPhi[i][reMapCel(indX,indY-1)+nCel] += -areaCur*absd(velNei)
            elif indFac[j] == 'e':
                if UBFCur[0] == 2: #Neumman
                    ASumPhi[i][reMapCel(indX,indY)] += areaCur
                elif UBFCur[0] == 0: #Internal
                    velCen = U[reMapCel(indX,indY)]
                    ASumPhi[i][reMapCel(indX,indY)] += areaCur*absd(velCen)
                    ASumPhi[i][reMapCel(indX+1,indY)] += areaCur*absd(-velCen)
            elif indFac[j] == 'n':
                if VBFCur[0] == 2: #Neumman
                    ASumPhi[i][reMapCel(indX,indY)+nCel] += areaCur
                elif VBFCur[0] == 0: #Internal
                    velCen = V[reMapCel(indX,indY)]
                    ASumPhi[i][reMapCel(indX,indY)+nCel] += areaCur*absd(velCen)
                    ASumPhi[i][reMapCel(indX,indY+1)+nCel] += areaCur*absd(-velCen)
            else: #West face
                if UBFCur[0] == 2: #Neumman
                    ASumPhi[i][reMapCel(indX,indY)] += -areaCur
                elif UBFCur[0] == 0: #Internal
                    velNei = U[reMapCel(indX-1,indY)]
                    ASumPhi[i][reMapCel(indX,indY)] += -areaCur*absd(-velNei)
                    ASumPhi[i][reMapCel(indX-1,indY)] += -areaCur*absd(velNei)
    return AAdv,bAdv,ASumPhi,bSumPhi

##Preprocess matrices
def preProcess(AAdv,bAdv,ADiff=ADiff,bDiff=bDiff,APres=APres,bPres=bPres):
    #### Derivation:
    #### AAdv*velVec + bAdv = APres*preVec + bPres + ADiff*velVec + bDiff
    #### (AAdv-ADiff)*velVec = APres*preVec + bPres + bDiff - bAdv
    #### ATot*velVec = APres*preVec + bTot
    #### ATot*velVec = APTot*preVec + bTot
    #### Derivation
    ATot   = AAdv-ADiff
    bTot   = bPres + bDiff - bAdv
    APTot  = APres*1.
    ##### Normalization
    for i in range(0,len(ATot)):
        divider  = ATot[i][i]
        bTot[i]  = bTot[i]/divider
        APTot[i] = APTot[i]/divider
        ATot[i]  = ATot[i]/divider
    #### ATot*velVec = APTot*preVec + bTot
    HTot = ATot*1.
    for i in range(0,len(HTot)):
        HTot[i][i] = 0.0
    #### velVec + HTot*velVec = APTot*preVec + bTot
    return ATot,bTot,APTot,HTot

##Steep Descent Solver
def SDSolve(ASD,bSD,guessIni=None,numSD=numSD):
    if guessIni is None:
        guessIni = np.zeros(len(ASD))
    #Iterate
    for i in range(0,numSD):
        resid = bSD - np.matmul(ASD,guessIni)
        alpha = np.dot(resid,resid)/np.dot(resid,np.matmul(ASD,resid))
        guessIni = guessIni+alpha*resid
        if np.mean(abs(alpha*resid))/np.mean(abs(guessIni)) <0.0000001:
            return guessIni
    return guessIni


##Solve the matrix
USave = []
VSave = []
preSave = []
conSave = 0 #Counter to save data
velVec = np.array(list(U)+list(V))
preVec = P*1.
for i in range(0,numOut):
    conSave += 1
    velOldOut = velVec*1. #Save old time velocity for outer loop relaxation
    ###Update advection terms for current iteration and will be fixed
    AAdv,bAdv,ASumPhi,bSumPhi = calAdvMat(velVec[0:len(U)],velVec[len(U):])
    ###Preprocess matrices
    ATot,bTot,APTot,HTot = preProcess(AAdv,bAdv)
    ###Predictor Step
    velVec = SDSolve(ATot , np.matmul(APTot,preVec)+bTot , guessIni=velVec) #v* from here
    ###Relaxation Step
    velVec = facRelOut*velVec + (1.-facRelOut)*velOldOut
    ###Corrector Loop
    AAP = np.matmul(ASumPhi,APTot)
    AH = np.matmul(ASumPhi,HTot)
    AB = np.matmul(ASumPhi,bTot)
    BS = bSumPhi

    #### Normalization
    for k in range(0,len(AAP)):
        divider  = AAP[k][k]
        if divider == 0:
            print('Divided by 0 warning')
        AB[k]  = AB[k]/divider
        AH[k]  = AH[k]/divider
        AAP[k] = AAP[k]/divider
        BS[k]  = BS[k]/divider

    for j in range(0,numIn):
        velOldIn = velVec*1. #Save old time velocity for inner loop relaxation
        preOldIn = preVec*1. #Save old time pressure for inner loop relaxation
        ### PPE to update new pressure
        preVec = SDSolve(AAP, np.matmul(AH,velVec)-AB-BS , guessIni=preVec) #P n+1 from V*
        ### Correct velocity
        velVec = -np.matmul(HTot,velVec) + np.matmul(APTot,preVec) + bTot #vel n+1 from P n+1
        ### Relaxation step
        velVec = facRelIn*velVec + (1.-facRelIn)*velOldIn
        preVec = facRelIn*preVec + (1.-facRelIn)*preOldIn
    
    #### Simple Progress Viewer
    print('Current Iteration: ' + str(i) + '\n')
    print('Pressure:\n')
    print(preVec)
    print('Velocity:\n')
    print(velVec)
    print('\n')

    #### Save Data
    if conSave == intSav:
        conSave = 0
        USave.append(velVec[0:len(U)]*1.)
        VSave.append(velVec[len(U):]*1.)
        preSave.append(preVec*1.)

##Print out result
coorSave = np.zeros([3,nCel])
for i in range(0,nCel):
    coorSave[0][i] = i
    coorSave[1][i] = centCellX[mapCel(i)[0]]
    coorSave[2][i] = centCellY[mapCel(i)[1]]

import csv
file = open('USave.csv','w')
writer = csv.writer(file)
for i in range(0,len(USave)):
    writer.writerow(USave[i])
file.close()

file = open('VSave.csv','w')
writer = csv.writer(file)
for i in range(0,len(VSave)):
    writer.writerow(VSave[i])
file.close()

file = open('PSave.csv','w')
writer = csv.writer(file)
for i in range(0,len(preSave)):
    writer.writerow(preSave[i])
file.close()

file = open('coSave.csv','w')
writer = csv.writer(file)
for i in range(0,len(coorSave)):
    writer.writerow(coorSave[i])
file.close()

print("end")

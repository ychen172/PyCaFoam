#MixingLayerTestCase########################################
############################################################################ User Input Section
#Domain Parameters
##User specify the mesh refinment in x and y direction
xRefine = list(np.arange(0.,0.031,0.001)) #Vertices along X direction
yRefine = list(np.arange(0.,0.0105,0.0005)) #Vertices along Y direction

#Boundary Parameters
##         [[xZoneLowBound, xZoneUpBound],[yZoneLowBound, yZoneUpBound],Boundary Type,[Boundary Values]]
##         Boundary Type 0 --> Internal Face (don't need to set this, default automatically)
##                       1 --> Dirichlet Face     Boundary Values --> The value fixed to (positive in +x and +y direction regardless inlet or outlet)
##                       2 --> Neumman Face       Boundary Values --> The gradient fixed to (positive in +x and +y direction regardless inlet or outlet)
##         Faces not covered will be automatically internal faces
UIniBFX = [[[0,0],[0.0066,0.0100],1,[1.5]],
           [[0,0],[0.0033,0.0066],1,[2.]],
           [[0,0],[0.0000,0.0033],1,[1.5]],
           [[0.03,0.03],[0,0.01],2,[0.]]] #X direction face boundary conditions
VIniBFX = [[[0,0],[0,0.01],1,[0.]],
           [[0.03,0.03],[0,0.01],2,[0.]]] #X direction face boundary conditions
PIniBFX = [[[0,0],[0.0066,0.0100],2,[0.]],
           [[0,0],[0.0033,0.0066],2,[0.]],
           [[0,0],[0.0000,0.0033],2,[0.]],
           [[0.03,0.03],[0,0.01],1,[78364.]]] #X direction face boundary conditions
UIniBFY = [[[0,0.03],[0,0],2,[0.]],
           [[0,0.03],[0.01,0.01],2,[0.]]]
VIniBFY = [[[0,0.03],[0,0],1,[0.]],
           [[0,0.03],[0.01,0.01],1,[0.]]]
PIniBFY = [[[0,0.03],[0,0],2,[0.]],
           [[0,0.03],[0.01,0.01],2,[0.]]]

#Initial Parameters
##      [ [[xZone1LowBound, xZone1UpBound],[yZone1LowBound, yZone1UpBound], valueForZone1]  ,  [[xZone2LowBound, xZone2UpBound],[yZone2LowBound, yZone2UpBound], valueForZone2]  , valueForAllOtherZone]
##      valueForAllOtherZone is always required as the last element, all other zones definition are optional
UIni  = [2]
VIni  = [0]
PIni  = [78364.]
nuIni  = [1.5e-1] 

#Numerical Parameters
numOut = 500    #Outer loop iteration
numIn  = 2    #Inner loop iteration
intSav = 50     #Save interval (iteration)
facRelOut = 0.005 #Relaxation factor for momentum prediction step [0:No progress] [1:Full Progress]
facRelIn = 0.005  #Relaxation factor for pressure correction step [0:No progress] [1:Full Progress]
###########################################################################

#Make nonuniform###########################################################
for i in range(1,len(xRefine)-1,2):
    xRefine[i] += 0.0005
for i in range(1,len(yRefine)-1,2):
    yRefine[i] += 0.00025
xRefine = list(xRefine)
yRefine = list(yRefine)
###########################################################################
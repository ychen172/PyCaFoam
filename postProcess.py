import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
U    = pd.read_csv('USave.csv')
V    = pd.read_csv('VSave.csv')
P    = pd.read_csv('PSave.csv')
Co   = pd.read_csv('coSave.csv')
U = U.values
V = V.values
P = P.values
Co = Co.values
X = Co[0]
Y = Co[1]
XU = np.unique(X)
YU = np.unique(Y)

k = 0
t = -1
UMat = np.zeros([len(XU),len(YU)])
VMat = np.zeros([len(XU),len(YU)])
PMat = np.zeros([len(XU),len(YU)])
for j in range(0,len(YU)):
    for i in range(0,len(XU)):
        PMat[i][j] = P[t][k]
        VMat[i][j] = V[t][k]
        UMat[i][j] = U[t][k]
        k += 1



fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(0,len(XU)):
    plt.plot(UMat[i],YU, label="X="+str(np.round(XU[i],3))+"m")
plt.ylabel('Y[m]')
plt.xlabel('U[m/s]')
ax.set_yticks(YU)
ax.set_xticks(np.arange(0,4.2,0.2))
plt.grid(True)
plt.legend(loc="upper right")
plt.savefig(os.getcwd() + '/' + 'Movie' + '/' + "Uvel.jpg")
plt.show()

VMatT = np.transpose(VMat)

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(0,len(YU)):
    plt.plot(XU,VMatT[i], label="Y="+str(np.round(YU[i],3))+"m")
plt.xlabel('X[m]')
plt.ylabel('V[m/s]')
ax.set_xticks(XU)
plt.ylim([-1,2])
plt.grid(True)
plt.legend(loc="upper right")
plt.savefig(os.getcwd() + '/' + 'Movie' + '/' + "Vvel.jpg")
plt.show()

[YMesh,XMesh] = np.meshgrid(YU,XU)
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
Color = plt.contourf(XMesh,YMesh,PMat,40)
cbar = fig.colorbar(Color)
cbar.ax.set_ylabel('   P/rho [m^2/s^2]', rotation=270)
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
ax.set_xticks(XU)
ax.set_yticks(YU)
plt.grid(True)
plt.savefig(os.getcwd() + '/' + 'Movie' + '/' + "P.jpg")
plt.show()

print("end")
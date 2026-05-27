#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1. 读入数据
data = np.loadtxt("mhd2D_results.dat")

# 假设文件的前两列是 x, y，第三列是 rho，第四列是 p
x   = data[:, 0]
y   = data[:, 1]
rho = data[:, 2]
p_  = data[:, 3]

# 2. 获取网格尺寸
x_unique = np.unique(x)
y_unique = np.unique(y)
Nx = len(x_unique)
Ny = len(y_unique)

# 3. 将一维数据重塑为二维网格
#   需要根据实际的输出循环顺序来决定 Nx, Ny 的位置。
#   如果外层循环是 x，内层循环是 y，则可以先尝试：
#       X   = x.reshape(Nx, Ny)
#       Y   = y.reshape(Nx, Ny)
#       RHO = rho.reshape(Nx, Ny)
#       P   = p_.reshape(Nx, Ny)
#   如果画图发现坐标对应混乱，可尝试转置或者将 Nx, Ny 调换。
try:
    X   = x.reshape(Nx, Ny)
    Y   = y.reshape(Nx, Ny)
    RHO = rho.reshape(Nx, Ny)
    P   = p_.reshape(Nx, Ny)
except ValueError:
    # 如果上面的 reshape 失败或图像坐标不正确，可以尝试调换 Nx 和 Ny
    X   = x.reshape(Ny, Nx)
    Y   = y.reshape(Ny, Nx)
    RHO = rho.reshape(Ny, Nx)
    P   = p_.reshape(Ny, Nx)

# 4. 绘制 3D 曲面图
fig = plt.figure(figsize=(12, 5))

# (1) 密度分布
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf1 = ax1.plot_surface(X, Y, RHO, cmap='viridis', edgecolor='none')
ax1.set_title('Density distribution')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('rho')
fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)

# (2) 压力分布
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax2.plot_surface(X, Y, P, cmap='viridis', edgecolor='none')
ax2.set_title('Pressure distribution')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('p')
fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)

plt.tight_layout()
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1. 读入数据
data = np.loadtxt("mhd2D_results.dat")

# 假设文件的前两列是 x, y，第三列是 rho，第四列是 p
x   = data[:, 0]
y   = data[:, 1]
rho = data[:, 2]
p_  = data[:, 3]

# 2. 获取网格尺寸
x_unique = np.unique(x)
y_unique = np.unique(y)
Nx = len(x_unique)
Ny = len(y_unique)

# 3. 将一维数据重塑为二维网格
#   需要根据实际的输出循环顺序来决定 Nx, Ny 的位置。
#   如果外层循环是 x，内层循环是 y，则可以先尝试：
#       X   = x.reshape(Nx, Ny)
#       Y   = y.reshape(Nx, Ny)
#       RHO = rho.reshape(Nx, Ny)
#       P   = p_.reshape(Nx, Ny)
#   如果画图发现坐标对应混乱，可尝试转置或者将 Nx, Ny 调换。
try:
    X   = x.reshape(Nx, Ny)
    Y   = y.reshape(Nx, Ny)
    RHO = rho.reshape(Nx, Ny)
    P   = p_.reshape(Nx, Ny)
except ValueError:
    # 如果上面的 reshape 失败或图像坐标不正确，可以尝试调换 Nx 和 Ny
    X   = x.reshape(Ny, Nx)
    Y   = y.reshape(Ny, Nx)
    RHO = rho.reshape(Ny, Nx)
    P   = p_.reshape(Ny, Nx)

# 4. 绘制 3D 曲面图
fig = plt.figure(figsize=(12, 5))

# (1) 密度分布
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf1 = ax1.plot_surface(X, Y, RHO, cmap='viridis', edgecolor='none')
ax1.set_title('Density distribution')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('rho')
fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)

# (2) 压力分布
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax2.plot_surface(X, Y, P, cmap='viridis', edgecolor='none')
ax2.set_title('Pressure distribution')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('p')
fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)

plt.tight_layout()
plt.show()

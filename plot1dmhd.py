#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def main():
    # 1. 读入数据
    # 假设文件格式为:
    # x   rho   p   v_x   v_y   v_z   U
    data = np.loadtxt("mhd2D_results.dat")
    y   = data[:, 0]
    rho = data[:, 1]
    p   = data[:, 2]
    vx  = data[:, 3]
    vy  = data[:, 4]
    vz  = data[:, 5]
    p_rho = data[:, 6]
    
    # 2. 创建 2x3 的子图布局
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 8))
    fig.suptitle("1D Brio & Wu test in y-direction", fontsize=16)

    # 为了让曲线看起来更像激波管的离散分布，可以加一些散点或阶跃线
    # 下面示例中采用带散点的折线图，你也可以改成仅散点或仅折线
    
    # (1) x vs rho
    axes[0, 0].plot(y, rho, 'o-', ms=1.5, lw=1)
    axes[0, 0].set_xlabel("y")
    axes[0, 0].set_ylabel(r"$\rho$")

    # (2) x vs p
    axes[0, 1].plot(y, p, 'o-', ms=1.5, lw=1)
    axes[0, 1].set_xlabel("y")
    axes[0, 1].set_ylabel("p")

    # (3) x vs v_x
    axes[1, 0].plot(y, vx, 'o-', ms=1.5, lw=1)
    axes[1, 0].set_xlabel("y")
    axes[1, 0].set_ylabel(r"$v_x$")

    # (4) x vs v_y
    axes[1, 1].plot(y, vy, 'o-', ms=1.5, lw=1)
    axes[1, 1].set_xlabel("y")
    axes[1, 1].set_ylabel(r"$v_y$")

    # (5) x vs v_z
    axes[2, 0].plot(y, vz, 'o-', ms=1.5, lw=1)
    axes[2, 0].set_xlabel("y")
    axes[2, 0].set_ylabel(r"$B_x$")

    # (6) x vs U
    axes[2, 1].plot(y, p_rho, 'o-', ms=1.5, lw=1)
    axes[2, 1].set_xlabel("y")
    axes[2, 1].set_ylabel(r"$\frac{p}{\rho}$")


    # 3. 布局与显示
    plt.tight_layout()
    # 让标题和子图不重叠
    plt.subplots_adjust(top=0.92)

    # 保存图片
    plt.savefig("1D_MHD_plot.png", dpi=300)

    # 显示图像
    plt.show()

if __name__ == "__main__":
    main()

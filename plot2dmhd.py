import numpy as np
import matplotlib.pyplot as plt

def main():
    # 1. 读取数据文件
    # 假设文件每行格式为： x  y  rho
    # 如果有其他列，请自行修改对应的索引
    data = np.loadtxt("DCmhd2D_results.dat")
    x   = data[:, 0]
    y   = data[:, 1]
    rho = data[:, 2]

    # 2. 获取网格信息
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx = len(x_unique)
    ny = len(y_unique)

    # 3. 将一维数组重塑为二维 (nx, ny)
    #    注意：一定要确认文件输出顺序是“先固定 x，再扫 y”
    #    如果实际是“先固定 y，再扫 x”，则需要改成 (ny, nx) 并调整 meshgrid
    rho_grid = rho.reshape((nx, ny))

    # 4. 构造网格
    X, Y = np.meshgrid(x_unique, y_unique, indexing='ij')

    # 5. 绘制等值线图
    plt.figure(figsize=(6, 6))  # 可根据需要调整画布大小

    # 使用自动范围，不设置 vmin、vmax
    # levels=40 表示等值线数量，你也可以改成其他值或直接写 contourf(X, Y, rho_grid, cmap='bone')
    contour = plt.contourf(X, Y, rho_grid, levels=40, cmap='bone')

    # 6. 添加颜色条
    cbar = plt.colorbar(contour)
    cbar.set_label(r"Density $\rho$")

    # 7. 标题、坐标轴等
    plt.title(r"Density $\rho$ at $t=0.50008\mathrm{s}$")
    plt.xlabel("x")
    plt.ylabel("y")

    # 8. 调整布局并显示
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

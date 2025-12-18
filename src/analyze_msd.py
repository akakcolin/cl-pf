#!/usr/bin/env python3
"""
MSD数据分析脚本
- 读取GROMACS xvg文件
- 绘制MSD曲线
- 线性拟合计算扩散系数
- 绘制双对数曲线
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 读取xvg文件
def read_xvg(filename):
    time, msd_li, msd_pf6 = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                time.append(float(parts[0]))
                msd_li.append(float(parts[1]))
                msd_pf6.append(float(parts[2]))
    return np.array(time), np.array(msd_li), np.array(msd_pf6)

# 线性函数用于拟合
def linear(t, D, b):
    return 6 * D * t + b

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(
        description='MSD数据分析: 读取xvg文件，计算扩散系数，绘制MSD曲线',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s msd.xvg
  %(prog)s msd.xvg -s 5000 -e 10000
  %(prog)s msd.xvg -o output.png -t "25°C"
  %(prog)s msd.xvg --labels "Li+" "TFSI-"
        '''
    )
    parser.add_argument('input', help='输入的xvg文件')
    parser.add_argument('-s', '--start', type=float, default=5000,
                        help='拟合起始时间 (ps)，默认: 5000')
    parser.add_argument('-e', '--end', type=float, default=10000,
                        help='拟合结束时间 (ps)，默认: 10000')
    parser.add_argument('-o', '--output', default='msd_analysis.png',
                        help='输出图像文件名，默认: msd_analysis.png')
    parser.add_argument('-t', '--title', default=None,
                        help='图像标题后缀 (如温度)，默认自动从文件名提取')
    parser.add_argument('--labels', nargs=2, default=['Li$^+$', 'PF$_6^-$'],
                        metavar=('LABEL1', 'LABEL2'),
                        help='两个物种的标签，默认: "Li+" "PF6-"')
    parser.add_argument('--no-show', action='store_true',
                        help='不显示图像窗口，仅保存文件')

    args = parser.parse_args()

    # 读取数据
    time, msd_1, msd_2 = read_xvg(args.input)

    # 拟合区间
    fit_start, fit_end = args.start, args.end
    fit_mask = (time >= fit_start) & (time <= fit_end)
    time_fit = time[fit_mask]

    # 拟合第一个物种
    msd_1_fit = msd_1[fit_mask]
    popt_1, _ = curve_fit(linear, time_fit, msd_1_fit)
    print(popt_1)
    D_1 = popt_1[0] * 1e7 / 1e5*100  # nm^2/ps -> 1e-5 cm^2/s
    print(f"{args.labels[0]} 扩散系数: {D_1:.4f} x 1e-7 cm^2/s")

    # 拟合第二个物种
    msd_2_fit = msd_2[fit_mask]
    popt_2, b_2 = curve_fit(linear, time_fit, msd_2_fit)
    D_2 = popt_2[0] * 1e7 / 1e5*100  # nm^2/ps -> 1e-5 cm^2/s
    print(f"{args.labels[1]} 扩散系数: {D_2:.4f} x 1e-7 cm^2/s")

    # 确定标题
    if args.title:
        title_suffix = f" ({args.title})"
    else:
        title_suffix = ""

    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 图1: MSD曲线及拟合
    ax1 = axes[0]
    ax1.plot(time, msd_1, 'b-', label=args.labels[0], linewidth=1.5)
    ax1.plot(time, msd_2, 'r-', label=args.labels[1], linewidth=1.5)

    # 拟合线
    time_plot = np.linspace(fit_start, fit_end, 100)
    ax1.plot(time_plot, linear(time_plot, *popt_1), 'b--',
             label=f'{args.labels[0]} fit: D={D_1:.4f}e-7 cm$^2$/s', linewidth=2)
    ax1.plot(time_plot, linear(time_plot, *popt_2), 'r--',
             label=f'{args.labels[1]} fit: D={D_2:.4f}e-7 cm$^2/s$', linewidth=2)

    # 标记拟合区间
    ax1.axvline(fit_start, color='gray', linestyle=':', alpha=0.5)
    ax1.axvline(fit_end, color='gray', linestyle=':', alpha=0.5)

    ax1.set_xlabel('Time (ps)', fontsize=12)
    ax1.set_ylabel('MSD (nm$^2$)', fontsize=12)
    ax1.set_title(f'Mean Squared Displacement{title_suffix}', fontsize=14)
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3)

    # 图2: 双对数曲线
    ax2 = axes[1]
    # 排除t=0
    mask = time > 0
    ax2.loglog(time[fit_mask]/1000, msd_1[fit_mask], 'b-', label=args.labels[0], linewidth=1.5)
    ax2.loglog(time[fit_mask]/1000, msd_2[fit_mask], 'r-', label=args.labels[1], linewidth=1.5)

    # 参考斜率线
    t_ref = np.logspace(2, 4, 50)
    ax2.loglog(time_plot/1000, linear(time_plot, *popt_1), 'k--', lw=3,alpha=0.5, label='Li$^+$ fit')
    ax2.loglog(time_plot/1000, linear(time_plot, *popt_2), 'g--', lw=3,alpha=0.5, label='PF6- fit')

    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.set_ylabel('MSD (nm$^2$)', fontsize=12)
    ax2.set_title('Log-Log MSD Plot', fontsize=14)
    ax2.legend(loc='upper left', fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(args.output, dpi=500, bbox_inches='tight')

    if not args.no_show:
        plt.show()

    print(f"\n图像已保存为 {args.output}")

if __name__ == '__main__':
    main()

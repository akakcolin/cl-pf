 # Gnuplot script for RDF plotting
# 用法:
# 1. 基本用法（仅 RDF）: gnuplot plot_rdf.gp
# 2. 自定义范围和输出: gnuplot -e "rmax=6; out='rdf_6A'" plot_rdf.gp
# 3. 包含配位数: gnuplot -e "show_cn=1; rmax=8; out='rdf_with_cn'" plot_rdf.gp
#
# 数据文件要求:
# - *.dat 文件格式: index r g(r) CN(r)
# - 必须包含: Si-C.dat Si-Si.dat Si-N.dat C-C.dat C-N.dat N-N.dat
if (!exists("out")) out = "all_rdf"
if (!exists("rmax")) rmax = 8
if (!exists("show_cn")) show_cn = 0

set terminal pngcairo enhanced font "Arial,12" size 800,600
set output sprintf('%s.png', out)
# 可选：生成 PDF 版本（手动启用）
# set terminal pdfcairo enhanced font "Arial,12" size 6.4,4.8
# set output sprintf('%s.pdf', out)

# 设置图表外观
set grid
set border linewidth 1.5

# 设置坐标轴范围
set xrange [0:rmax]
set yrange [0:*]

# 设置坐标轴标签
set xlabel "径向距离 r (Å)" font "Arial,14"
set ylabel "径向分布函数 g(r)" font "Arial,14"
if (show_cn) set y2label "配位数 CN(r)" font "Arial,14"

# 标题与多轴设置
if (show_cn) {
  set y2tics nomirror
  set title sprintf("SiCN 体系的径向分布函数与配位数 (r≤%.1f Å)", rmax) font "Arial,16" offset 0,-1
} else {
  set title sprintf("SiCN 体系的径向分布函数 (r≤%.1f Å)", rmax) font "Arial,16" offset 0,-1
}

set key top right font "Arial,11"
set key box linestyle 1 linewidth 1

set style line 1 lc rgb '#FF4444' lw 2 pt 7 ps 0.5  # Si-C - 红色
set style line 2 lc rgb '#4444FF' lw 2 pt 7 ps 0.5  # Si-Si - 蓝色
set style line 3 lc rgb '#44FF44' lw 2 pt 7 ps 0.5  # Si-N - 绿色
set style line 4 lc rgb '#FFAA00' lw 2 pt 7 ps 0.5  # C-C - 橙色
set style line 5 lc rgb '#AA00FF' lw 2 pt 7 ps 0.5  # C-N - 紫色
set style line 6 lc rgb '#00AAAA' lw 2 pt 7 ps 0.5  # N-N - 青色


FILE_LIST = "Si-C Si-Si Si-N C-C C-N N-N"
# 可通过命令行覆盖：gnuplot -e "files='Si-C Si-N'" plot_rdf.gp
if (exists("files")) FILE_LIST = files
PLOT_CMD = ""

# 当显示配位数时启用 y2 轴范围
if (show_cn) {
  set y2range [0:*]
}

do for [i=1:words(FILE_LIST)] {
  name = word(FILE_LIST, i)
  PLOT_CMD = PLOT_CMD . (strlen(PLOT_CMD) > 0 ? ", " : "")
  extra = show_cn ? sprintf(", '%s.dat' using 2:4 with lines lw 1.5 lc rgb 'gray40' axes x1y2 notitle", name) : ""
  PLOT_CMD = PLOT_CMD . sprintf("'%s.dat' using 2:3 with lp ls %d title '%s'%s", name, i, name, extra)
}

# 执行绘图
eval("plot " . PLOT_CMD)

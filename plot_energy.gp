set terminal pngcairo size 600, 800 enhanced font 'Times,16'
set output 'energy.png'

set title "Energy vs Time: E(t)"
set xlabel "t"
set ylabel "E"

set grid

set key inside top right box font ",14"
#unset key
plot \
  "ke.txt" using 2:3 with lines lw 2 lc rgb "blue" title "CPU Code", \
  "ke_gpu.txt" using 2:3 with lines lw 2 lc rgb "red" dt 2  title "GPU Code"

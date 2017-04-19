set terminal postscript enhanced color
set output "position.eps"

p 'position.dat' u 1:2 w l not,\
'' u 1:3 w l not,\
'' u 1:4 w l not,\
'' u 1:5 w l not,\
'' u 1:6 w l not,\
'' u 1:7 w l not,\
'' u 1:8 w l not,\
'' u 1:9 w l not,\
'' u 1:10 w l not,\
'' u 1:11 w l not
#'' u 12:($10*0) w p pt 7 not,\
#'' u 13:($10*0) w p pt 7 not,\
#'' u 14:($10*0) w p pt 7 not,\
#'' u 15:($10*0) w p pt 7 not,\
#'' u 16:($10*0) w p pt 7 not,\
#'' u 17:($10*0) w p pt 7 not,\
#'' u 18:($10*0) w p pt 7 not,\
#'' u 19:($10*0) w p pt 7 not,\
3'' u 20:($10*0) w p pt 7 not,\

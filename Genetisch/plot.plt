#==================================================
#         Author : Thomas Rometsch
#      Plot name : data fit
#==================================================

reset
clear

#==================================================
#             Definitions
#==================================================

writepdf = 0;

file1 = "lichtkurve.dat"

label0 = "lichtkurve"

output_name = "fits.pdf"

#==================================================
#             Set terminal
#==================================================

if (writepdf!=1) {
  set term qt enhanced
}
if (writepdf==1) {
  set term pdfcairo enhanced size 7in,4in dashed
  set output output_name
}

#==================================================
#             Style
#==================================================

set style line 1 lt 7 ps 0.3 lc rgb "red" lw 1
set style line 2 lt 7 ps 0.3 lc rgb "orange" lw 1
set style line 3 lt 7 ps 0.3 lc rgb "#6495ED" lw 1
set style line 4 lt 4 lc rgb "black" lw 1
set style line 5 lt 1 lc rgb "forest-green" lw 1 ps 1

#set format y "%2.1tx10^{%L}"
#==================================================
#             Axis Labels
#==================================================

set ylabel "y"
set xlabel "x"


#==================================================
#             Analytic Functions
#==================================================
#Fitting P1
#status = 0
#x = 0.49999
#y = 0.49999

#Fitting P2
#status = 0
#x = 0.59986
#y = 0.10055

chisq_reduced_0 = 1279.22
sin0(x) = 1.0418*x + 17.186

chisq_reduced_1 = 771.724
sin1(x) = 1.004*x + 18.991 + 11.428*sin( 2*pi*( x/9.98023 + 0.99998) ) 

chisq_reduced_3 = 412.629
sin3(x) = 1.0215*x + 17.989 + 3.963*sin( 2*pi*( x/4.33151 + 0.99916) )  + 7.868*sin( 2*pi*( x/4.626 + 0.00016) )  + 9.998*sin( 2*pi*( x/10.0072 + 0.9984) ) 

chisq_reduced_5 = 445.777
sin5(x) = 1.0054*x + 18.509 + 4.961*sin( 2*pi*( x/4.33151 + 0.99368) )  + 9.981*sin( 2*pi*( x/9.8151 + 0.89756) )  + 12.066*sin( 2*pi*( x/29.8546 + 0.99706) )  + 7.247*sin( 2*pi*( x/4.6162 + 0.00258) )  + 11.486*sin( 2*pi*( x/30.449 + 0.01162) ) 



set samples 100000

#==================================================
#             Labels
#==================================================

label1 = "m=0 , {/Symbol c}^2 = ".sprintf("%f",chisq_reduced_0);
label2 = "m=1 , {/Symbol c}^2 = ".sprintf("%f",chisq_reduced_1);
label3 = "m=3 , {/Symbol c}^2 = ".sprintf("%f",chisq_reduced_3);
label4 = "m=5 , {/Symbol c}^2 = ".sprintf("%f",chisq_reduced_5);



#==================================================
#             Plot
#==================================================

set key top right
set multiplot layout 2,2

set title label1
plot  file1 u 1:2 notitle ls 1 w p, \
      sin0(x) ls 2 notitle

set title label2
plot  file1 u 1:2 notitle ls 1 w p, \
      sin1(x) ls 3 notitle

set title label3
plot  file1 u 1:2 notitle ls 1 w p, \
      sin3(x) ls 4 notitle

set title label4
plot  file1 u 1:2 notitle ls 1 w p, \
      sin5(x) ls 5 notitle

if (writepdf!=1) {
  pause -1
}

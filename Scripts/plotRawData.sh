echo "set terminal png; set output '~/_all_growth_curves.png'; set xlabel 'time / days'; set ylabel 'spheroid radius / micrometer'; plot '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_I_GC.dat'w e t'[G]=1mM, [O]=0.28mM', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_II_GC.dat'w e t'[G]=5mM, [O]=0.28mM', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_III_GC.dat'w e t'[G]=25mM, [O]=0.28mM', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_IV_GC.dat'w e t'[G]=25mM, [O]=0.07mM'" | gnuplot

echo "set terminal png; set output '~/_all_ECM_curves.png'; set xlabel 'time / days'; set ylabel 'spheroid radius / micrometer'; plot '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_II_T3_ECM.dat'w e t'[G]=5mM, [O]=0.28mM, T3', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_II_T4_ECM.dat'w e t'[G]=5mM, [O]=0.28mM, T4', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_III_T3_ECM.dat'w e t'[G]=25mM, [O]=0.28mM, T3', '/Users/jagiella/Documents/MATLAB/qsub_demo/SK-MES1_III_T4_ECM.dat'w e t'[G]=25mM, [O]=0.07mM, '" | gnuplot
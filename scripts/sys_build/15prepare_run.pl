	system("mkdir -p ./traj_1");
	open OUTPUT, ">./traj_1/run_gpu.sh" or die "can not create!\n";
	print OUTPUT "#!/bin/bash
np_n=`nproc`
np_n2=\$((\$np_n/2-4))
#np_n2=28
source ~/software/amber24/amber.sh
export CUDA_VISIBLE_DEVICES=0

export prmtop=system_start.prmtop
export name=system_start
mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 01_Min.in -o 01_Min.out -p \$prmtop -c \${name}.inpcrd -r 01_Min.rst -x 01_Min.nc -ref  \${name}.inpcrd
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 02_Heat.in -o 02_Heat.out -p \$prmtop -c 01_Min.rst -r 02_Heat.rst -x 02_Heat.nc -ref 01_Min.rst
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 03_Heat2.in -o 03_Heat2.out -p \$prmtop -c 02_Heat.rst -r 03_Heat2.rst -x 03_Heat2.nc -ref 02_Heat.rst
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 04_NPT.in -o 04_NPT.out -p \$prmtop -c 03_Heat2.rst -r 04_NPT.rst -x 04_NPT.nc -ref 03_Heat2.rst
sleep 3s

\$AMBERHOME/bin/pmemd.cuda -O -i 05_Hold.in -o 05_Hold_1.out -p \$prmtop -c 04_NPT.rst -r 05_Hold_1.rst -x 05_Hold_1.nc -ref 04_NPT.rst -inf 05_Hold_1.mdinfo
sleep 3s

for i in {2..10};
do	
	let h=\$i-1
	\$AMBERHOME/bin/pmemd.cuda -O -i 05_Hold.in -o 05_Hold_\$i.out -p \$prmtop -c 05_Hold_\$h.rst -r 05_Hold_\$i.rst -x 05_Hold_\$i.nc -ref 05_Hold_\$h.rst -inf 05_Hold_\$i.mdinfo
	sleep 3s
done

export name=1
\$AMBERHOME/bin/pmemd.cuda -O -i 06_Prod_\$name.in -o 06_Prod_\$name.out -p \$prmtop -c 05_Hold_10.rst -r 06_Prod_\$name.rst -x 06_Prod_\$name.nc -inf 06_Prod_\$name.mdinfo
	";
	close OUTPUT;



	open OUTPUT, ">./traj_1/01_Min.in" or die "can not create!\n";
	print OUTPUT "Minimization
 &cntrl
  imin=1,
  maxcyc=10000,
  ncyc=5000,
  ntb=1,
  ntp=0,
  ntf=1,
  ntc=1,
  ntpr=50,
  ntwr=2000,
  cut=12.0,
  fswitch=10.0,
  ntr=1,
  restraint_wt=5.0, restraintmask='\@C=',
 /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/02_Heat.in" or die "can not create!\n";
	print OUTPUT "Heating 100K
 &cntrl
  imin=0,
  ntx=1,
  irest=0,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=10000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=1000,
  ntwr=1000,
  ntwx=1000,
  dt=0.001,
  nmropt=1,
  ntb=1,
  ntp=0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=2.0, restraintmask='\@C=',
  tempi = 0, temp0 = 100,
  iwrap = 1,
 /
 &wt
  type='TEMP0',
  istep1=0,
  istep2=10000,
  value1=0.0,
  value2=100.0 /
 &wt type='END' /\n";
 	close OUTPUT;


	open OUTPUT, ">./traj_1/03_Heat2.in" or die "can not create!\n";
	print OUTPUT "Heating 310K
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=50000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=10000,
  ntwr=10000,
  ntwx=10000,
  dt=0.001,
  nmropt=1,
  ntb=1,
  ntp=0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=2.0, restraintmask='\@C='
  tempi =100, temp0 = 310,
  iwrap = 1,
 /
 &wt type='TEMP0', istep1=0, istep2=40000, value1=100.0, value2=310.0 /
 &wt type='TEMP0', istep1=40001, istep2=50000, value1=310.0, value2=310.0 /
 &wt type='END' /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/04_NPT.in" or die "can not create!\n";
	print OUTPUT "NPT 310K
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=50000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=5000,
  ntwr=5000,
  ntwx=5000,
  dt=0.002,
  temp0=310.0,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  csurften=3,
  gamma_ten=0.0,
  ninterface=2,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=1.0, restraintmask='\@C=',
  iwrap = 1,
 /
 /
 &ewald
  skinnb=3.0,
 /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/05_Hold.in" or die "can not create!\n";
	print OUTPUT "Production 310K 500ps
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=250000,
  ntt=3,
  gamma_ln=1.0,
  temp0=310.0,
  ntpr=5000,
  ntwr=5000,
  ntwx=5000,
  dt=0.002,
  ig=-1,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  csurften=3,
  gamma_ten=0.0,
  ninterface=2,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  iwrap = 1,
 /
 /
 &ewald
  skinnb=5.0,
 /\n";
 	close OUTPUT;


	open OUTPUT, ">./traj_1/06_Prod_1.in" or die "can not create!\n";
	print OUTPUT "Production 310K 500ns
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=250000000,
  ntt=3,
  gamma_ln=1.0,
  temp0=310.0,
  ntpr=50000,
  ntwr=50000,
  ntwx=50000,
  dt=0.002,
  ig=-1,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  csurften=3,
  gamma_ten=0.0,
  ninterface=2,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  iwrap = 1,
 /
 /
 &ewald
  skinnb=5.0,
 /\n";
	close OUTPUT;



	open OUTPUT, ">./traj_1/run_gpu_old.sh" or die "can not create!\n";
	print OUTPUT "#!/bin/bash
np_n=`nproc`
np_n2=\$((\$np_n/2-4))
#np_n2=28
source ~/software/amber24/amber.sh
export CUDA_VISIBLE_DEVICES=0

export prmtop=system_start.prmtop
export name=system_start
mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 01_Min_old.in -o 01_Min.out -p \$prmtop -c \${name}.inpcrd -r 01_Min.rst -x 01_Min.nc -ref  \${name}.inpcrd
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 02_Heat_old.in -o 02_Heat.out -p \$prmtop -c 01_Min.rst -r 02_Heat.rst -x 02_Heat.nc -ref 01_Min.rst
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 03_Heat2_old.in -o 03_Heat2.out -p \$prmtop -c 02_Heat.rst -r 03_Heat2.rst -x 03_Heat2.nc -ref 02_Heat.rst
sleep 3s

mpirun --use-hwthread-cpus -np \$np_n2 \$AMBERHOME/bin/pmemd.MPI -O -i 04_NPT_old.in -o 04_NPT.out -p \$prmtop -c 03_Heat2.rst -r 04_NPT.rst -x 04_NPT.nc -ref 03_Heat2.rst
sleep 3s

\$AMBERHOME/bin/pmemd.cuda -O -i 05_Hold_old.in -o 05_Hold_1.out -p \$prmtop -c 04_NPT.rst -r 05_Hold_1.rst -x 05_Hold_1.nc -ref 04_NPT.rst -inf 05_Hold_1.mdinfo
sleep 3s

for i in {2..10};
do	
	let h=\$i-1
	\$AMBERHOME/bin/pmemd.cuda -O -i 05_Hold_old.in -o 05_Hold_\$i.out -p \$prmtop -c 05_Hold_\$h.rst -r 05_Hold_\$i.rst -x 05_Hold_\$i.nc -ref 05_Hold_\$h.rst -inf 05_Hold_\$i.mdinfo
	sleep 3s
done

export name=1
\$AMBERHOME/bin/pmemd.cuda -O -i 06_Prod_old_\$name.in -o 06_Prod_\$name.out -p \$prmtop -c 05_Hold_10.rst -r 06_Prod_\$name.rst -x 06_Prod_\$name.nc -inf 06_Prod_\$name.mdinfo
	";
	close OUTPUT;



	open OUTPUT, ">./traj_1/01_Min_old.in" or die "can not create!\n";
	print OUTPUT "Minimization
 &cntrl
  imin=1,
  maxcyc=10000,
  ncyc=5000,
  ntb=1,
  ntp=0,
  ntf=1,
  ntc=1,
  ntpr=50,
  ntwr=2000,
  cut=12.0,
  fswitch=10.0,
  ntr=1,
  restraint_wt=5.0, restraintmask='\@C=',
 /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/02_Heat_old.in" or die "can not create!\n";
	print OUTPUT "Heating 100K
 &cntrl
  imin=0,
  ntx=1,
  irest=0,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=10000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=1000,
  ntwr=1000,
  ntwx=1000,
  dt=0.001,
  nmropt=1,
  ntb=1,
  ntp=0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=2.0, restraintmask='\@C=',
  tempi = 0, temp0 = 100,
  iwrap = 1,
 /
 &wt
  type='TEMP0',
  istep1=0,
  istep2=10000,
  value1=0.0,
  value2=100.0 /
 &wt type='END' /\n";
 	close OUTPUT;


	open OUTPUT, ">./traj_1/03_Heat2_old.in" or die "can not create!\n";
	print OUTPUT "Heating 310K
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=50000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=10000,
  ntwr=10000,
  ntwx=10000,
  dt=0.001,
  nmropt=1,
  ntb=1,
  ntp=0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=2.0, restraintmask='\@C='
  tempi =100, temp0 = 310,
  iwrap = 1,
 /
 &wt type='TEMP0', istep1=0, istep2=40000, value1=100.0, value2=310.0 /
 &wt type='TEMP0', istep1=40001, istep2=50000, value1=310.0, value2=310.0 /
 &wt type='END' /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/04_NPT_old.in" or die "can not create!\n";
	print OUTPUT "NPT 310K
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=50000,
  ntt=3,
  gamma_ln=1.0,
  ntr=1,
  ig=-1,
  ntpr=5000,
  ntwr=5000,
  ntwx=5000,
  dt=0.002,
  temp0=310.0,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  restraint_wt=1.0, restraintmask='\@C=',
  iwrap = 1,
 /
 /
 &ewald
  skinnb=3.0,
 /\n";
	close OUTPUT;


	open OUTPUT, ">./traj_1/05_Hold_old.in" or die "can not create!\n";
	print OUTPUT "Production 310K 500ps
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=250000,
  ntt=3,
  gamma_ln=1.0,
  temp0=310.0,
  ntpr=5000,
  ntwr=5000,
  ntwx=5000,
  dt=0.002,
  ig=-1,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  iwrap = 1,
 /
 /
 &ewald
  skinnb=5.0,
 /\n";
 	close OUTPUT;


	open OUTPUT, ">./traj_1/06_Prod_old_1.in" or die "can not create!\n";
	print OUTPUT "Production 310K 500ns
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=250000000,
  ntt=3,
  gamma_ln=1.0,
  temp0=310.0,
  ntpr=50000,
  ntwr=50000,
  ntwx=50000,
  dt=0.002,
  ig=-1,
  ntb=2,
  ntp=3,
  barostat=2,
  pres0=1.0,
  taup=2.0,
  cut=12.0,
  fswitch=10.0,
  ioutfm=1,
  ntxo=2,
  iwrap = 1,
 /
 /
 &ewald
  skinnb=5.0,
 /\n";
	close OUTPUT;




	system("cp system_start.prmtop ./traj_1/");
	system("cp system_start.inpcrd ./traj_1/");
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

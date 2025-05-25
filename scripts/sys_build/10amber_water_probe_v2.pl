
	if ($ARGV[-1] == 0){
		die "not add waters into the pocket!";
	}
	elsif ($ARGV[-1] == 1){
		print "generate the waters into pocket!";
	}
	else {
		system("cp $ARGV[-1] wats_inhole_del_2.pdb ");
		die "water pdb is provided!\n";
	}
	
	@tm_ids = ();
	for ($ag_i = 0; $ag_i < @ARGV-1; $ag_i++){
		$tm_ids[0] .= $ARGV[$ag_i];
	}
      open INPUT, "complex_prep_hs.pdb" or die "can not open!\n";
      open OUTPUT, ">chain_$tm_ids[0]_prep_hs.pdb" or die "can not create!\n";
      while(chomp($line=<INPUT>)){
	      @items = split /\s+/, $line;
	      @gezis = split //, $line;
	      if ($items[0] eq "ATOM"){
	            for ($ag_i = 0; $ag_i < @ARGV-1; $ag_i++){
			    if ($gezis[21] eq $ARGV[$ag_i]){
				print OUTPUT "$line\n";
			    }
	            		
	            }
		}
	      if ($items[0] eq "TER"){
			print OUTPUT "$line\n";
	      }
      }
      close INPUT;
      close OUTPUT;


	system("pdb4amber -i chain_$tm_ids[0]_prep_hs.pdb -o chain_$tm_ids[0]_for_rism.pdb -y -d --reduce");
	open OUTPUT, ">leap_for_rism.in" or die "can not create leap file!\n";
	print OUTPUT "
source leaprc.gaff
source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
rec = loadpdb chain_$tm_ids[0]_for_rism.pdb
solvatebox rec TIP3PBOX 11.0  
saveAmberParm rec chain_$tm_ids[0]_for_rism.prmtop chain_$tm_ids[0]_for_rism.inpcrd
quit\n";
	close OUTPUT;

	system("tleap -s -f leap_for_rism.in");
	system(" ambpdb -p chain_$tm_ids[0]_for_rism.prmtop -c chain_$tm_ids[0]_for_rism.inpcrd > chain_$tm_ids[0]_for_rism_box.pdb");


	open OUTPUT, ">mdin.rism" or die "can not create rism file!\n";
	print OUTPUT " single-point 3D-RISM calculation using the sander interface
&cntrl
   ntx=1, nstlim=0, irism=1,
/
&rism
   periodic='pme',
   closure='kh', tolerance=1e-6,
   grdspc=0.35,0.35,0.35, centering=0,
   mdiis_del=0.4, mdiis_nvec=20, maxstep=5000, mdiis_restart=50,
   solvcut=9.0,
   verbose=2, npropagate=0,
   apply_rism_force=0,
   volfmt='dx', ntwrism=1,
/\n";
	close OUTPUT;

	system("sander -O -i mdin.rism -o chain_$tm_ids[0]_for_rism.kh.r3d -p chain_$tm_ids[0]_for_rism.prmtop -c chain_$tm_ids[0]_for_rism.inpcrd -xvv cSPCE_kh.xvv -guv chain_$tm_ids[0]_for_rism.kh");
	system("metatwist --dx chain_$tm_ids[0]_for_rism.kh.O.0.dx --species O --convolve 4 --sigma 1.0 --odx chain_$tm_ids[0]_for_rism.kh.O.dx > chain_$tm_ids[0]_for_rism.lp");
	system("metatwist --dx chain_$tm_ids[0]_for_rism.kh.O.0.dx --ldx chain_$tm_ids[0]_for_rism.kh.O.dx --map blobsper --species O WAT --bulk 55.55 --threshold 0.5 > chain_$tm_ids[0]_for_rism.blobs");
	system("grep -v TER chain_$tm_ids[0]_for_rism.kh.O.0-chain_$tm_ids[0]_for_rism.kh.O-blobs-centroid.pdb > wats.pdb");


	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

	$mid_itv = 0.1;#mid flux range
	$rmsd_itv = 0.5; #rmsd flux range
	$plau_st = 0.90;#percent of rmsd in flux range
	$density_st = 0.5; #rmsd density threshold
	$frame_itv = 100; #frame num threshold to check rmsd flux
	$thick_z_default = 40;
	@mem_name_list = ("CHL","OL","PA", "PC");
	#input : complex start and end resi num, lig start and end resi num, start system prmtop and inpcrd, production trajectories.....
	
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

	@receptor_num_list = ();
	@ligand_num_list = ();
	@traj_list = ();
	@ragion_list = ();
	$itv_num = 0;
	$mem_mark = 0;
	for ($argv_i = 0; $argv_i < @ARGV ; $argv_i++){
		if ($ARGV[$argv_i] eq "-h"){
			print "run_mmpbsa.pl requires Amber22 R and MDAnalysis. please deactivate conda\n";
			print "	-r	start and end residue numbers of receptor, seperated by spaces.
	-l	start and end residue numbers of ligand, seperated by spaces.
	-y	trajectories for mmpbsa/mmgbsa, seperated by spaces.
	-f	start and stop frame num of trajectories if need to set. This set will overlap auto start frame.
	-p	start topology file (prmtop) of MD system.
	-c	start structure file (inpcrd) of MD system.
	-t	frame interval number; default is 1; high interval time save calculation time, but may reduce accuracy.
	-n	threads number for MMPBSA.py
	-m	whether it is a membrane protein? default is 0(no).
       	";
		die "	-h	print above information\n";
		}
		if ($ARGV[$argv_i] eq "-r"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					#print "$ARGV[$argv_i2]\n";
					last;
				}
				else{
					#print "$ARGV[$argv_i2]\n";
					push @receptor_num_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-l"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @ligand_num_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-y"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @traj_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-f"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @ragion_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-p"){
			$topo_file = $ARGV[$argv_i+1];
		}
		if ($ARGV[$argv_i] eq "-c"){
			$crd_file = $ARGV[$argv_i+1];
		}
		if ($ARGV[$argv_i] eq "-t"){
			$itv_num = $ARGV[$argv_i+1];
		}
		if ($ARGV[$argv_i] eq "-n"){
			$thread_num = $ARGV[$argv_i+1];
		}
		if ($ARGV[$argv_i] eq "-m"){
			$mem_mark = $ARGV[$argv_i+1];
		}
	}
	system("sleep 5s");
	system("ambpdb -p $topo_file -c $crd_file > system_start.pdb");
	system("sleep 5s");
		
	@atom_count_list = ();
	@resi_num_list = ();
	@atom_num_list = ();
	open INPUT, "system_start.pdb" or die "can not open!\n";
	$atom_count = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i<=5; $gezis_i++){
			$atom_mark.= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$atom_count++;
			push @atom_count_list, $atom_count;

			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i<=25; $gezis_i++){
				$resi_num.= $gezis[$gezis_i];
			}
			push @resi_num_list, $resi_num;

			$atom_num = undef;
			for ($gezis_i = 6; $gezis_i<=10; $gezis_i++){
				$atom_num.= $gezis[$gezis_i];
			}
			push @atom_num_list, $atom_num;
		}
	}
	close INPUT;


	$receptor_start_mark = 0;
	$receptor_end_mark = 0;
	$ligand_start_mark = 0;
	$ligand_end_mark = 0;

	@receptor_range = ();
	for ($num_i = 0; $num_i < @receptor_num_list; $num_i+=2){
		for ($count_i = 0; $count_i < @atom_count_list; $count_i++){
			if ($receptor_end_mark == 1){
				last;
			}
			if (($resi_num_list[$count_i] == $receptor_num_list[$num_i]) && ($receptor_start_mark == 0)) {
				if ($atom_count_list[$count_i] == $atom_num_list[$count_i]){
					push @receptor_range, $atom_count_list[$count_i];  #not to use atom_num_list, or add int()!!!!!!!
					$receptor_start_mark = 1;
				}
				else{
					die "atom_count wrong!\n";
				}
			}
			if (($resi_num_list[$count_i] == $receptor_num_list[$num_i+1] + 1 ) && ($receptor_start_mark == 1)) {
				if ($atom_count_list[$count_i] == $atom_num_list[$count_i]){
					push @receptor_range, $atom_count_list[$count_i] - 1;
					$receptor_start_mark = 0;
					$receptor_end_mark = 1;
				}
				else{
					die "atom_count wrong!\n";
				}
			}
		}
	}
	if ($receptor_end_mark == 0){
		push @receptor_range, $atom_count_list[-1];
	}
	@ligand_range = ();
	for ($num_i = 0; $num_i < @ligand_num_list; $num_i+=2){
		for ($count_i = 0; $count_i < @atom_count_list; $count_i++){
			if ($ligand_end_mark == 1){
				last;
			}
			if (($resi_num_list[$count_i] == $ligand_num_list[$num_i]) && ($ligand_start_mark == 0)) {
				if ($atom_count_list[$count_i] == $atom_num_list[$count_i]){
					push @ligand_range, $atom_count_list[$count_i];
					$ligand_start_mark = 1;
				}
				else{
					die "atom_count wrong!\n";
				}
			}
			if (($resi_num_list[$count_i] == $ligand_num_list[$num_i+1] + 1 ) && ($ligand_start_mark == 1)) {
				if ($atom_count_list[$count_i] == $atom_num_list[$count_i]){
					push @ligand_range, $atom_count_list[$count_i] - 1;
					$ligand_start_mark = 0;
					$ligand_end_mark = 1;
				}
				else{
					die "atom_count wrong!\n";
				}
			}
		}
	}
	if ($ligand_end_mark == 0){
		push @ligand_range, $atom_count_list[-1];
	}



	open OUTPUT, ">cpptraj_parm.in" or die "can not create!\n";

	print OUTPUT "parm $topo_file\n";
	print OUTPUT "parmstrip !\@$receptor_range[0]-$receptor_range[1]";
	for ($rec_i = 2; $rec_i < @receptor_range; $rec_i += 2){
		print OUTPUT ",$receptor_range[$rec_i]-$receptor_range[$rec_i+1]";
	}
	print OUTPUT "\n";
	print OUTPUT "parmwrite out receptor_dry.prmtop\nrun\nclear all\n";

	print OUTPUT "parm $topo_file\n";
	print OUTPUT "parmstrip !\@$ligand_range[0]-$ligand_range[1]";
	for ($lig_i = 2; $lig_i < @ligand_range; $lig_i += 2){
		print OUTPUT ",$ligand_range[$lig_i]-$ligand_range[$lig_i+1]";
	}
	print OUTPUT "\n";
	print OUTPUT "parmwrite out ligand_dry.prmtop\nrun\nclear all\n";

	print OUTPUT "parm $topo_file\n";
	print OUTPUT "parmstrip !\@$receptor_range[0]-$receptor_range[1]";
	for ($rec_i = 2; $rec_i < @receptor_range; $rec_i += 2){
		print OUTPUT ",$receptor_range[$rec_i]-$receptor_range[$rec_i+1]";
	}
	for ($lig_i = 0; $lig_i < @ligand_range; $lig_i += 2){
		print OUTPUT ",$ligand_range[$lig_i]-$ligand_range[$lig_i+1]";
	}
	print OUTPUT "\n";
	print OUTPUT "parmwrite out complex_dry.prmtop\nrun\nclear all\n";

	close OUTPUT;


	system("sleep 5s");
	system("cpptraj -i cpptraj_parm.in");

	@traj_energy_list=();
	
	open OUTMID, ">top_mid_list.txt" or die "can not create!\n";

	for ($input_i = 0; $input_i < @traj_list; $input_i++){
		$true_traj_num = $input_i + 1;
		open OUTPUT1, ">rmsd_hist_$true_traj_num.r" or die "can not create!\n";
		open OUTPUT, ">cpptraj_rmsd_$true_traj_num.in" or die "can not create!\n";
		print OUTPUT "parm $topo_file\n";
		print OUTPUT "trajin $traj_list[$input_i]\n";
		print OUTPUT "autoimage anchor !\@H=\nrun\n";
		print OUTPUT "strip !\@$receptor_range[0]-$receptor_range[1]";
		for ($rec_i = 2; $rec_i < @receptor_range; $rec_i += 2){
			print OUTPUT ",$receptor_range[$rec_i]-$receptor_range[$rec_i+1]";
		}
		for ($lig_i = 0; $lig_i < @ligand_range; $lig_i += 2){
			print OUTPUT ",$ligand_range[$lig_i]-$ligand_range[$lig_i+1]";
		}
		print OUTPUT "\n";

		print OUTPUT "trajout complex_dry_$true_traj_num.nc\nrun\nclear all\n";
		print OUTPUT "parm complex_dry.prmtop\n";
		print OUTPUT "trajin complex_dry_$true_traj_num.nc\n";
		#print OUTPUT "reference complex_start_dry.pdb name complex_start_dry\n";
		#print OUTPUT "rms !\@H= out rmsd_fit_start.dat ref complex_start_dry\n";
		print OUTPUT "rms !\@H= out rmsd_fit_$true_traj_num.dat\nrun\nclear all\n";
		print OUTPUT1 "rmsd_cr <- read.table\(\"rmsd_fit_$true_traj_num.dat\"\)\n";

	
		print OUTPUT1 "rmsd_ht<-hist\(rmsd_cr\[,2\]\)\n";

		print OUTPUT1 "top_mid<-rmsd_ht\$mids\[order\(rmsd_ht\$density,decreasing=TRUE\)\[1\]\]\n";
		print OUTPUT1 "write.table(rmsd_ht\$density, file=\"rmsd_density_$true_traj_num.dat\",col.names=FALSE)\n";
		print OUTPUT1 "write.table(rmsd_ht\$mids, file=\"rmsd_mids_$true_traj_num.dat\",col.names=FALSE)\n";
		print OUTPUT1 "write.table(top_mid, file=\"top_mid_$true_traj_num.dat\",col.names=FALSE)\n";
		close OUTPUT;
		close OUTPUT1;
        
		system("sleep 5s");
		system("cpptraj -i cpptraj_rmsd_$true_traj_num.in");
		system("sleep 5s");
		system("Rscript rmsd_hist_$true_traj_num.r");
        
		
		open INPUT1, "top_mid_$true_traj_num.dat" or die "can not create!\n";
		while(chomp($line=<INPUT1>)){
			@items = split /\s+/, $line;
			$top_mid = $items[1]; 
		}
		close INPUT1;
        
		@rmsd_list = ();
		open INPUT, "rmsd_fit_$true_traj_num.dat" or die "can not create!\n";
		$line_c = 0;
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 1){
				@items = split /\s+/, $line;
				push @rmsd_list, $items[2];
			}
		}
		close INPUT;

        
		$start_cut_mark = 0;
		for ($rmsd_i = 0; $rmsd_i < @rmsd_list - $frame_itv; $rmsd_i++){
			if ($start_cut_mark == 0){
				if (($rmsd_list[$rmsd_i] >= $top_mid - $mid_itv) && ($rmsd_list[$rmsd_i] <= $top_mid + $mid_itv) ){
					$plau_count = 0;
					for ($rmsd_i2 = $rmsd_i + 1; $rmsd_i2 <= $rmsd_i + $frame_itv; $rmsd_i2++  ){
						if ( ($rmsd_list[$rmsd_i2] >= $rmsd_list[$rmsd_i] - $rmsd_itv) &&  ($rmsd_list[$rmsd_i2] <= $rmsd_list[$rmsd_i] + $rmsd_itv)){
							$plau_count++;
						}
					}
					if (($plau_count/$frame_itv) >= $plau_st){
						$traj_cut_point = $rmsd_i+1;	
						$start_cut_mark = 1;
					}
					
				}
			}
		}
		if (@ragion_list){
			print OUTMID "$true_traj_num\t$top_mid\t$traj_cut_point\n";
			print "find traj $true_traj_num cut point!\n";
			open OUTPUT, ">cpptraj_cutforenergy_$true_traj_num.in" or die "can not create!\n";
			print OUTPUT "parm complex_dry.prmtop\n";
			print OUTPUT "trajin complex_dry_$true_traj_num.nc\n";

			print OUTPUT "trajout complex_cut_$true_traj_num.nc start $ragion_list[0] stop $ragion_list[1]\n";
			print OUTPUT "\nrun\nclear all\n";
			close OUTPUT;
			system("sleep 5s");
			system("cpptraj -i cpptraj_cutforenergy_$true_traj_num.in");
			push @traj_energy_list, $true_traj_num;
		}
        
        
		elsif ($start_cut_mark == 1){
			print OUTMID "$true_traj_num\t$top_mid\t$traj_cut_point\n";
			print "find traj $true_traj_num cut point!\n";
			open OUTPUT, ">cpptraj_cutforenergy_$true_traj_num.in" or die "can not create!\n";
			print OUTPUT "parm complex_dry.prmtop\n";
			print OUTPUT "trajin complex_dry_$true_traj_num.nc\n";

			print OUTPUT "trajout complex_cut_$true_traj_num.nc start $traj_cut_point\n";
			print OUTPUT "\nrun\nclear all\n";
			close OUTPUT;
			system("sleep 5s");
			system("cpptraj -i cpptraj_cutforenergy_$true_traj_num.in");
			push @traj_energy_list, $true_traj_num;
		}
        
        

	}
	
	close OUTMID; 

	if (!$itv_num){
		$itv_a = 1;
	}
	else{
		$itv_a = $itv_num;
	}

	if (!$mem_mark){
		
		open OUTPUT, ">mmgbsa.in" or die "can not create!\n";
		print OUTPUT "&general
   interval=$itv_a,
   verbose=2, keep_files=0,
/
&gb
   igb=5, saltcon=0.150,
/
&decomp
  idecomp=1,
  dec_verbose=1,
/\n";
		close OUTPUT;

		open OUTPUT, ">mmpbsa.in" or die "can not create!\n";
		print OUTPUT "&general
   interval=$itv_a,
   verbose=2, keep_files=0,
/
&pb
   istrng=0.15, fillratio=4.0,
/
&decomp
  idecomp=1,
  dec_verbose=1,
/\n";
		close OUTPUT;
        
		open OUTPUT, ">run.sh" or die "can not create!\n";
		if ($thread_num == 1){
			print OUTPUT "MMPBSA.py -O -i mmgbsa.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmgbsa.dat -do deco_mmgbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";

			print OUTPUT "MMPBSA.py -O -i mmpbsa.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmpbsa.dat -do deco_mmpbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";
		
		}
		else{
			print OUTPUT "mpirun --use-hwthread-cpus  -np $thread_num MMPBSA.py.MPI -O -i mmgbsa.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmgbsa.dat -do deco_mmgbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";

			print OUTPUT "mpirun --use-hwthread-cpus  -np $thread_num MMPBSA.py.MPI -O -i mmpbsa.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmpbsa.dat -do deco_mmpbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";
			
		}
		close OUTPUT;
		system("sleep 5s");
		system("bash run.sh");
		#system("rm rmsd_density_*.dat rmsd_mids_*.dat top_mid_*.dat rmsd_hist_*.r cpptraj_cutforenergy_*.in cpptraj_rmsd_*.in");
	}
	else{
		open INPUT, "system_start.pdb" or die "can not open!\n";
		$pos_z_up = 0;
		$pos_z_down = 0;
		$pos_z_up_2 = 0;
		$pos_z_down_2 = 0;
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i<=5; $gezis_i++){
				$atom_mark.= $gezis[$gezis_i];
			}
			if ($atom_mark eq "ATOM  ") {
				$resi_name_clear = undef;
				for ($gezis_i = 17; $gezis_i<=19; $gezis_i++){
					if ($gezis[$gezis_i] ne " "){
						$resi_name_clear.= $gezis[$gezis_i];
					}
				}
				$atom_num = undef;
				for ($gezis_i = 6; $gezis_i<=10; $gezis_i++){
					$atom_num.= $gezis[$gezis_i];
				}

				for ($mem_i = 0; $mem_i < @mem_name_list; $mem_i++){
					if ($resi_name_clear eq "$mem_name_list[$mem_i]"){
						$pos_z = undef;
						for ($gezis_i = 46; $gezis_i<=53; $gezis_i++){
							$pos_z.= $gezis[$gezis_i];
						}
						if ($pos_z_up < $pos_z){
							$pos_z_up = $pos_z;
						}
						if ($pos_z_down > $pos_z){
							$pos_z_down = $pos_z;
						}
						last;
					}
				}
				if (($atom_num >= $receptor_range[0] ) && ($atom_num <= $receptor_range[1])){
					$pos_z = undef;
					for ($gezis_i = 46; $gezis_i<=53; $gezis_i++){
						$pos_z.= $gezis[$gezis_i];
					}
					if ($pos_z_up_2 < $pos_z){
						$pos_z_up_2 = $pos_z;
					}
					if ($pos_z_down_2 > $pos_z){
						$pos_z_down_2 = $pos_z;
					}
					last;
				
				}

			}
		}
		close INPUT;

		$center_z_1 = ($pos_z_up + $pos_z_down)/2;
		$center_z_2 = ($pos_z_up_2 + $pos_z_down_2)/2;
		if ($center_z_1 == 0){
			$center_z = $center_z_2;
		}
		else{
			$center_z = $center_z_1;
		}

		$thick_z_calc = $pos_z_up - $pos_z_down;
		if (($thick_z_calc > $thick_z_default) || ($thick_z_calc == 0)) {
			$thick_z = $thick_z_default;
		}
		else{
			$thick_z = $thick_z_calc;
		}
		open OUTPUT, ">mmpbsa_mem.in" or die "can not create!\n";
		print OUTPUT "&general
   interval=$itv_a,
   verbose=2, keep_files=0, use_sander=1, debug_printlevel=2,
/
&pb
   radiopt=0, indi=20.0, istrng=0.150,
   fillratio=1.25, ipb=1, nfocus=1,
   bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0,
   npbverb=1, solvopt=2, inp=1,
   memopt=1, emem=7.0, mctrdz=$center_z, mthick=$thick_z, poretype=1,
   maxarcdot=15000
/
&decomp
  idecomp=1,
  dec_verbose=1,
/\n";
		
		open OUTPUT, ">cpptraj_start_frame.in" or die "can not create!\n";
		print OUTPUT "parm $topo_file\n";
		print OUTPUT "trajin $crd_file \n";
		print OUTPUT "strip !\@$receptor_range[0]-$receptor_range[1]";
		for ($rec_i = 2; $rec_i < @receptor_range; $rec_i += 2){
			print OUTPUT ",$receptor_range[$rec_i]-$receptor_range[$rec_i+1]";
		}
		for ($lig_i = 0; $lig_i < @ligand_range; $lig_i += 2){
			print OUTPUT ",$ligand_range[$lig_i]-$ligand_range[$lig_i+1]";
		}
		print OUTPUT "\n";
		print OUTPUT "trajout complex_start.pdb pdb\nrun\n";
		close OUTPUT;
		system("sleep 5s");
		system("cpptraj -i cpptraj_start_frame.in");

		for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
			open OUTPUT, ">2align_traj_$traj_energy_list[$traj_i].py" or die "can not create!\n";
			print OUTPUT "import MDAnalysis as mda
from MDAnalysis.analysis import align
u = mda.Universe('complex_dry.prmtop','complex_cut_$traj_energy_list[$traj_i].nc')
rf = mda.Universe('complex_dry.prmtop','complex_start.pdb')
align.AlignTraj(u, rf, select = 'name CA', filename='complex_cut_aligned_$traj_energy_list[$traj_i].nc',match_atom=True).run()\n";

			close OUTPUT;
			system("sleep 5s");
			system("echo wait trajectory alignment");
			system("python3 2align_traj_$traj_energy_list[$traj_i].py");
		}
		

		open OUTPUT, ">run.sh" or die "can not create!\n";
		if ($thread_num == 1){

			print OUTPUT "MMPBSA.py -O -i mmpbsa_mem.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmpbsa.dat -do deco_mmpbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_aligned_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";
		
		}
		else{

			print OUTPUT "mpirun --use-hwthread-cpus  -np $thread_num MMPBSA.py.MPI -O -i mmpbsa_mem.in -cp complex_dry.prmtop -rp receptor_dry.prmtop -o final_mmpbsa.dat -do deco_mmpbsa.dat   -lp ligand_dry.prmtop -y ";
			for ($traj_i = 0; $traj_i < @traj_energy_list; $traj_i++){
				print OUTPUT "complex_cut_aligned_$traj_energy_list[$traj_i].nc ";
			}
			print OUTPUT "\n";
			
		}
		close OUTPUT;
		system("sleep 5s");
		system("bash run.sh");
		#system("rm cpptraj_start_frame.in  2align_traj_*.py rmsd_density_*.dat rmsd_mids_*.dat top_mid_*.dat rmsd_hist_*.r cpptraj_cutforenergy_*.in cpptraj_rmsd_*.in");
	}

#########################for lig
	open OUTPUT, ">mmpbsa_extr_results_lig.txt" or die "can not create!\n";
	print OUTPUT "States\tTotal\t";
	for ($lig_i = 0; $lig_i < @ligand_num_list; $lig_i += 2){
		for ($res_i = $ligand_num_list[$lig_i]; $res_i <= $ligand_num_list[$lig_i+1]; $res_i++){
			print OUTPUT "$res_i\t";
		}
	}
	print OUTPUT "\n";

	open INPUT1, "./final_mmpbsa.dat\n" or die "can not open!\n";
	$mark_1 = 0;
	while(chomp($line=<INPUT1>)){
		if ($line=~ /^POISSON BOLTZMANN:.*/){
		#if ($line=~ /^GENERALIZED BORN:.*/){
			$mark_1 = 1;
			next;
		}
		if (($line=~ /^DELTA TOTAL.*/) && ($mark_1 == 1)) {
			@items = split /\s+/, $line;
			print OUTPUT "combine\t$items[2]\t";
		}

	}
	close INPUT1;

	open INPUT2, "./deco_mmpbsa.dat\n" or die "can not open!\n";
	$mark_2 = 0;
	$line_c = 0;
	@sem_list = ();
	@resi_num_name = ();
	while(chomp($line=<INPUT2>)){
		if ($line=~ /^Total Energy Decomposition:.*/){
			$mark_2 = 1;
			next;
		}
		if (($line=~ /^Sidechain Energy Decomposition:.*/) || ($line=~ /^Backbone Energy Decomposition:.*/)) {
			$mark_2 = 0;
			last;
		}
		if ($mark_2 == 1){
			$line_c++;
			if ($line_c >= 3){
				@gezis = split //, $line;
				$gbsa_resi_num = undef;
				for ($gezis_i = 3; $gezis_i<=6; $gezis_i++){
					$gbsa_resi_num.= $gezis[$gezis_i];
				}
				$gbsa_resi_name = undef;
				for ($gezis_i = 0; $gezis_i<=2; $gezis_i++){
					$gbsa_resi_name.= $gezis[$gezis_i];
				}

				for ($lig_i = 0; $lig_i < @ligand_num_list; $lig_i += 2){
					if (($gbsa_resi_num >= $ligand_num_list[$lig_i]) && ($gbsa_resi_num <= $ligand_num_list[$lig_i+1])){
						@items = split /\,/, $line;
						print OUTPUT "$items[-3]\t";
						push @resi_num_name, $gbsa_resi_name.int($gbsa_resi_num);
						@gezis_2 = split //, $items[-1];
						$resi_sem = undef;
						for ($gezis_i = 0; $gezis_i<@gezis_2-2; $gezis_i++){
							$resi_sem.= $gezis_2[$gezis_i];
						}
						push @sem_list , $resi_sem;
					}
				}
			
			}
		}
	}
	print OUTPUT "\n";
	close INPUT2;
	print OUTPUT "sem\tnon\t";
	for ($resi_i = 0; $resi_i < @sem_list; $resi_i++){
		print OUTPUT "$sem_list[$resi_i]\t";
	}
	print OUTPUT "\n";
	print OUTPUT "name\tnum\t";
	for ($resi_i = 0; $resi_i < @resi_num_name; $resi_i++){
		print OUTPUT "$resi_num_name[$resi_i]\t";
	}
	print OUTPUT "\n";
	close OUTPUT;

	open OUTPUT, ">3decom_bar_forai_lig_mmpbsa.r";
	print OUTPUT "args = commandArgs(trailingOnly=TRUE)
b_ener<-read.table(args[1],header=FALSE)
b_ener<-t(b_ener)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


colforbar<-c(\"red\",\"green\",\"violet\",\"blue\",\"orange\",\"limegreen\",\"deeppink4\",\"chartreuse4\",\"orangered\",\"darkblue\",\"deepskyblue2\",\"darkorange2\",\"yellowgreen\",\"springgreen3\",\"royalblue2\",\"slateblue1\",\"magenta2\",\"lightsalmon3\",\"cyan\",\"purple3\",\"plum1\",\"tan2\",\"palegreen4\",\"turquoise2\")

lig_per_energy<-as.numeric(b_ener[,2][3:length(b_ener[,1])])
lig_per_sem<-as.numeric(b_ener[,3][3:length(b_ener[,1])])
lig_per_name<-b_ener[,4][3:length(b_ener[,1])]

png(\"mmpbsa_result_top10_lig.png\",height=500,width=1000)
b_plot<-barplot(lig_per_energy[order(lig_per_energy)][1:10],col=\"orange\",names.arg=lig_per_name[order(lig_per_energy)][1:10],ylab=\"Binding energy per residue (kcal/mol)\",xlab=\"Residue\")
text(b_plot,lig_per_energy[order(lig_per_energy)][1:10],labels = round(lig_per_energy[order(lig_per_energy)][1:10],2),pos = 3)
error.bar(b_plot,lig_per_energy[order(lig_per_energy)][1:10],lig_per_sem[order(lig_per_energy)][1:10])
dev.off()\n";
	close OUTPUT;
	system("sleep 5s");
	system("Rscript 3decom_bar_forai_lig_mmpbsa.r mmpbsa_extr_results_lig.txt");


	if (!$mem_mark){
		open OUTPUT, ">mmgbsa_extr_results_lig.txt" or die "can not create!\n";
		print OUTPUT "States\tTotal\t";
		for ($lig_i = 0; $lig_i < @ligand_num_list; $lig_i += 2){
			for ($res_i = $ligand_num_list[$lig_i]; $res_i <= $ligand_num_list[$lig_i+1]; $res_i++){
				print OUTPUT "$res_i\t";
			}
		}
		print OUTPUT "\n";
        
		open INPUT1, "./final_mmgbsa.dat\n" or die "can not open!\n";
		$mark_1 = 0;
		while(chomp($line=<INPUT1>)){
			#if ($line=~ /^POISSON BOLTZMANN:.*/){
			if ($line=~ /^GENERALIZED BORN:.*/){
				$mark_1 = 1;
				next;
			}
			if (($line=~ /^DELTA TOTAL.*/) && ($mark_1 == 1)) {
				@items = split /\s+/, $line;
				print OUTPUT "combine\t$items[2]\t";
			}
        
		}
		close INPUT1;
        
		open INPUT2, "./deco_mmgbsa.dat\n" or die "can not open!\n";
		$mark_2 = 0;
		$line_c = 0;
		@sem_list = ();
		@resi_num_name = ();
		while(chomp($line=<INPUT2>)){
			if ($line=~ /^Total Energy Decomposition:.*/){
				$mark_2 = 1;
				next;
			}
			if (($line=~ /^Sidechain Energy Decomposition:.*/) || ($line=~ /^Backbone Energy Decomposition:.*/)) {
				$mark_2 = 0;
				last;
			}
			if ($mark_2 == 1){
				$line_c++;
				if ($line_c >= 3){
					@gezis = split //, $line;
					$gbsa_resi_num = undef;
					for ($gezis_i = 3; $gezis_i<=6; $gezis_i++){
						$gbsa_resi_num.= $gezis[$gezis_i];
					}
					$gbsa_resi_name = undef;
					for ($gezis_i = 0; $gezis_i<=2; $gezis_i++){
						$gbsa_resi_name.= $gezis[$gezis_i];
					}
        
					for ($lig_i = 0; $lig_i < @ligand_num_list; $lig_i += 2){
						if (($gbsa_resi_num >= $ligand_num_list[$lig_i]) && ($gbsa_resi_num <= $ligand_num_list[$lig_i+1])){
							@items = split /\,/, $line;
							print OUTPUT "$items[-3]\t";
							push @resi_num_name, $gbsa_resi_name.int($gbsa_resi_num);
							@gezis_2 = split //, $items[-1];
							$resi_sem = undef;
							for ($gezis_i = 0; $gezis_i<@gezis_2-2; $gezis_i++){
								$resi_sem.= $gezis_2[$gezis_i];
							}
							push @sem_list , $resi_sem;
						}
					}
				
				}
			}
		}
		print OUTPUT "\n";
		close INPUT2;
		print OUTPUT "sem\tnon\t";
		for ($resi_i = 0; $resi_i < @sem_list; $resi_i++){
			print OUTPUT "$sem_list[$resi_i]\t";
		}
		print OUTPUT "\n";
		print OUTPUT "name\tnum\t";
		for ($resi_i = 0; $resi_i < @resi_num_name; $resi_i++){
			print OUTPUT "$resi_num_name[$resi_i]\t";
		}
		print OUTPUT "\n";
		close OUTPUT;

		open OUTPUT, ">3decom_bar_forai_lig_mmgbsa.r";
		print OUTPUT "args = commandArgs(trailingOnly=TRUE)
b_ener<-read.table(args[1],header=FALSE)
b_ener<-t(b_ener)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


colforbar<-c(\"red\",\"green\",\"violet\",\"blue\",\"orange\",\"limegreen\",\"deeppink4\",\"chartreuse4\",\"orangered\",\"darkblue\",\"deepskyblue2\",\"darkorange2\",\"yellowgreen\",\"springgreen3\",\"royalblue2\",\"slateblue1\",\"magenta2\",\"lightsalmon3\",\"cyan\",\"purple3\",\"plum1\",\"tan2\",\"palegreen4\",\"turquoise2\")

lig_per_energy<-as.numeric(b_ener[,2][3:length(b_ener[,1])])
lig_per_sem<-as.numeric(b_ener[,3][3:length(b_ener[,1])])
lig_per_name<-b_ener[,4][3:length(b_ener[,1])]

png(\"mmgbsa_result_top10_lig.png\",height=500,width=1000)
b_plot<-barplot(lig_per_energy[order(lig_per_energy)][1:10],col=\"orange\",names.arg=lig_per_name[order(lig_per_energy)][1:10],ylab=\"Binding energy per residue (kcal/mol)\",xlab=\"Residue\")
text(b_plot,lig_per_energy[order(lig_per_energy)][1:10],labels = round(lig_per_energy[order(lig_per_energy)][1:10],2),pos = 3)
error.bar(b_plot,lig_per_energy[order(lig_per_energy)][1:10],lig_per_sem[order(lig_per_energy)][1:10])
dev.off()\n";
		close OUTPUT;
		system("sleep 5s");
		system("Rscript 3decom_bar_forai_lig_mmgbsa.r mmgbsa_extr_results_lig.txt");

	}

##############################for rec
	open OUTPUT, ">mmpbsa_extr_results_rec.txt" or die "can not create!\n";
	print OUTPUT "States\tTotal\t";
	for ($rec_i = 0; $rec_i < @receptor_num_list; $rec_i += 2){
		for ($res_i = $receptor_num_list[$rec_i]; $res_i <= $receptor_num_list[$rec_i+1]; $res_i++){
			print OUTPUT "$res_i\t";
		}
	}
	print OUTPUT "\n";

	open INPUT1, "./final_mmpbsa.dat\n" or die "can not open!\n";
	$mark_1 = 0;
	while(chomp($line=<INPUT1>)){
		if ($line=~ /^POISSON BOLTZMANN:.*/){
		#if ($line=~ /^GENERALIZED BORN:.*/){
			$mark_1 = 1;
			next;
		}
		if (($line=~ /^DELTA TOTAL.*/) && ($mark_1 == 1)) {
			@items = split /\s+/, $line;
			print OUTPUT "combine\t$items[2]\t";
		}

	}
	close INPUT1;

	open INPUT2, "./deco_mmpbsa.dat\n" or die "can not open!\n";
	$mark_2 = 0;
	$line_c = 0;
	@sem_list = ();
	@resi_num_name = ();
	while(chomp($line=<INPUT2>)){
		if ($line=~ /^Total Energy Decomposition:.*/){
			$mark_2 = 1;
			next;
		}
		if (($line=~ /^Sidechain Energy Decomposition:.*/) || ($line=~ /^Backbone Energy Decomposition:.*/)) {
			$mark_2 = 0;
			last;
		}
		if ($mark_2 == 1){
			$line_c++;
			if ($line_c >= 3){
				@gezis = split //, $line;
				$gbsa_resi_num = undef;
				for ($gezis_i = 3; $gezis_i<=6; $gezis_i++){
					$gbsa_resi_num.= $gezis[$gezis_i];
				}
				$gbsa_resi_name = undef;
				for ($gezis_i = 0; $gezis_i<=2; $gezis_i++){
					$gbsa_resi_name.= $gezis[$gezis_i];
				}

				for ($rec_i = 0; $rec_i < @receptor_num_list; $rec_i += 2){
					if (($gbsa_resi_num >= $receptor_num_list[$rec_i]) && ($gbsa_resi_num <= $receptor_num_list[$rec_i+1])){
						@items = split /\,/, $line;
						print OUTPUT "$items[-3]\t";
						push @resi_num_name, $gbsa_resi_name.int($gbsa_resi_num);
						@gezis_2 = split //, $items[-1];
						$resi_sem = undef;
						for ($gezis_i = 0; $gezis_i<@gezis_2-2; $gezis_i++){
							$resi_sem.= $gezis_2[$gezis_i];
						}
						push @sem_list , $resi_sem;
					}
				}
			
			}
		}
	}
	print OUTPUT "\n";
	close INPUT2;
	print OUTPUT "sem\tnon\t";
	for ($resi_i = 0; $resi_i < @sem_list; $resi_i++){
		print OUTPUT "$sem_list[$resi_i]\t";
	}
	print OUTPUT "\n";
	print OUTPUT "name\tnum\t";
	for ($resi_i = 0; $resi_i < @resi_num_name; $resi_i++){
		print OUTPUT "$resi_num_name[$resi_i]\t";
	}
	print OUTPUT "\n";
	close OUTPUT;

	open OUTPUT, ">3decom_bar_forai_rec_mmpbsa.r";
	print OUTPUT "args = commandArgs(trailingOnly=TRUE)
b_ener<-read.table(args[1],header=FALSE)
b_ener<-t(b_ener)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


colforbar<-c(\"red\",\"green\",\"violet\",\"blue\",\"orange\",\"limegreen\",\"deeppink4\",\"chartreuse4\",\"orangered\",\"darkblue\",\"deepskyblue2\",\"darkorange2\",\"yellowgreen\",\"springgreen3\",\"royalblue2\",\"slateblue1\",\"magenta2\",\"lightsalmon3\",\"cyan\",\"purple3\",\"plum1\",\"tan2\",\"palegreen4\",\"turquoise2\")

rec_per_energy<-as.numeric(b_ener[,2][3:length(b_ener[,1])])
rec_per_sem<-as.numeric(b_ener[,3][3:length(b_ener[,1])])
rec_per_name<-b_ener[,4][3:length(b_ener[,1])]

png(\"mmpbsa_result_top10_rec.png\",height=500,width=1000)
b_plot<-barplot(rec_per_energy[order(rec_per_energy)][1:10],col=\"orange\",names.arg=rec_per_name[order(rec_per_energy)][1:10],ylab=\"Binding energy per residue (kcal/mol)\",xlab=\"Residue\")
text(b_plot,rec_per_energy[order(rec_per_energy)][1:10],labels = round(rec_per_energy[order(rec_per_energy)][1:10],2),pos = 3)
error.bar(b_plot,rec_per_energy[order(rec_per_energy)][1:10],rec_per_sem[order(rec_per_energy)][1:10])
dev.off()\n";
	close OUTPUT;
	system("sleep 5s");
	system("Rscript 3decom_bar_forai_rec_mmpbsa.r mmpbsa_extr_results_rec.txt");

	if (!$mem_mark){
		open OUTPUT, ">mmgbsa_extr_results_rec.txt" or die "can not create!\n";
		print OUTPUT "States\tTotal\t";
		for ($rec_i = 0; $rec_i < @receptor_num_list; $rec_i += 2){
			for ($res_i = $receptor_num_list[$rec_i]; $res_i <= $receptor_num_list[$rec_i+1]; $res_i++){
				print OUTPUT "$res_i\t";
			}
		}
		print OUTPUT "\n";
        
		open INPUT1, "./final_mmgbsa.dat\n" or die "can not open!\n";
		$mark_1 = 0;
		while(chomp($line=<INPUT1>)){
			#if ($line=~ /^POISSON BOLTZMANN:.*/){
			if ($line=~ /^GENERALIZED BORN:.*/){
				$mark_1 = 1;
				next;
			}
			if (($line=~ /^DELTA TOTAL.*/) && ($mark_1 == 1)) {
				@items = split /\s+/, $line;
				print OUTPUT "combine\t$items[2]\t";
			}
        
		}
		close INPUT1;
        
		open INPUT2, "./deco_mmgbsa.dat\n" or die "can not open!\n";
		$mark_2 = 0;
		$line_c = 0;
		@sem_list = ();
		@resi_num_name = ();
		while(chomp($line=<INPUT2>)){
			if ($line=~ /^Total Energy Decomposition:.*/){
				$mark_2 = 1;
				next;
			}
			if (($line=~ /^Sidechain Energy Decomposition:.*/) || ($line=~ /^Backbone Energy Decomposition:.*/)) {
				$mark_2 = 0;
				last;
			}
			if ($mark_2 == 1){
				$line_c++;
				if ($line_c >= 3){
					@gezis = split //, $line;
					$gbsa_resi_num = undef;
					for ($gezis_i = 3; $gezis_i<=6; $gezis_i++){
						$gbsa_resi_num.= $gezis[$gezis_i];
					}
					$gbsa_resi_name = undef;
					for ($gezis_i = 0; $gezis_i<=2; $gezis_i++){
						$gbsa_resi_name.= $gezis[$gezis_i];
					}
        
					for ($rec_i = 0; $rec_i < @receptor_num_list; $rec_i += 2){
						if (($gbsa_resi_num >= $receptor_num_list[$rec_i]) && ($gbsa_resi_num <= $receptor_num_list[$rec_i+1])){
							@items = split /\,/, $line;
							print OUTPUT "$items[-3]\t";
							push @resi_num_name, $gbsa_resi_name.int($gbsa_resi_num);
							@gezis_2 = split //, $items[-1];
							$resi_sem = undef;
							for ($gezis_i = 0; $gezis_i<@gezis_2-2; $gezis_i++){
								$resi_sem.= $gezis_2[$gezis_i];
							}
							push @sem_list , $resi_sem;
						}
					}
				
				}
			}
		}
		print OUTPUT "\n";
		close INPUT2;
		print OUTPUT "sem\tnon\t";
		for ($resi_i = 0; $resi_i < @sem_list; $resi_i++){
			print OUTPUT "$sem_list[$resi_i]\t";
		}
		print OUTPUT "\n";
		print OUTPUT "name\tnum\t";
		for ($resi_i = 0; $resi_i < @resi_num_name; $resi_i++){
			print OUTPUT "$resi_num_name[$resi_i]\t";
		}
		print OUTPUT "\n";
		close OUTPUT;

		open OUTPUT, ">3decom_bar_forai_rec_mmgbsa.r";
		print OUTPUT "args = commandArgs(trailingOnly=TRUE)
b_ener<-read.table(args[1],header=FALSE)
b_ener<-t(b_ener)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


colforbar<-c(\"red\",\"green\",\"violet\",\"blue\",\"orange\",\"limegreen\",\"deeppink4\",\"chartreuse4\",\"orangered\",\"darkblue\",\"deepskyblue2\",\"darkorange2\",\"yellowgreen\",\"springgreen3\",\"royalblue2\",\"slateblue1\",\"magenta2\",\"lightsalmon3\",\"cyan\",\"purple3\",\"plum1\",\"tan2\",\"palegreen4\",\"turquoise2\")

rec_per_energy<-as.numeric(b_ener[,2][3:length(b_ener[,1])])
rec_per_sem<-as.numeric(b_ener[,3][3:length(b_ener[,1])])
rec_per_name<-b_ener[,4][3:length(b_ener[,1])]

png(\"mmgbsa_result_top10_rec.png\",height=500,width=1000)
b_plot<-barplot(rec_per_energy[order(rec_per_energy)][1:10],col=\"orange\",names.arg=rec_per_name[order(rec_per_energy)][1:10],ylab=\"Binding energy per residue (kcal/mol)\",xlab=\"Residue\")
text(b_plot,rec_per_energy[order(rec_per_energy)][1:10],labels = round(rec_per_energy[order(rec_per_energy)][1:10],2),pos = 3)
error.bar(b_plot,rec_per_energy[order(rec_per_energy)][1:10],rec_per_sem[order(rec_per_energy)][1:10])
dev.off()\n";
		close OUTPUT;
		system("sleep 5s");
		system("Rscript 3decom_bar_forai_rec_mmgbsa.r mmgbsa_extr_results_rec.txt");

	}


	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

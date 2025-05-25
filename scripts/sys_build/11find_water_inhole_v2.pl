	

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
	#$rec_pdb = $tm_ids[0];
	#$wat_pdb = $ARGV[1];
	$wat_pdb = "wats.pdb";
	$out_pdb = "wats_inhole.pdb";
	$hry_itv = 2;
	$atom_itv = 6;
	


	sub min_list {
		@swr_list = @_;
		$min_val = $swr_list[0];
		for ($swr_i = 1; $swr_i < @swr_list; $swr_i++){
			if ($swr_list[$swr_i] < $min_val){
				$min_val = $swr_list[$swr_i];
			}
		}
		return $min_val;
	}
	
	sub max_list {
		@swr_list = @_;
		$max_val = $swr_list[0];
		for ($swr_i = 1; $swr_i < @swr_list; $swr_i++){
			if ($swr_list[$swr_i] > $max_val){
				$max_val = $swr_list[$swr_i];
			}
		}
		return $max_val;
	}
	
	open INPUT, "chain_$tm_ids[0]_for_rism_box.pdb" or die "can not open!\n";
	open OUTPUT, ">chain_$tm_ids[0]_for_rism_dry.pdb" or die "can not create!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
			$atom_mark .= "$gezis[$gezis_i]";
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
			$resi_name .= "$gezis[$gezis_i]";
		}
		if ($atom_mark eq "ATOM  "){
			if ($resi_name eq "WAT"){
				last;
			}
			else{
				print OUTPUT "$line\n";
			}
		
		}
	}
	close INPUT;
	close OUTPUT;

	#open INPUT, "$rec_pdb" or die "can not open!\n";
	open INPUT, "chain_$tm_ids[0]_for_rism_dry.pdb" or die "can not open!\n";
	@pos_x_list = ();
	@pos_y_list = ();
	@pos_z_list = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
			$atom_mark .= "$gezis[$gezis_i]";
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
			$resi_name .= "$gezis[$gezis_i]";
		}
		if ($resi_name eq "WAT"){
			last;
		}
		if ($atom_mark eq "ATOM  "){
			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37; $gezis_i++){
				$pos_x .= "$gezis[$gezis_i]";
			}
			push @pos_x_list, $pos_x;

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45; $gezis_i++){
				$pos_y .= "$gezis[$gezis_i]";
			}
			push @pos_y_list, $pos_y;

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53; $gezis_i++){
				$pos_z .= "$gezis[$gezis_i]";
			}
			push @pos_z_list, $pos_z;
			
		}
	}
	close INPUT;

	$min_z = &min_list(@pos_z_list);
	$max_z = &max_list(@pos_z_list);

	print "min $min_z  max  $max_z\n";

	open OUTPUT, ">$out_pdb" or die "can create!\n";
	open INPUT, "$wat_pdb" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
			$atom_mark .= "$gezis[$gezis_i]";
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
			$resi_name .= "$gezis[$gezis_i]";
		}
		if (($atom_mark eq "ATOM  ")  && ($resi_name eq "WAT")){
			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37; $gezis_i++){
				$pos_x .= "$gezis[$gezis_i]";
			}

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45; $gezis_i++){
				$pos_y .= "$gezis[$gezis_i]";
			}

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53; $gezis_i++){
				$pos_z .= "$gezis[$gezis_i]";
			}
			
			$east_mark = 0;
			$west_mark = 0;
			$north_mark = 0;
			$south_mark = 0;
			for ($pos_i = 0; $pos_i < @pos_x_list; $pos_i++){
				if (($pos_z_list[$pos_i] <= $pos_z + $hry_itv) && ($pos_z_list[$pos_i] >= $pos_z - $hry_itv)) {
					if (($pos_x_list[$pos_i] <= $pos_x + $hry_itv) && ($pos_x_list[$pos_i] >= $pos_x - $hry_itv)) {
						if ($pos_y_list[$pos_i] > $pos_y + $atom_itv){
							$east_mark = 1;
							#print "east match $pos_x_list[$pos_i]\n";
						}
						#elsif(($pos_y_list[$pos_i] > $pos_y) && ($pos_y_list[$pos_i] < $pos_y + $atom_itv)){
						#	$east_mark = 0;
						#	last;
						#}
						elsif($pos_y_list[$pos_i] < $pos_y - $atom_itv){
							$west_mark = 1;
							
						}
						#elsif(($pos_y_list[$pos_i] < $pos_y) && ($pos_y_list[$pos_i] > $pos_y - $atom_itv)){
						#	$west_mark = 0;
						#	last;
						#}
					}
					if (($pos_y_list[$pos_i] <= $pos_y + $hry_itv) || ($pos_y_list[$pos_i] >= $pos_y - $hry_itv)) {
						if ($pos_x_list[$pos_i] > $pos_x + $atom_itv){
							$north_mark = 1;
						}
						#elsif(($pos_x_list[$pos_i] > $pos_x) && ($pos_x_list[$pos_i] < $pos_x + $atom_itv)){
						#	$north_mark = 0;
						#	last;
						#}
						elsif($pos_x_list[$pos_i] < $pos_x - $atom_itv){
							$south_mark = 1;
							
						}
						#elsif(($pos_x_list[$pos_i] < $pos_x) && ($pos_x_list[$pos_i] > $pos_x - $atom_itv)){
						#	$south_mark = 0;
						#	last;
						#}
					}
				}
			}
			
			print "($east_mark == 1) && ($west_mark == 1) && ($north_mark == 1) && ($south_mark == 1)\n";
			if (($east_mark == 1) && ($west_mark == 1) && ($north_mark == 1) && ($south_mark == 1)){
				print OUTPUT "$line\n";
			}
		}
	}
	close INPUT;

	close OUTPUT;

	open OUTPUT, ">zshow_wats_inhole.pml" or die "can not create!\n";
	print OUTPUT "load chain_$tm_ids[0]_for_rism_dry.pdb\n";
	print OUTPUT "load wats_inhole.pdb\n";
	print OUTPUT "save wats_inhole_del.pdb, wats_inhole\n";
	print OUTPUT "set sphere_scale, 0.6\n";
	print OUTPUT "show sphere, wats_inhole\n";
	print OUTPUT "bg white\nhide line\nshow cartoon\n";
	print OUTPUT "quit\n";
	close OUTPUT;

	system("pymol -qc zshow_wats_inhole.pml");
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

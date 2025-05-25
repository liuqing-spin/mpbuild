		

	if ($ARGV[0] == 0){
		die "not add waters into the pocket!";
	}
	elsif ($ARGV[0] == 1){
		print "generate the waters into pocket!";
	}
	else {
		system("cp $ARGV[0] wats_inhole_del_2.pdb ");
		die "water pdb is provided!\n";
	}
	
	sub calc_dist{
		@inpos_pair = @_;
		$dis_sq_sum = 0;
		for ($inpos_i = 0; $inpos_i < 3; $inpos_i++){
			$dist_swr = $inpos_pair[$inpos_i+3] - $inpos_pair[$inpos_i];
			$dist_sq_sum += $dist_swr**2;
		}
		$dist_out = sqrt($dist_sq_sum);
		return $dist_out;
	}

	open INPUT, "bilayer_complex_prep_hs_dry.pdb" or die "can not open!\n";
	@prot_posX_list = ();
	@prot_posY_list = ();
	@prot_posZ_list = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0;$gezis_i <= 5; $gezis_i++){
			$atom_mark.=$gezis[$gezis_i];
		}

		if ($atom_mark eq "ATOM  "){
			$posX = undef;
			for ($gezis_i = 30;$gezis_i <= 37; $gezis_i++){
				$posX.=$gezis[$gezis_i];
			}
			push @prot_posX_list, $posX;

			$posY = undef;
			for ($gezis_i = 38;$gezis_i <= 45; $gezis_i++){
				$posY.=$gezis[$gezis_i];
			}
			push @prot_posY_list, $posY;

			$posZ = undef;
			for ($gezis_i = 46;$gezis_i <= 53; $gezis_i++){
				$posZ.=$gezis[$gezis_i];
			}
			push @prot_posZ_list, $posZ;
		}
	}
	close INPUT;

	open INPUT, "wats_inhole_del_align.pdb" or die "can not open!\n";
	open OUTPUT, ">wats_inhole_del_2.pdb" or die "can not open!\n";
	@wat_posZ_list = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0;$gezis_i <= 5; $gezis_i++){
			$atom_mark.=$gezis[$gezis_i];
		}

		if ($atom_mark eq "ATOM  "){
			$posX = undef;
			for ($gezis_i = 30;$gezis_i <= 37; $gezis_i++){
				$posX.=$gezis[$gezis_i];
			}

			$posY = undef;
			for ($gezis_i = 38;$gezis_i <= 45; $gezis_i++){
				$posY.=$gezis[$gezis_i];
			}

			$posZ = undef;
			for ($gezis_i = 46;$gezis_i <= 53; $gezis_i++){
				$posZ.=$gezis[$gezis_i];
			}

			$forbid_mark = 0;
			for ($pos_i = 0; $pos_i < @prot_posX_list; $pos_i++){
				@list_for_calc = ($posX, $posY, $posZ, $prot_posX_list[$pos_i], $prot_posY_list[$pos_i], $prot_posZ_list[$pos_i]);
				$dist_check = &calc_dist(@list_for_calc);
				if ($dist_check < 3){
					$forbid_mark = 1;
					last;
				}
			}
			if ($forbid_mark == 0){
				print OUTPUT "$line\n";
			}
		}
	}
	close INPUT;
	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

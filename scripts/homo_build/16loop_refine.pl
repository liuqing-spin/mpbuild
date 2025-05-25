	@items = split /_/, $ARGV[1];

	while($ARGV[1]){
		system("sleep 10s");
		open INPUT, "$ARGV[1]_graft_conect.log" or die "can not open!\n" ;
		while(chomp($line=<INPUT>)){
			$log_line = $line;
		}
		close INPUT;
		if ($log_line=~/^DONE.*/){
			print "schrodinger prep is completed!\n";
			last;
		}

	}

	open INPUT, "$ARGV[1]_graft_conect_prep.pdb" or die "can not open!\n" ;
	open OUTPUT, ">$ARGV[1]_fe_model_align_graft_amber_relax.pdb" or die "can not open!\n" ;
	while(chomp($line=<INPUT>)){
		push @line_list, $line;
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$resi_name = undef;
			for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
				$resi_name .= $gezis[$gezis_i];
			}
			if ($resi_name eq "WAT"){
				last;
			}
			elsif ($resi_name eq "HIE"){
				$part_A = undef;
				for ($gezis_i = 0; $gezis_i <= 16 ; $gezis_i++){
					$part_A .= $gezis[$gezis_i];
				}
				$part_B = undef;
				for ($gezis_i = 22; $gezis_i < @gezis ; $gezis_i++){
					$part_B .= $gezis[$gezis_i];
				}
				printf OUTPUT "$part_A%5s$part_B\n", "HIS A";
			}
			else{
				$part_A = undef;
				for ($gezis_i = 0; $gezis_i <= 20 ; $gezis_i++){
					$part_A .= $gezis[$gezis_i];
				}
				$part_B = undef;
				for ($gezis_i = 22; $gezis_i < @gezis ; $gezis_i++){
					$part_B .= $gezis[$gezis_i];
				}
				printf OUTPUT "$part_A%1s$part_B\n", "A";
			}
		}
	}
	close OUTPUT;
	close INPUT;
	
	open OUTPUT, ">align_relax.pml" or die "can not create!\n" ;
	print OUTPUT "load $ARGV[1]_fe_model_align_graft_amber_relax.pdb\n";
	print OUTPUT "load $ARGV[1].pdb\n";
	print OUTPUT "align $ARGV[1]_fe_model_align_graft_amber_relax, $ARGV[1]\n";
	print OUTPUT "save $ARGV[1]_fe_model_relax_align.pdb,   $ARGV[1]_fe_model_align_graft_amber_relax\n";
	print OUTPUT "quit\n";
	close OUTPUT;

	system("pymol -qc align_relax.pml");



	system("cp $ARGV[1]_fe_model_relax_align.pdb  $ARGV[1]_fe_model_align_fix_refine.pdb");
	open INPUT, "$ARGV[1]_fe_model_relax_align.pdb" or die "can not open!\n" ;
	@line_list = ();
	@pos_x_list = ();
	@pos_y_list = ();
	@pos_z_list = ();
	@resi_num_list = ();
	while(chomp($line=<INPUT>)){
		push @line_list, $line;
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$atom_name = undef;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				$atom_name .= $gezis[$gezis_i];
			}

			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
				$pos_x .= $gezis[$gezis_i];
			}

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
				$pos_y .= $gezis[$gezis_i];
			}

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
				$pos_z .= $gezis[$gezis_i];
			}

			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			push @pos_x_list, $pos_x;
			push @pos_y_list, $pos_y;
			push @pos_z_list, $pos_z;
			push @resi_num_list, $resi_num;

		}
	}
	close INPUT;


	open INPUT, "$ARGV[0]" or die "can not open!\n" ;
	@line_list_2 = ();
	@pos_x_list_2 = ();
	@pos_y_list_2 = ();
	@pos_z_list_2 = ();
	@resi_num_list_2 = ();
	while(chomp($line=<INPUT>)){
		push @line_list_2, $line;
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$chain_match = 0;
		for ($cid_i = 2; $cid_i < @ARGV; $cid_i++){
			if ($gezis[21] eq $ARGV[$cid_i]){
				$chain_match = 1;
				last;
			}
		}
		if (($atom_mark eq "ATOM  ") && ($gezis[21] ne $items[1]) && ($chain_match == 1))  {
			$atom_name = undef;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				$atom_name .= $gezis[$gezis_i];
			}

			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
				$pos_x .= $gezis[$gezis_i];
			}

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
				$pos_y .= $gezis[$gezis_i];
			}

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
				$pos_z .= $gezis[$gezis_i];
			}

			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			push @pos_x_list_2, $pos_x;
			push @pos_y_list_2, $pos_y;
			push @pos_z_list_2, $pos_z;
			push @resi_num_list_2, $resi_num;

		}
	}
	close INPUT;


	$min_dist = 999999999999;
	$cut_off_dist = 0.5;
	@resi_num_ster = ();
	for ($pos_i = 0; $pos_i < @pos_x_list; $pos_i++){
		for ($pos_i_2 = 0; $pos_i_2 < @pos_x_list_2; $pos_i_2++){
			$atom_dist = sqrt(($pos_x_list[$pos_i] - $pos_x_list_2[$pos_i_2])**2 + ($pos_y_list[$pos_i] - $pos_y_list_2[$pos_i_2])**2 + ($pos_z_list[$pos_i] - $pos_z_list_2[$pos_i_2])**2 );
			if ($atom_dist < $min_dist){
				$min_dist = $atom_dist;
				$min_resi = $resi_num_list[$pos_i];
			}
			if ($atom_dist < $cut_off_dist){
				push @resi_num_ster, $resi_num_list[$pos_i];
			}
		}
		
	}
	print "min_dist; $min_dist; resi $min_resi\n";

	$em_seq = undef;
	$md_seq = undef;
	open INPUT, "chain_$items[1]_mm.aln" or die "can not open 2!\n";
	while(chomp($line=<INPUT>)){
		@re_item = split /\s+/, $line;
		if ($re_item[0] eq "chain_$items[1]_em"){
			$em_seq.=$re_item[1];
		}
		if ($re_item[0] eq "chain_$items[1]_md"){
			$md_seq.=$re_item[1];
		}
	}
	
	@md_seq_list = split //, $md_seq;
	@em_seq_list = split //, $em_seq;
	@resi_num_ster_fix = ();
	for ($em_i = 0; $em_i < @em_seq_list; $em_i++){
		if ($em_seq_list[$em_i] eq "-") {
			for ($ster_i = 0; $ster_i < @resi_num_ster; $ster_i++){
				if (($em_i + 1) == $resi_num_ster[$ster_i]){
					push @resi_num_ster_fix, $resi_num_ster[$ster_i];
					last;
				}
			}
		}
	}
	close INPUT;

	@refine_itv = ();
	push @refine_itv, $resi_num_ster_fix[0];
	#$itv_a = $resi_num_ster_fix[0];
	for ($resi_i = 1; $resi_i < @resi_num_ster_fix; $resi_i++){
		if ($resi_num_ster_fix[$resi_i] > $resi_num_ster_fix[$resi_i-1] + 2){
			push @refine_itv, $resi_num_ster_fix[$resi_i-1];
			push @refine_itv, $resi_num_ster_fix[$resi_i];
			#$itv_a = $resi_num_ster_fix[$resi_i];
		}
	}
	push @refine_itv, $resi_num_ster_fix[-1];

	print "resi_ster_list: @resi_num_ster_fix\n";
	print "refine_itv: @refine_itv\n";


	open OUTPUT, ">16output.txt" or die "can not open 2!\n";
	print OUTPUT "resi_ster: @resi_num_ster\n";
	print OUTPUT "resi_ster_list: @resi_num_ster_fix\n";
	print OUTPUT "refine_itv: @refine_itv\n";
	print OUTPUT "min_dist; $min_dist; resi $min_resi\n";
	close OUTPUT;

	$model_file = $ARGV[1]."_fe";
	for ($resi_i = 0; $resi_i < @refine_itv; $resi_i += 2){
		if ($refine_itv[$resi_i+1]){
			$loop_start = $refine_itv[$resi_i] - 2;
			$loop_end = $refine_itv[$resi_i+1] + 2;
			system("python3 16loop_refine.py $ARGV[1] $model_file $loop_start $loop_end ssbond_filter_model_$items[1].txt ");
			system("sleep 5s");
			$loop_check_file = $model_file."_loop_check.pml";
			system("pymol -qc $loop_check_file");

			system("sleep 5s");
			for ($file_i = 1; $file_i <= 20; $file_i++){
				open INPUT, "loop_refine_align_$file_i.pdb" or die "can not open!\n" ;
				@line_list_3 = ();
				@pos_x_list_3 = ();
				@pos_y_list_3 = ();
				@pos_z_list_3 = ();
				@resi_num_list_3 = ();
				while(chomp($line=<INPUT>)){
					push @line_list_3, $line;
					@gezis = split //, $line;
					$atom_mark = undef;
					for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
						$atom_mark .= $gezis[$gezis_i];
					}
					$resi_num = undef;
					for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
						$resi_num .= $gezis[$gezis_i];
					}
					if (($atom_mark eq "ATOM  ") && ($resi_num >= $loop_start) && ($resi_num <= $loop_end)){
						$atom_name = undef;
						for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
							$atom_name .= $gezis[$gezis_i];
						}
                        
						$pos_x = undef;
						for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
							$pos_x .= $gezis[$gezis_i];
						}
                        
						$pos_y = undef;
						for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
							$pos_y .= $gezis[$gezis_i];
						}
                        
						$pos_z = undef;
						for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
							$pos_z .= $gezis[$gezis_i];
						}
                        
						push @pos_x_list_3, $pos_x;
						push @pos_y_list_3, $pos_y;
						push @pos_z_list_3, $pos_z;
						push @resi_num_list_3, $resi_num;
                        
					}
				}
				close INPUT;
				$cut_off_dist = 2;
				@resi_num_ster = ();
				$dist_fail = 0;
				for ($pos_i_2 = 0; $pos_i_2 < @pos_x_list_2; $pos_i_2++){
					for ($pos_i_3 = 0; $pos_i_3 < @pos_x_list_3; $pos_i_3++){
						$atom_dist = sqrt(($pos_x_list_2[$pos_i_2] - $pos_x_list_3[$pos_i_3])**2 + ($pos_y_list_2[$pos_i_2] - $pos_y_list_3[$pos_i_3])**2 + ($pos_z_list_2[$pos_i_2] - $pos_z_list_3[$pos_i_3])**2 );
						if ($atom_dist < $cut_off_dist){
							$dist_fail = 1;
							print "loop_refine_align_$file_i.pdb loop refine fail!\n";
							last;
						}
					}
					if ($dist_fail == 1){
						last;
					}
					
				}
				if ($dist_fail == 0){
					print "loop_refine_align_$file_i.pdb loop refine ok!\n";
					system("cp loop_refine_align_$file_i.pdb loop_refine_result_$resi_i.pdb");
					system("cp loop_refine_result_$resi_i.pdb  $ARGV[1]_fe_model_align_fix_refine.pdb");
					system("mv $ARGV[1]_fe_model_relax_align.pdb  $ARGV[1]_fe_model_relax_align.pdb_$resi_i");
					system("cp loop_refine_result_$resi_i.pdb  $ARGV[1]_fe_model_relax_align.pdb");
					system("rm loop_refine_align*.pdb");
					system("rm $model_file.BL*.pdb");
					system("rm $model_file.DL*");
					#$model_file = "loop_refine_result_$resi_i" ;
					last;
				}


				
			}
		}
	}

	@items = split /_/, $ARGV[0];
	$C_N_dist_max = 1.6;  #this is 1.3 measured by pymol. 
	$C_N_dist_min = 1.4;  #this is 1.3 measured by pymol. 
	
	system("mv $ARGV[0].pdb $ARGV[0]_00.pdb");

 %aa_name = (
     "R"  =>   "ARG",
     "H"  =>   "HIS",
     "K"  =>   "LYS",
     "D"  =>   "ASP",
     "E"  =>   "GLU",
     "S"  =>   "SER",
     "T"  =>   "THR",
     "N"  =>   "ASN",
     "Q"  =>   "GLN",
     "C"  =>   "CYS",
     "G"  =>   "GLY",
     "P"  =>   "PRO",
     "A"  =>   "ALA",
     "V"  =>   "VAL",
     "I"  =>   "ILE",
     "L"  =>   "LEU",
     "M"  =>   "MET",
     "F"  =>   "PHE",
     "Y"  =>   "TYR",
     "W"  =>   "TRP",
 );
	$raw_seq = undef;
	$em_seq = undef;
	$raw_mark = 0;
	$em_mark = 0;
	@conv_list = ();
	@em_seq_list = ();
	@raw_seq_list = ();
	open INPUT, "chain_$items[1]_re_02.aln" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		@re_item = split /\s+/, $line;
		if ($re_item[0] eq "chain_$items[1]_raw"){
			$raw_seq=$re_item[1];
			$raw_mark = 1;
			@raw_chr = split //, $raw_seq;
			push @raw_seq_list, @raw_chr;
			$raw_chr_count = @raw_chr;
			next;
		}
		if ($re_item[0] eq "chain_$items[1]_em"){
			$em_seq=$re_item[1];
			$em_mark = 1;
			@em_chr = split //, $em_seq;
			push @em_seq_list, @em_chr;
			$em_chr_count = @em_chr;
			if ($em_chr_count == $raw_chr_count){
				$wp_mark = 0;
				for ($gezis_i = 0; $gezis_i < @gezis; $gezis_i++){
					if (($gezis[$gezis_i] eq " ") && ($wp_mark == 0)){
						$wp_mark = 1;
						next;
					}
					if (($gezis[$gezis_i] ne " ") && ($wp_mark == 1)){
						$cth_start_num = $gezis_i;
						last;
					}
				}
			}
			else{
				die "em num not eq raw num!\n";
			}	
			next;
		}
		if (($raw_mark == 1) && ($em_mark == 1)){
			#for ($gezis_i = 17; $gezis_i < $raw_chr_count + 17; $gezis_i++ ){
			for ($gezis_i = $cth_start_num; $gezis_i < $cth_start_num + $raw_chr_count; $gezis_i++ ){
				push @conv_list, $gezis[$gezis_i];
			}
		}
		$raw_mark = 0;
		$em_mark = 0;
	}
	close INPUT;

	$raw_chr_count_all = @raw_seq_list;
	$em_chr_count_all = @em_seq_list;
	$conv_count_all = @conv_list;

	print "raw count $raw_chr_count_all; em count $em_chr_count_all; conv count $conv_count_all\n";

	open INPUT, "$ARGV[0]_00.pdb" or die "can not open 4!\n" ;
	@C_pos_x=();
	@C_pos_y=();
	@C_pos_z=();
	@N_pos_x=();
	@N_pos_y=();
	@N_pos_z=();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
			$resi_num .= $gezis[$gezis_i];
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
			$resi_name .= $gezis[$gezis_i];
		}
		$atom_name_clear = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			if ($gezis[$gezis_i] ne " "){
				$atom_name_clear .= $gezis[$gezis_i];
			}
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
		if ($atom_name_clear eq "C"){
			push @C_pos_x, $pos_x;
			push @C_pos_y, $pos_y;
			push @C_pos_z, $pos_z;
		}
		if ($atom_name_clear eq "N"){
			push @N_pos_x, $pos_x;
			push @N_pos_y, $pos_y;
			push @N_pos_z, $pos_z;
		}

	}
	close INPUT;




	$em_resi_count = 0;
	@mut_resi_list = ();
	@mut_raw_name_list = ();
	@mut_em_name_list = ();
	@raw_corres_resi_num = ();
	@em_corres_resi_num = ();
	@em_resi_dele = ();
	@em_resi_dele_itr = ();
	$raw_resi_count=0;

	@resi_num_C_cnc = ();
	@resi_num_N_cnc = ();
	$itv_dele_index=-999;

	$dele_itr_start = 0;
	$dele_itr_end = 1;
	@resi_dele_itv_list = ();
	if (($raw_chr_count_all == $conv_count_all) && ($em_chr_count_all == $conv_count_all)){
		for ($raw_i = 0; $raw_i < @raw_seq_list; $raw_i++){
			if ($raw_seq_list[$raw_i] ne "-"){
				$raw_resi_count++;
			}
			if ($em_seq_list[$raw_i] ne "-"){
				$em_resi_count++;
				if ($conv_list[$raw_i] ne "*"){
					push @mut_resi_list, $em_resi_count;
					push @mut_raw_name_list, $raw_seq_list[$raw_i]; 
					push @mut_em_name_list, $em_seq_list[$raw_i]; 
				}

				if ($raw_seq_list[$raw_i] eq "-"){
					push @em_resi_dele, $em_resi_count;
					push @em_resi_dele_itr, $em_resi_count;
					if (($raw_seq_list[$raw_i-1] ne "-") && ($dele_itr_start==0) && ($dele_itr_end==1)){
						push @resi_dele_itv_list, $raw_resi_count;
						$dele_itr_start=1;
						$dele_itr_end=0;
					}
					if (($raw_seq_list[$raw_i+1] ne "-") && ($dele_itr_start==1) && ($dele_itr_end==0)){
						push @resi_dele_itv_list, $raw_resi_count+1;
						$dele_itr_start=0;
						$dele_itr_end=1;
					}
					$em_seq_list[$raw_i] = "-";
					$conv_list[$raw_i] = " ";
				}
				elsif($raw_i == $itv_dele_index){
					push @em_resi_dele, $em_resi_count;
					$em_seq_list[$raw_i] = "-";
					$conv_list[$raw_i] = " ";
					push @raw_corres_resi_num, $raw_resi_count;
					push @em_corres_resi_num, $em_resi_count;
				}
				else{
					push @raw_corres_resi_num, $raw_resi_count;
					push @em_corres_resi_num, $em_resi_count;
				}

				if (($em_seq_list[$raw_i+1] ne "-") && ($raw_i+1 < @raw_seq_list))   {
					$cndist = sqrt(($C_pos_x[$em_resi_count-1] - $N_pos_x[$em_resi_count])**2 + ($C_pos_y[$em_resi_count-1] - $N_pos_y[$em_resi_count])**2 + ($C_pos_z[$em_resi_count-1] - $N_pos_z[$em_resi_count])**2) ;
					if ($cndist > $C_N_dist_max ){
						if ($em_seq_list[$raw_i-1] eq "-"){
							push @em_resi_dele, $em_resi_count;
							$em_seq_list[$raw_i] = "-";
							$conv_list[$raw_i] = " ";
						}
						else{
							push @resi_num_C_cnc, $raw_resi_count;
							push @resi_num_N_cnc, $raw_resi_count+1;
						}
					}
				}
				else{
					for($raw_i2=$raw_i+2; $raw_i2<@raw_seq_list; $raw_i2++){
						if ($em_seq_list[$raw_i2] ne "-"){
							$cndist = sqrt(($C_pos_x[$em_resi_count-1] - $N_pos_x[$em_resi_count])**2 + ($C_pos_y[$em_resi_count-1] - $N_pos_y[$em_resi_count])**2 + ($C_pos_z[$em_resi_count-1] - $N_pos_z[$em_resi_count])**2) ;
							if ($cndist < $C_N_dist_min ){
								push @em_resi_dele, $em_resi_count;
								$em_seq_list[$raw_i] = "-";
								$conv_list[$raw_i] = " ";
								$itv_dele_index=$raw_i2;
								#push @em_resi_dele, $em_resi_count+1;
								#$em_seq_list[$raw_i2] = "-";
								#$conv_list[$raw_i2] = " ";
							}
							last;
						}
					}
					
				}

			}
		}



		#for ($raw_i = 0; $raw_i < @raw_seq_list; $raw_i++){
		#	if ($raw_seq_list[$raw_i] ne "-"){
		#		$raw_resi_count++;
		#	}
		#	if ($em_seq_list[$raw_i] ne "-"){
		#		$em_resi_count++;
		#		if ($raw_seq_list[$raw_i] eq "-"){
		#			push @em_resi_dele, $em_resi_count;
		#		}
		#		else{
		#			push @raw_corres_resi_num, $raw_resi_count;
		#			push @em_corres_resi_num, $em_resi_count;
                #
		#		}
		#		if ($conv_list[$raw_i] ne "*"){
		#			push @mut_resi_list, $em_resi_count;
		#			push @mut_raw_name_list, $raw_seq_list[$raw_i]; 
		#			push @mut_em_name_list, $em_seq_list[$raw_i]; 
		#		}
		#	}
		#}
	}
	else{
		die "total count not eq!\n";
	}
	
	print "@conv_list\n";
	#print "@mut_em_name_list\n";
	
	@em_pdb_corres_resi_num = ();

	open INPUT, "$ARGV[0]_00.pdb" or die "can not open 4!\n" ;
	open OUTPUT, ">$ARGV[0].pdb" or die "can not open 4!\n" ;
	$crt_resi_num = -9999999999999;
	$pdb_resi_count = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
			$resi_num .= $gezis[$gezis_i];
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
			$resi_name .= $gezis[$gezis_i];
		}
		$atom_name = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			$atom_name .= $gezis[$gezis_i];
		}
		$new_resi_mark = 0;
		if ($resi_num != $crt_resi_num){
			$crt_resi_num = $resi_num;
			$pdb_resi_count++;
			$new_resi_mark = 1;
		}
		$dele_match_itr = 0;
		for ($resi_i = 0; $resi_i < @em_resi_dele_itr; $resi_i++){
			if ($pdb_resi_count == $em_resi_dele_itr[$resi_i]) {
				$dele_match_itr = 1;
				last;
			}
		}
		if ($dele_match_itr == 1){
			next;
		}
		if ($new_resi_mark == 1){
			push @em_pdb_corres_resi_num, $crt_resi_num;
		}

		$dele_match = 0;
		for ($resi_i = 0; $resi_i < @em_resi_dele; $resi_i++){
			if ($pdb_resi_count == $em_resi_dele[$resi_i]) {
				$dele_match = 1;
				last;
			}
		}
		if ($dele_match == 1){
			next;
		}
		$mut_match = 0;
		$raw_resi_num = -9999;
		for ($resi_i = 0; $resi_i < @mut_resi_list; $resi_i++){
			if (($pdb_resi_count == $mut_resi_list[$resi_i]) && ($resi_name eq $aa_name{$mut_em_name_list[$resi_i]})) {
				$mut_match = 1;
				$corres_match = 0;
				for ($renu_i = 0; $renu_i < @em_corres_resi_num; $renu_i++){
					if ($em_corres_resi_num[$renu_i] == $pdb_resi_count){
						$raw_resi_num = $raw_corres_resi_num[$renu_i];
						$corres_match = 1;
						last;
					}
				}
				if ($corres_match == 0){
					die "residue number corresponding wrong!\n";
				}
				if (($atom_name eq " C  ") || ($atom_name eq " O  ") || ($atom_name eq " CA ") || ($atom_name eq " N  ")){
					$part_A = undef;
					for ($gezis_i = 0; $gezis_i <= 16 ; $gezis_i++){
						$part_A .= $gezis[$gezis_i];
					}
					$part_B = undef;
					for ($gezis_i = 20; $gezis_i <= 21 ; $gezis_i++){
						$part_B .= $gezis[$gezis_i];
					}
					$part_C = undef;
					for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
						$part_C .= $gezis[$gezis_i];
					}
					printf OUTPUT "$part_A%3s$part_B%4g$part_C\n", $aa_name{$mut_raw_name_list[$resi_i]}, $raw_resi_num;
				}
				last;
			}
		}
		if ($mut_match == 0){
			$corres_match = 0;
			for ($renu_i = 0; $renu_i < @em_corres_resi_num; $renu_i++){
				if ($em_corres_resi_num[$renu_i] == $pdb_resi_count){
					$raw_resi_num = $raw_corres_resi_num[$renu_i];
					$corres_match = 1;
					last;
				}
			}
			if ($corres_match == 0){
				die "residue number corresponding wrong!\n";
			}
			$part_A = undef;
			for ($gezis_i = 0; $gezis_i <= 21 ; $gezis_i++){
				$part_A .= $gezis[$gezis_i];
			}
			$part_C = undef;
			for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
				$part_C .= $gezis[$gezis_i];
			}
			printf OUTPUT "$part_A%4g$part_C\n",  $raw_resi_num;
		}
	}
	close INPUT;
	close OUTPUT;

	$pdb_resi_num_count = @em_pdb_corres_resi_num;
	$aln_resi_num_count = @em_corres_resi_num;
	if ($pdb_resi_num_count != $aln_resi_num_count){
		print  "pdb resi num is $pdb_resi_num_count, not equal aln resi num $aln_resi_num_count!\n"
	}
	open OUTPUT, ">$ARGV[0]_raw_vs_em_resi_num.txt" or die "can not open 4!\n" ;
	print OUTPUT "em_count\tem_pdb\traw\n";
	for ($raw_i = 0; $raw_i < @raw_corres_resi_num; $raw_i++){
		print OUTPUT "$em_corres_resi_num[$raw_i]\t$em_pdb_corres_resi_num[$raw_i]\t$raw_corres_resi_num[$raw_i]\n";
	}
	close OUTPUT;

	open OUTPUT, ">$ARGV[0]_resi_CN_cnc.txt" or die "can not open 4!\n" ;
	print OUTPUT "C\tN\n";
	for ($raw_i = 0; $raw_i < @resi_num_C_cnc; $raw_i++){
		print OUTPUT "$resi_num_C_cnc[$raw_i]\t$resi_num_N_cnc[$raw_i]\n";
	}
	print "resi_dele_itv_list: @resi_dele_itv_list\n";
	for ($resi_i = 1; $resi_i < @resi_dele_itv_list; $resi_i+=2){
		print OUTPUT "$resi_dele_itv_list[$resi_i-1]\t$resi_dele_itv_list[$resi_i]\n";
	}
	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

	#open OUTPUT, ">5resave_renumber.pml" or die "can not open 3!\n";
	#print OUTPUT "load $ARGV[0].pdb\n";
	#print OUTPUT "save $ARGV[0]_atomre.pdb, $ARGV[0]\n";
	#close OUTPUT;
	#system("pymol -qc 5resave_renumber.pml");

	open INPUT, "./$ARGV[0].pdb" or die "can not open!\n";
	open OUTPUT, ">./$ARGV[0]_atomre.pdb" or die "can not create!\n";
	$atom_c = 0;
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
		$atom_num = undef;
		for ($gezis_i = 6; $gezis_i <= 10 ; $gezis_i++){
			$atom_num .= $gezis[$gezis_i];
		}
		$atom_name = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			$atom_name .= $gezis[$gezis_i];
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
			$resi_name .= $gezis[$gezis_i];
		}
		$part_A = undef;
		for ($gezis_i = 11; $gezis_i < @gezis ; $gezis_i++){
			$part_A .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$atom_c++;
			printf OUTPUT "ATOM  %5d$part_A\n", $atom_c;
		}
		if ($line=~/^TER.*/){
			$atom_c++;
			printf OUTPUT "TER   %5d      $last_resi_name $last_chain_ID%4s\n",$atom_c, $last_resi_num ;
		}
		$last_resi_name = $resi_name;
		$last_resi_num = $resi_num;
		$last_chain_ID = $gezis[21];
	}
	close OUTPUT;
	close INPUT;




	@lig_pep_list = ();
	open INPUT, "./4lig_pep_list.txt" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		push @lig_pep_list, $line;
	}
	close INPUT;

	open INPUT, "./$ARGV[0]_atomre.pdb" or die "can not open!\n";
	open OUTPUT, ">$ARGV[0]_ssinfo.pdb" or die "can not create!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			print OUTPUT "$line\n";
		}
		if ($atom_mark eq "TER   "){
			print OUTPUT "$line\n";
		}

	}
	close INPUT;


	@cnb_a = ();
	@cnb_b = ();
	for ($lig_i = 0; $lig_i < @lig_pep_list; $lig_i++){

		open INPUT, "$lig_pep_list[$lig_i].pdb" or die "can not open!\n";
		@resi_num_list = ();
		@resi_name_list = ();
		@atom_num_list = ();
		@atom_name_list = ();
		@cnc_a_list = ();
		@cnc_b_list = ();
       		while(chomp($line=<INPUT>)){
                        @gezis = split //, $line;
                        $atom_mark = undef;
                        for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
                        	$atom_mark .= "$gezis[$gezis_i]";
			}
			if (($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM")) {
				$resi_num = undef;
				for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
					$resi_num .= $gezis[$gezis_i];
				}
				push @resi_num_list, $resi_num;
                        
				$atom_num = undef;
				for ($gezis_i = 6; $gezis_i <= 10 ; $gezis_i++){
					$atom_num .= $gezis[$gezis_i];
				}
				push @atom_num_list, $atom_num;
                        
				$atom_name = undef;
				for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
					$atom_name .= $gezis[$gezis_i];
				}
				push @atom_name_list, $atom_name;
                        
				$resi_name = undef;
				for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
					$resi_name .= $gezis[$gezis_i];
				}
				push @resi_name_list, $resi_name;
			}
			if ($atom_mark eq "CONECT") {
				$cnc_a = undef;
				for ($gezis_i = 6; $gezis_i <= 10; $gezis_i++){
					$cnc_a .= "$gezis[$gezis_i]";
				}
				push @cnc_a_list, $cnc_a;
				$cnc_b = undef;
				for ($gezis_i = 11; $gezis_i <= 15; $gezis_i++){
					$cnc_b .= "$gezis[$gezis_i]";
				}
				push @cnc_b_list, $cnc_b;
			}
                }
		close INPUT;

		@resi_num_cnc_a_list = ();
		@resi_name_cnc_a_list = ();
		@atom_name_cnc_a_list = ();
		@resi_num_cnc_b_list = ();
		@resi_name_cnc_b_list = ();
		@atom_name_cnc_b_list = ();
		for ($cnc_i = 0; $cnc_i < @cnc_a_list; $cnc_i++){
			for ($atom_i = 0; $atom_i < @atom_num_list; $atom_i++){
				if ($cnc_a_list[$cnc_i] == $atom_num_list[$atom_i]){
					push @resi_num_cnc_a_list, $resi_num_list[$atom_i];
					push @resi_name_cnc_a_list, $resi_name_list[$atom_i];
					push @atom_name_cnc_a_list, $atom_name_list[$atom_i];
				}
				if ($cnc_b_list[$cnc_i] == $atom_num_list[$atom_i]){
					push @resi_num_cnc_b_list, $resi_num_list[$atom_i];
					push @resi_name_cnc_b_list, $resi_name_list[$atom_i];
					push @atom_name_cnc_b_list, $atom_name_list[$atom_i];
				}
			}
		}



		open INPUT, "$lig_pep_list[$lig_i]_align.pdb" or die "can not open!\n";
		@resi_num_list = ();
		@resi_name_list = ();
		@atom_num_list = ();
		@atom_name_list = ();
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$atom_mark .= $gezis[$gezis_i];
			}
			if (($atom_mark eq "ATOM  ") || ($atom_mark eq "HETATM")) {
				$part_A = undef;
				for ($gezis_i = 11; $gezis_i < @gezis ; $gezis_i++){
					$part_A .= $gezis[$gezis_i];
				}
				$atom_c++;
				printf OUTPUT "$atom_mark%5d$part_A\n", $atom_c;

				$resi_num = undef;
				for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
					$resi_num .= $gezis[$gezis_i];
				}
				push @resi_num_list, $resi_num;
                        
				$atom_num = undef;
				for ($gezis_i = 6; $gezis_i <= 10 ; $gezis_i++){
					$atom_num .= $gezis[$gezis_i];
				}
				push @atom_num_list, $atom_num;
                        
				$atom_name = undef;
				for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
					$atom_name .= $gezis[$gezis_i];
				}
				push @atom_name_list, $atom_name;
                        
				$resi_name = undef;
				for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
					$resi_name .= $gezis[$gezis_i];
				}
				push @resi_name_list, $resi_name;
			}
                }
		close INPUT;
		print OUTPUT "TER\n";

		for ($cnc_i = 0; $cnc_i < @resi_num_cnc_a_list; $cnc_i++){
			$match_ss_mark = 0;
			for ($atom_i = 0; $atom_i < @atom_num_list; $atom_i++){
				if (($resi_num_cnc_a_list[$cnc_i] == $resi_num_list[$atom_i]) &&  ($resi_name_cnc_a_list[$cnc_i] eq $resi_name_list[$atom_i])  &&  ($atom_name_cnc_a_list[$cnc_i] eq $atom_name_list[$atom_i]) ){
					$match_ss_mark = 1;
					push @cnb_a, $atom_num_list[$atom_i];
					last;
				}
			}
			if ($match_ss_mark == 0){
				die "lig or pep structure error!\n";
			}
		}
		for ($cnc_i = 0; $cnc_i < @resi_num_cnc_b_list; $cnc_i++){
			$match_ss_mark = 0;
			for ($atom_i = 0; $atom_i < @atom_num_list; $atom_i++){
				if (($resi_num_cnc_b_list[$cnc_i] == $resi_num_list[$atom_i]) &&  ($resi_name_cnc_b_list[$cnc_i] eq $resi_name_list[$atom_i])  &&  ($atom_name_cnc_b_list[$cnc_i] eq $atom_name_list[$atom_i]) ){
					$match_ss_mark = 1;
					push @cnb_b, $atom_num_list[$atom_i];
					last;
				}
			}
			if ($match_ss_mark == 0){
				die "lig or pep structure error!\n";
			}
		}



	}


	open INPUT, "./$ARGV[0]_atomre.pdb" or die "can not open!\n";
	@atomre_lines = ();
	while(chomp($line=<INPUT>)){
		push @atomre_lines, $line;
	}
	close INPUT;
	$C_N_dist_max = 1.6;  #this is 1.3 measured by pymol. 
	open OUTPUT1, ">7ssinfo_match.txt" or die "can not create!\n";

	$domain_count = 0;
	$line_ini = 0;
	$itv_a = 0;
	for ($ag_i = 2; $ag_i < @ARGV; $ag_i++){
		@resi_num_list_2 = (); 
		@atom_num_list_2 = (); 
		@atom_name_list_2 = ();
		@chain_ID_list_2 = (); 
		@pos_x_list_2 = ();
		@pos_y_list_2 = ();
		@pos_z_list_2 = ();
		$crt_resi_num = -99999;
		$resi_count = 0;
		for ($line_i = $line_ini; $line_i < @atomre_lines; $line_i++){
			
			@gezis = split //, $atomre_lines[$line_i];
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$atom_mark .= $gezis[$gezis_i];
			}
			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			if ($resi_num > $crt_resi_num){
				$resi_count++;
				$crt_resi_num = $resi_num;
			}
			$atom_num = undef;
			for ($gezis_i = 6; $gezis_i <= 10 ; $gezis_i++){
				$atom_num .= $gezis[$gezis_i];
			}
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
			if ($atom_mark eq "ATOM  "){
				push @resi_num_list_2, $resi_num;
				push @atom_num_list_2, $atom_num;
				push @atom_name_list_2, $atom_name;
				push @chain_ID_list_2, $gezis[21];
				push @pos_x_list_2, $pos_x;
				push @pos_y_list_2, $pos_y;
				push @pos_z_list_2, $pos_z;
			}
			if ($atom_mark eq "TER   "){
				$domain_count++;
				$line_ini = $line_i+1;
				last;
			}
		}
		print OUTPUT1 "domain $domain_count\nresi_count: $resi_count\nlabel: $chain_$ARGV[$ag_i]\n";
		#$itv_a+=$resi_count;
		

		@ss_a = ();
		@ss_b = ();
		@ss_label = ();
		@conect_list = ();
		@conect_label_list = ();
		@conect_list_2 = ();
		@conect_label_list_2 = ();
		$em_seq = undef;
		$md_seq = undef;
		open INPUT, "./chain_$ARGV[$ag_i]/chain_$ARGV[$ag_i]_mm.aln" or die "can not open 2!\n";
		while(chomp($line=<INPUT>)){
			@re_item = split /\s+/, $line;
			if ($re_item[0] eq "chain_$ARGV[$ag_i]_em"){
				$em_seq.=$re_item[1];
			}
			if ($re_item[0] eq "chain_$ARGV[$ag_i]_md"){
				$md_seq.=$re_item[1];
			}
		}
        
		@md_seq_list = split //, $md_seq;
		@em_seq_list = split //, $em_seq;
		@fix_region = ();
		for ($em_i = 0; $em_i < @em_seq_list; $em_i++){
			if ($em_seq_list[$em_i] eq "-") {
				push @fix_region, $em_i + 1;
			}
		}
		close INPUT; 
		
		$model_resi_count = @md_seq_list;
		print OUTPUT1 "model_resi_count: $model_resi_count\nresi_count: $resi_count\n";
		print "fix_region: @fix_region\n";
		print OUTPUT1 "fix_region: @fix_region\n";
        
		#open INPUT, "./chain_$ARGV[$ag_i]/chain_$ARGV[$ag_i]_itv.txt" or die "can not open 5!\n";
		#while(chomp($line=<INPUT>)){
		#	@itv_input = split /\s+/, $line;
		#	$itv_a = $itv_input[0];
		#}
		#close INPUT;

		push @conect_list, $fix_region[0] + $itv_a;
		push @conect_label_list, $ARGV[$ag_i];
		for ($fix_i = 1; $fix_i < @fix_region; $fix_i++){
			if ($fix_region[$fix_i] > $fix_region[$fix_i-1] + 1){
				push @conect_list, $fix_region[$fix_i-1] + $itv_a;
				push @conect_label_list, $ARGV[$ag_i];
				push @conect_list, $fix_region[$fix_i] + $itv_a;
				push @conect_label_list, $ARGV[$ag_i];
			}
		}
		push @conect_list, $fix_region[-1] + $itv_a;
		push @conect_label_list, $ARGV[$ag_i];
		print OUTPUT1 "fix region itv : @conect_list\n";
		
		$input_cnc_file = "./chain_$ARGV[$ag_i]/chain_$ARGV[$ag_i]"."_resi_CN_cnc.txt";
		open INPUT, "$input_cnc_file" or die "can not open 4!\n" ;
		$line_c = 0;
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 1){
				@cncs = split /\s+/, $line;
				push @conect_list_2, $cncs[0] + $itv_a;
				push @conect_label_list_2, $ARGV[$ag_i];
				push @conect_list_2, $cncs[1] + $itv_a;
				push @conect_label_list_2, $ARGV[$ag_i];
			}
		}
		close INPUT;
		print OUTPUT1 "CN cnc : @conect_list_2\n";


		open INPUT, "./chain_$ARGV[$ag_i]/ssbond_filter_model_$ARGV[$ag_i].txt" or die "can not open!\n";
		$line_c = 0;
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 0){
				@ss_info = split /\s+/, $line;
				push @ss_a, $ss_info[4] + $itv_a;
				push @ss_b, $ss_info[7] + $itv_a;
				push @ss_label, $ss_info[3];
			}
			
		}
		close INPUT;
		print OUTPUT1 "s-s bond a: @ss_a\n";
		print OUTPUT1 "s-s bond b: @ss_b\n";

		print "@ss_a\n@ss_b\n@ss_label\n";
		print "@conect_list\n@conect_label_list\n";
		print "@conect_list_2\n@conect_label_list_2\n";


		for ($cne_i = 0; $cne_i <= 1; $cne_i+=2){
			@conect_pair_b = ();
			$mark_3 = 0;
			$mark_4 = 0;
			$mark_5 = 0;
			$C_atom_idx = undef;
			$N_atom_idx = undef;
			$cndist = -9999;
			for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
				if (($conect_list[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ") && ($conect_label_list[$cne_i+1] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_3 = 1;
					$C_atom_idx = $ret_i;
					print OUTPUT1 "$pos_x_list_2[$ret_i]\t$pos_y_list_2[$ret_i]\t$pos_z_list_2[$ret_i]\n";
				}
				if (($conect_list[$cne_i + 1] + 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")  && ($conect_label_list[$cne_i+1] eq $chain_ID_list_2[$ret_i])  ){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_4 = 1;
					$N_atom_idx = $ret_i;
					print OUTPUT1 "$pos_x_list_2[$ret_i]\t$pos_y_list_2[$ret_i]\t$pos_z_list_2[$ret_i]\n";
				}
				if (($mark_3 == 1)  && ($mark_4 == 1) ){
					$cndist = sqrt(($pos_x_list_2[$C_atom_idx] - $pos_x_list_2[$N_atom_idx])**2 + ($pos_y_list_2[$C_atom_idx] - $pos_y_list_2[$N_atom_idx])**2 + ($pos_z_list_2[$C_atom_idx] - $pos_z_list_2[$N_atom_idx])**2) ;
					if ($cndist > $C_N_dist_max ){
						$mark_5 = 1;
					}
				}
			}
			if (($mark_3 == 1)  && ($mark_4 == 1)  && ($mark_5 == 1) ){
				print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
				print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]! $cndist $C_N_dist_max $C_atom_idx $N_atom_idx\n";
				printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
				}
			else{
				print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]! mark : $mark_3 $mark_4 $mark_5\n";
				print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]! mark : $mark_3 $mark_4 $mark_5\n";
			}
		}
		for ($cne_i = @conect_list-2; $cne_i < @conect_list; $cne_i+=2){
			@conect_pair_a = ();
			$mark_1 = 0;
			$mark_2 = 0;
			$mark_5 = 0;
			$C_atom_idx = undef;
			$N_atom_idx = undef;
			$cndist = -9999;
			for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
				if (($conect_list[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ") && ($conect_label_list[$cne_i] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_a, $atom_num_list_2[$ret_i];
					$mark_1 = 1;
					$N_atom_idx = $ret_i;
				}
				if (($conect_list[$cne_i] - 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ") && ($conect_label_list[$cne_i] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_a, $atom_num_list_2[$ret_i];
					$mark_2 = 1;
					$C_atom_idx = $ret_i;
				}
				if (($mark_1 == 1)  && ($mark_2 == 1) ){
					$cndist = sqrt(($pos_x_list_2[$C_atom_idx] - $pos_x_list_2[$N_atom_idx])**2 + ($pos_y_list_2[$C_atom_idx] - $pos_y_list_2[$N_atom_idx])**2 + ($pos_z_list_2[$C_atom_idx] - $pos_z_list_2[$N_atom_idx])**2) ;
					if ($cndist > $C_N_dist_max ){
						$mark_5 = 1;
					}
				}
			}
			if (($mark_1 == 1) && ($mark_2 == 1)   && ($mark_5 == 1) ){
				print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
				print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
				printf OUTPUT "CONECT%5d%5d\n", $conect_pair_a[0],$conect_pair_a[1];
				}
			else{
				print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]! mark : $mark_1 $mark_2 $mark_5\n";
				print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]! mark : $mark_1 $mark_2 $mark_5\n";
			}
		}
		for ($cne_i = 2; $cne_i < @conect_list-2; $cne_i+=2){
			@conect_pair_a = ();
			@conect_pair_b = ();
			$mark_1 = 0;
			$mark_2 = 0;
			$mark_3 = 0;
			$mark_4 = 0;
			$mark_5 = 0;
			$mark_6 = 0;
			$C_atom_idx_a = undef;
			$N_atom_idx_a = undef;
			$C_atom_idx_b = undef;
			$N_atom_idx_b = undef;
			$cndist_1 = -9999;
			$cndist_2 = -9999;
			for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
				if (($conect_list[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ") && ($conect_label_list[$cne_i] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_a, $atom_num_list_2[$ret_i];
					$mark_1 = 1;
					$N_atom_idx_a = $ret_i;
				}
				if (($conect_list[$cne_i] - 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ") && ($conect_label_list[$cne_i] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_a, $atom_num_list_2[$ret_i];
					$mark_2 = 1;
					$C_atom_idx_a = $ret_i;
				}
				if (($mark_1 == 1)  && ($mark_2 == 1) ){
					$cndist_1 = sqrt(($pos_x_list_2[$C_atom_idx_a] - $pos_x_list_2[$N_atom_idx_a])**2 + ($pos_y_list_2[$C_atom_idx_a] - $pos_y_list_2[$N_atom_idx_a])**2 + ($pos_z_list_2[$C_atom_idx_a] - $pos_z_list_2[$N_atom_idx_a])**2) ;
					if ($cndist_1 > $C_N_dist_max ){
						$mark_5 = 1;
					}
				}
				if (($conect_list[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ") && ($conect_label_list[$cne_i+1] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_3 = 1;
					$C_atom_idx_b = $ret_i;
				}
				if (($conect_list[$cne_i + 1] + 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ") && ($conect_label_list[$cne_i+1] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_4 = 1;
					$N_atom_idx_b = $ret_i;
				}
				if (($mark_3 == 1)  && ($mark_4 == 1) ){
					$cndist_2 = sqrt(($pos_x_list_2[$C_atom_idx_b] - $pos_x_list_2[$N_atom_idx_b])**2 + ($pos_y_list_2[$C_atom_idx_b] - $pos_y_list_2[$N_atom_idx_b])**2 + ($pos_z_list_2[$C_atom_idx_b] - $pos_z_list_2[$N_atom_idx_b])**2) ;
					if ($cndist_2 > $C_N_dist_max ){
						$mark_6 = 1;
					}
				}
			}
			if (($mark_1 == 1) && ($mark_2 == 1)  && ($mark_5 == 1)){
				print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i]-1!\n";
				print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i]-1!\n";
				printf OUTPUT "CONECT%5d%5d\n", $conect_pair_a[0],$conect_pair_a[1];
				}
			else{
				print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i]-1! mark : $mark_1 $mark_2 $mark_5\n\n";
				print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i]-1! mark : $mark_1 $mark_2 $mark_5\n\n";
			}
			if (($mark_3 == 1) && ($mark_4 == 1)  && ($mark_6 == 1)){
				print "loop terminals match: $conect_list[$cne_i+1]  $conect_list[$cne_i+1]+1!\n";
				print OUTPUT1 "loop terminals match: $conect_list[$cne_i+1]  $conect_list[$cne_i+1]+1!\n";
				printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
				}
			else{
				print "NO loop terminals match: $conect_list[$cne_i+1]  $conect_list[$cne_i+1]+1! mark : $mark_3 $mark_4 $mark_6\n\n";
				print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i+1]  $conect_list[$cne_i+1]+1! mark : $mark_3 $mark_4 $mark_6\n\n";
			}
		}
		for ($cne_i = 0; $cne_i < @conect_list_2; $cne_i+=2){
			@conect_pair_b = ();
			$mark_3 = 0;
			$mark_4 = 0;
			$mark_5 = 0;
			$C_atom_idx = undef;
			$N_atom_idx = undef;
			$cndist = -9999;
			for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
				if (($conect_list_2[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ") && ($conect_label_list[$cne_i] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_3 = 1;
					$C_atom_idx = $ret_i;
				}
				if (($conect_list_2[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ") && ($conect_label_list[$cne_i+1] eq $chain_ID_list_2[$ret_i])){
					push @conect_pair_b, $atom_num_list_2[$ret_i];
					$mark_4 = 1;
					$N_atom_idx = $ret_i;
				}
				if (($mark_3 == 1)  && ($mark_4 == 1) ){
					$cndist = sqrt(($pos_x_list_2[$C_atom_idx] - $pos_x_list_2[$N_atom_idx])**2 + ($pos_y_list_2[$C_atom_idx] - $pos_y_list_2[$N_atom_idx])**2 + ($pos_z_list_2[$C_atom_idx] - $pos_z_list_2[$N_atom_idx])**2) ;
					if ($cndist > $C_N_dist_max ){
						$mark_5 = 1;
					}
				}
			}
			if (($mark_3 == 1)  && ($mark_4 == 1)   && ($mark_5 == 1) ){
				print "loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
				print OUTPUT1 "loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
				printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
				}
			else{
				print "NO loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]! mark : $mark_3 $mark_4 $mark_5\n\n";
				print OUTPUT1 "NO loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]! mark : $mark_3 $mark_4 $mark_5\n\n";
			}
		}
        
        
        
        
		#@ssbond_list_w = @ssbond_list;
		#open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not open!\n";
		for ($ss_i = 0; $ss_i < @ss_a; $ss_i++){
			@ss_pair = ();
			$mark_1 = 0;
			$mark_2 = 0;
			for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
				if (($ss_a[$ss_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " SG ") && ($chain_ID_list_2[$ret_i] eq $ss_label[$ss_i]) ){
					push @ss_pair, $atom_num_list_2[$ret_i];
					$mark_1 = 1;
				}
				if (($ss_b[$ss_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " SG ") && ($chain_ID_list_2[$ret_i] eq $ss_label[$ss_i]) ){
					push @ss_pair, $atom_num_list_2[$ret_i];
					$mark_2 = 1;
				}
			}
			if (($mark_1 == 1) && ($mark_2 == 1) ){
				print "s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
				print OUTPUT1 "s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
				printf OUTPUT "CONECT%5d%5d\n", $ss_pair[0],$ss_pair[1];
			}
			else{
				print "NO s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]! mark : $mark_1 $mark_2 \n";
				print OUTPUT1 "NO s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]! mark : $mark_1 $mark_2 \n";
			}
			
		}
        
	}


	for ($cnc_i = 0; $cnc_i < @cnb_b; $cnc_i++){
		printf OUTPUT "CONECT%5d%5d\n", $cnb_a[$cnc_i],  $cnb_b[$cnc_i];
	}
	close OUTPUT;
	close OUTPUT1;

	#system("\$SCHRODINGER/utilities/prepwizard -rehtreat -disulfides -propka_pH 7.0 -rmsd 2 -fillsidechains $ARGV[0]_ssinfo.pdb $ARGV[0]_ssinfo_prep.pdb");
	system("$ARGV[1]/utilities/prepwizard -rehtreat -disulfides -propka_pH 7.0 -rmsd 20 -fillsidechains $ARGV[0]_ssinfo.pdb $ARGV[0]_ssinfo_prep.pdb");

	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

	@items = split /_/, $ARGV[0];

	system("../pdb2fasta chain_$items[1]_fe.B99990001.pdb > chain_$items[1]_fe_model.fasta");

	open INPUT, "chain_$items[1]_fe_model.fasta" or die "can not open 1!\n";
	open OUTPUT, ">chain_$items[1]_fe_model_2.fasta" or die "can not create!\n";
	$line_c = 0;
	print OUTPUT ">chain_$items[1]_md\n";
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c > 1){
			print OUTPUT "$line\n";
		}
	}
	close INPUT;
	close OUTPUT;

	system("cat chain_$items[1]_em_2.fasta chain_$items[1]_fe_model_2.fasta > chain_$items[1]_mm.fasta");
	system("../clustalw2 chain_$items[1]_mm.fasta");
	
	#start fix the aln file
	$em_seq = undef;
	$tm_seq = undef;
	$cth_seq = undef;

	$seq_mark = 0;
	system("mv chain_$items[1]_mm.aln chain_$items[1]_mm.aln_bak");
	open INPUT, "chain_$items[1]_mm.aln_bak" or die "can not open!\n";
	open OUTPUT, ">chain_$items[1]_mm.aln" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@re_item = split /\s+/, $line;
		if ($re_item[0] eq "chain_$items[1]_em"){
			$seq_mark++;
			$em_seq.=$re_item[1];

			$wp_mark = 0;
			@gezis = split //, $line;
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
			next;
		}
		if ($re_item[0] eq "chain_$items[1]_md"){
			$seq_mark++;
			$tm_seq.=$re_item[1];
			next;
		}
		if ($seq_mark == 2){
			@gezis = split //, $line;
			for ($gezis_i = $cth_start_num; $gezis_i < $cth_start_num+60; $gezis_i++){
				$cth_seq.=$gezis[$gezis_i];
			}
			$seq_mark = 0;
		}
	}
	
	@em_seq_list = split //, $em_seq;
	@tm_seq_list = split //, $tm_seq;
	@cth_seq_list = split //, $cth_seq;
	for ($em_i = 0; $em_i < @em_seq_list; $em_i++){
		if ($em_seq_list[$em_i] ne "-"){
			$start_index = $em_i ;
			last;
		}
	}
	for ($em_i = -1; $em_i > -@em_seq_list; $em_i--){
		if ($em_seq_list[$em_i] ne "-"){
			$end_index = @em_seq_list + ($em_i + 1) - 1;
			last;
		}
	}

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


	open INPUT, "chain_$items[1].pdb" or die "can not open 4!\n" ;
	@C_pos_x=();
	@C_pos_y=();
	@C_pos_z=();
	@N_pos_x=();
	@N_pos_y=();
	@N_pos_z=();
	@resi_num_list = ();
	@resi_name_list = ();
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
			push @resi_num_list, $resi_num;
			push @resi_name_list, $resi_name;
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

	$input_cnc_file = "chain_$items[1]"."_resi_CN_cnc.txt";
	open INPUT, "$input_cnc_file" or die "can not open 4!\n" ;
	$line_c = 0;
	@C_cnc_list = ();
	@N_cnc_list = ();
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c > 1){
			@cncs = split /\s+/, $line;
			push @C_cnc_list, $cncs[0];
			push @N_cnc_list, $cncs[1];
		}
	}
	close INPUT;





	$C_N_dist_min = 1.4;  #this is 1.3 measured by pymol. 
	$itv_start_mark = 0;
	$itv_end_mark = 1;
	$itv_mark = 0;
	$em_resi_count=0;
	for ($tm_i = $start_index; $tm_i <= $end_index; $tm_i++){
		if ($em_seq_list[$tm_i] ne "-") {
			$em_resi_count++;
		}
		if (($em_seq_list[$tm_i] eq "-") && ($itv_start_mark == 0) && ($itv_end_mark == 1)) {
			#$em_itv_start_b = $em_seq_list[$tm_i];
			$tm_itv_start = $tm_seq_list[$tm_i];
			$tm_itv_start_idx = $tm_i;
			$em_itv_start_a = $em_seq_list[$tm_i-1];
			$em_itv_start_motif = $em_seq_list[$tm_i-3].$em_seq_list[$tm_i-2].$em_seq_list[$tm_i-1];
			$tm_itv_start_motif = $tm_seq_list[$tm_i-3].$tm_seq_list[$tm_i-2].$tm_seq_list[$tm_i-1];
			$itv_start_mark = 1;
			$itv_end_mark = 0;
			$itv_mark = 1;
		}
		elsif (($em_seq_list[$tm_i] ne "-") && ($itv_start_mark == 1) && ($itv_end_mark == 0)) {
			$em_itv_end_b = $em_seq_list[$tm_i];
			$tm_itv_end = $tm_seq_list[$tm_i-1];
			$tm_itv_end_idx = $tm_i -1;
			$em_itv_end_motif = $em_seq_list[$tm_i].$em_seq_list[$tm_i+1].$em_seq_list[$tm_i+2];
			$tm_itv_end_motif = $tm_seq_list[$tm_i].$tm_seq_list[$tm_i+1].$tm_seq_list[$tm_i+2];
			#$em_itv_end_a = $em_seq_list[$tm_i-1];
			$itv_start_mark = 0;
			$itv_end_mark = 1;
			if (($itv_mark == 1) && ($em_itv_end_b eq $tm_itv_start) && ($em_itv_end_motif ne $tm_itv_end_motif) ){
				$itv_mark = 0;
				$em_seq_list[$tm_itv_start_idx] = $em_itv_end_b;
				$em_seq_list[$tm_itv_end_idx+1] = "-";
				$cth_seq_list[$tm_itv_start_idx] = "\*";
				$cth_seq_list[$tm_itv_end_idx+1] = " ";
			}
			elsif (($itv_mark == 1) && ($em_itv_start_a eq $tm_itv_end) && ($em_itv_start_motif ne $tm_itv_start_motif) ){
				$itv_mark = 0;
				$em_seq_list[$tm_itv_end_idx] = $em_itv_start_a;
				$em_seq_list[$tm_itv_start_idx-1] = "-";
				$cth_seq_list[$tm_itv_end_idx] = "\*";
				$cth_seq_list[$tm_itv_start_idx-1] = " ";
			}
			elsif (($itv_mark == 1) && ($em_itv_end_b eq $tm_itv_start) && ($em_itv_end_motif eq $tm_itv_end_motif) ){
				if ($resi_name_list[$em_resi_count-1] ne $aa_name{$em_itv_end_b}){
					print "warning: aln resi anme not same as pdb resi name!\n";
				}
				$cnc_match=0;
				for ($cnc_i = 0; $cnc_i< @C_cnc_list; $cnc_i++){
					if (($C_cnc_list[$cnc_i] == $resi_num_list[$em_resi_count-2]) && ($N_cnc_list[$cnc_i] == $resi_num_list[$em_resi_count-1])){
						$cnc_match = 1;
						last;
					}
				}
				$cndist = sqrt(($C_pos_x[$em_resi_count-1] - $N_pos_x[$em_resi_count])**2 + ($C_pos_y[$em_resi_count-1] - $N_pos_y[$em_resi_count])**2 + ($C_pos_z[$em_resi_count-1] - $N_pos_z[$em_resi_count])**2) ;
				$itv_mark = 0;
				if (($cndist<$C_N_dist_min) || ($cnc_match==1))   {
					$em_seq_list[$tm_itv_start_idx] = $em_itv_end_b;
					$em_seq_list[$tm_itv_end_idx+1] = "-";
					$cth_seq_list[$tm_itv_start_idx] = "\*";
					$cth_seq_list[$tm_itv_end_idx+1] = " ";
				}
			}
			elsif (($itv_mark == 1) && ($em_itv_start_a eq $tm_itv_end) && ($em_itv_start_motif eq $tm_itv_start_motif) ){
				if ($resi_name_list[$em_resi_count-1] ne $aa_name{$em_itv_end_b}){
					print "warning: aln resi anme not same as pdb resi name!\n";
				}
				$cnc_match=0;
				for ($cnc_i = 0; $cnc_i< @C_cnc_list; $cnc_i++){
					if (($C_cnc_list[$cnc_i] == $resi_num_list[$em_resi_count-2]) && ($N_cnc_list[$cnc_i] == $resi_num_list[$em_resi_count-1])){
						$cnc_match = 1;
						last;
					}
				}
				$cndist = sqrt(($C_pos_x[$em_resi_count-1] - $N_pos_x[$em_resi_count])**2 + ($C_pos_y[$em_resi_count-1] - $N_pos_y[$em_resi_count])**2 + ($C_pos_z[$em_resi_count-1] - $N_pos_z[$em_resi_count])**2) ;
				$itv_mark = 0;
				if (($cndist<$C_N_dist_min) || ($cnc_match==1))   {
					$em_seq_list[$tm_itv_end_idx] = $em_itv_start_a;
					$em_seq_list[$tm_itv_start_idx-1] = "-";
					$cth_seq_list[$tm_itv_end_idx] = "\*";
					$cth_seq_list[$tm_itv_start_idx-1] = " ";
				}
			}
		}
	}




	print OUTPUT "CLUSTAL 2.1 multiple sequence alignment


";	
	$resi_count=0;
	$tm_line=undef;
	$em_line=undef;
	$cth_line=undef;
	for ($tm_i = 0; $tm_i < @tm_seq_list; $tm_i++){
		$resi_count++;
		$tm_line.=$tm_seq_list[$tm_i];
		$em_line.=$em_seq_list[$tm_i];
		$cth_line.=$cth_seq_list[$tm_i];
		if ($resi_count==60){
			print OUTPUT "chain_$items[1]_em      $em_line\n";
			print OUTPUT "chain_$items[1]_md      $tm_line\n";
			print OUTPUT "                $cth_line\n\n";
			$resi_count = 0;
			$tm_line=undef;
			$em_line=undef;
			$cth_line=undef;
		}
	}
	if (($resi_count<60) && ($resi_count>0)){
		print OUTPUT "chain_$items[1]_em      $em_line\n";
		print OUTPUT "chain_$items[1]_md      $tm_line\n";
		print OUTPUT "                $cth_line\n";
		$resi_count = 0;
		$tm_line=undef;
		$em_line=undef;
		$cth_line=undef;
	}


	close INPUT;
	close OUTPUT;

	#end of the fix of aln file




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
	$start_index = 0;
	$end_index = @md_seq_list - 1;
	for ($md_i = 0; $md_i < @md_seq_list; $md_i++){
		if ($md_seq_list[$md_i] ne "-"){
			$start_index = $md_i ;
			last;
		}
	}
	for ($md_i = -1; $md_i > -@md_seq_list; $md_i--){
		if ($md_seq_list[$md_i] ne "-"){
			$end_index = @md_seq_list + ($md_i + 1) - 1;
			last;
		}
	}

	for ($md_i = $start_index; $md_i <= $end_index; $md_i++){
		if ($md_seq_list[$md_i] eq "-"){
			print "the model may miss resiues!\n";
			die;
		}
	}
	#print "start index: $start_index\nend_index: $end_index\n";
	#print OUTPUT "\n";
	#close OUTPUT;
	close INPUT;


	open INPUT, "chain_$items[1].pdb" or die "can not open 3!\n";
	$line_c = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$line_c++;
			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			if ($line_c == 1){
				$start_seq_num = $resi_num;
			}
			else{
				$end_seq_num = $resi_num;
			}
		}
	}
	close INPUT;

	$seq_len = $end_seq_num - ($start_seq_num - 1);
	
	open OUTPUT, ">error_info_9.txt" or die "can not create!\n";
	if ($seq_len != @em_seq_list){
		print "the residue serial num in original pdb may be wrong!\n ";
		print OUTPUT "the residue serial num in original pdb may be wrong!\n ";
		#die;
	}
	close OUTPUT;

	open OUTPUT, ">chain_$items[1]_itv.txt" or die "can not create!\n";
	#$itv_seq = $start_seq_num + $start_index - 1;
	$itv_seq=0;
	print OUTPUT "$itv_seq\t$seq_len\n";
	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

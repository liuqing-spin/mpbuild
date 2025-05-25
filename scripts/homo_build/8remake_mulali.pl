	#input info: input_for5
	
	@items = split /_/, $ARGV[2];

	system("../pdb2fasta $ARGV[0]_fit.pdb > $ARGV[0]_fit.fasta");

	open INPUT, "$ARGV[0]_fit.fasta" or die "can not open!\n";
	open OUTPUT, ">$ARGV[0]_fit_2.fasta" or die "can not open!\n";
	$line_c = 0;
	print OUTPUT ">$ARGV[0]\n";
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c > 1){
			print OUTPUT "$line\n";
		}
	}
	close INPUT;
	close OUTPUT;

	system("cat $ARGV[0]_fit_2.fasta $ARGV[2]_em_2.fasta  > template_em.fasta ");
	system("../clustalw2 template_em.fasta");
	system("mv template_salign.ali template_salign.ali_bak");

	open INPUT, "template_salign.ali_bak" or die "can not open!\n";
	$line_c = 0;
	while(chomp($line=<INPUT>)){
		if ($line=~">P1.*"){
			$line_c++;
			push @line_forali, $line;
			next;
		}
		if ($line_c==1){
			push @line_forali, $line;
			$line_c = 0;
			next;
		}
	}

	close INPUT;


	$em_seq = undef;
	$tm_seq = undef;
	$cth_seq = undef;

	$seq_mark = 0;
	system("mv template_em.aln template_em.aln_bak");
	open INPUT, "template_em.aln_bak" or die "can not open!\n";
	open OUTPUT, ">template_em.aln" or die "can not open!\n";
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
		if ($re_item[0] eq "$ARGV[0]"){
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


	open INPUT, "$ARGV[2].pdb" or die "can not open 4!\n" ;
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

	$input_cnc_file = "$ARGV[2]"."_resi_CN_cnc.txt";
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
			print OUTPUT "$ARGV[0]            $tm_line\n";
			print OUTPUT "chain_$items[1]_em      $em_line\n";
			print OUTPUT "                $cth_line\n\n";
			$resi_count = 0;
			$tm_line=undef;
			$em_line=undef;
			$cth_line=undef;
		}
	}
	if (($resi_count<60) && ($resi_count>0)){
		print OUTPUT "$ARGV[0]            $tm_line\n";
		print OUTPUT "chain_$items[1]_em      $em_line\n";
		print OUTPUT "                $cth_line\n";
		$resi_count = 0;
		$tm_line=undef;
		$em_line=undef;
		$cth_line=undef;
	}


	close INPUT;
	close OUTPUT;

	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";



	
	open OUTPUT, ">template_salign.ali" or die "can not open!\n";
	print OUTPUT "\n";
	for ($line_i = 0; $line_i < @line_forali; $line_i+=2){
		print OUTPUT "$line_forali[$line_i]\n";
		print OUTPUT "$line_forali[$line_i+1]\n";
		if ($line_forali[$line_i] =~ ">P1;$ARGV[0].*"){
			$gezis_c = 0;
			$end_mark = 0;
			for ($gezis_i = 0; $gezis_i < @tm_seq_list; $gezis_i++){
				$gezis_c++;
				if ($gezis_c == 75){
					if ($gezis_i == @tm_seq_list - 1){
						print OUTPUT "$tm_seq_list[$gezis_i]*\n";
						$end_mark = 1;
					}
					else{
						print OUTPUT "$tm_seq_list[$gezis_i]\n";
						$gezis_c = 0;
					}
				}
				else{
					print OUTPUT "$tm_seq_list[$gezis_i]";
				}
			}
			if ($end_mark == 0){
				print OUTPUT "*\n";
			}
			print OUTPUT "\n";
		}	
		elsif ($line_forali[$line_i] =~ ">P1;$ARGV[2]$items[1].*"){
			$gezis_c = 0;
			$end_mark = 0;
			for ($gezis_i = 0; $gezis_i < @em_seq_list; $gezis_i++){
				$gezis_c++;
				if ($gezis_c == 75){
					if ($gezis_i == @em_seq_list - 1){
						print OUTPUT "$em_seq_list[$gezis_i]*\n";
						$end_mark = 1;
					}
					else{
						print OUTPUT "$em_seq_list[$gezis_i]\n";
						$gezis_c = 0;
					}
				}
				else{
					print OUTPUT "$em_seq_list[$gezis_i]";
				}
			}
			if ($end_mark == 0){
				print OUTPUT "*\n";
			}
			#print OUTPUT "\n";
		}	
	}
	close OUTPUT;
	close INPUT;



	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

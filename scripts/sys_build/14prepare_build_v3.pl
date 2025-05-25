
	
	@cid_list = ();
	@lig_list = ();
	@pep_list = ();
	$wt_inh = undef;
	for ($argv_i = 0; $argv_i < @ARGV ; $argv_i++){
		if ($ARGV[$argv_i] eq "-h"){
			print "		-ci	the list of chain IDs.
	-lg	the list of ligand pdbs.
	-lp	the list of peptide pdbs.
	-wt	the water inhole signal, 1 or 0.
       	";
		die "	-h	print above information\n";
		}
		if ($ARGV[$argv_i] eq "-ci"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @cid_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-lg"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @lig_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-lp"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @pep_list, $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-wt"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$wt_inh = $ARGV[$argv_i2];
				}
			}
		}
	}
	
	@lig_pep_pdb_list = ();
	@lig_name_list = ();
	for ($lig_i = 0; $lig_i < @lig_list; $lig_i++){
		@lig_name_temp = split /\./, $lig_list[$lig_i];
		push @lig_pep_pdb_list, $lig_name_temp[0];
		push @lig_name_list, $lig_name_temp[0];
	}
	@pep_name_list = ();
	for ($pep_i = 0; $pep_i < @pep_list; $pep_i++){
		@pep_name_temp = split /\./, $pep_list[$pep_i];
		push @lig_pep_pdb_list, $pep_name_temp[0];
		push @pep_name_list, $pep_name_temp[0];
	}


	sub near_int{
		$num1 = shift;
		if ($num1 >= 0){
			$num2 = int($num1);
			$num3 = $num1 - $num2;
			if ($num3 >= 0.5){
				$near_int_num = $num2+1;
			}
			else{
				$near_int_num = $num2;
			}
		}
		else{
			$num2 = int($num1);
			$num3 = $num2 - $num1;
			if ($num3 >= 0.5){
				$near_int_num = $num2-1;
			}
			else{
				$near_int_num = $num2;
			}
		}
		return $near_int_num;
	}


	#$dfr = 0;
	#@ss_a = ();
	#@ss_b = ();
	#@chain_ID_list = ();
	#open OUTPUT, ">ssbond_list_in_bilayer_system.txt" or die "can not create!\n";
	#for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
	#	$input_file = "./chain_$cid_list[$chain_i]/ssbond_filter_model_$cid_list[$chain_i].txt";
	#	open INPUT, "$input_file" or die "can not open!\n";
	#	while(chomp($line=<INPUT>)){
	#		@items = split /\s+/, $line;
	#		if ($items[4]){
	#			$index_A = $items[4] + $dfr;
	#			$index_B = $items[7] + $dfr;
	#			print OUTPUT "$cid_list[$chain_i]\t$index_A\t$index_B\n";
	#			print "$cid_list[$chain_i]\t$index_A\t$index_B\n";
	#			push @ss_a, $index_A;
	#			push @ss_b, $index_B;
	#			push @chain_ID_list, $cid_list[$chain_i];
	#		}
	#	}
	#	close INPUT;
        #
	#	$input_file = "./chain_$cid_list[$chain_i]/chain_$cid_list[$chain_i]_itv.txt";
	#	open INPUT1, "$input_file" or die "can not open!\n";
	#	while(chomp($line=<INPUT1>)){
	#		@items = split /\s+/, $line;
	#		$seq_len = $items[1];
	#	}
	#	$dfr += $seq_len;
	#	close INPUT1;
	#}
        
	open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not open!\n";
	$domain_c = 0;
	$domain_mark = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$seg_mark = undef;
		for ($gezis_i = 71; $gezis_i <= 74 ; $gezis_i++){
			$seg_mark .= $gezis[$gezis_i];
		}
		if ($seg_mark eq "MEMB"){
			last;
		}
		if (($gezis[12] eq "H")  ||  ($gezis[13] eq "H") ) {
			next;
		}
		$atom_name_clear = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			if ($gezis[$gezis_i] ne " "){
				$atom_name_clear .= $gezis[$gezis_i];
			}
		}
		if ($atom_name_clear eq "OXT"){
			next;
		}
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
		        $atom_mark .= "$gezis[$gezis_i]";
		}
		if ((($atom_mark eq "TER   ") ||  ($atom_mark eq "TER"))  && ($domain_mark == 1))  {
			print OUTPUT "TER\n";
			close OUTPUT;
			$domain_mark = 0;
			next;
		}
		if ((($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM")) && ($domain_mark == 0))  {
			$domain_c++;
			open OUTPUT, ">bilayer_domain_$domain_c.pdb" or die "can not create!\n";
			$domain_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if ((($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM")) && ($domain_mark == 1))  {
			print OUTPUT "$line\n";
			next;
		}
	}
	close INPUT;

	@prep_lines = ();
	open INPUT, "complex_prep_hs.pdb" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		if (($gezis[12] eq "H")  || ($gezis[13] eq "H") ) {
			next;
		}
		$atom_name_clear = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			if ($gezis[$gezis_i] ne " "){
				$atom_name_clear .= $gezis[$gezis_i];
			}
		}
		if ($atom_name_clear eq "OXT"){
			next;
		}
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
		        $atom_mark .= "$gezis[$gezis_i]";
		}
		if (($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM"))  {
			push @prep_lines, $line;
		}
	}
	close INPUT;


	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		$output_file = "./prep_$cid_list[$chain_i].pdb";
		open OUTPUT, ">$output_file" or die "can not open!\n";
		for ($line_i = 0; $line_i < @prep_lines; $line_i++){
                        @gezis = split //, $prep_lines[$line_i];
			if ($gezis[21] eq $cid_list[$chain_i]){
				print OUTPUT "$prep_lines[$line_i]\n";
			}
		}
		print OUTPUT "TER\n";
		close OUTPUT;
	}

	open OUTPUT, ">15align_prep_prot.pml" or die "can not create!\n";
	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		print OUTPUT "set retain_order, 1\n";
		print OUTPUT "load prep_$cid_list[$chain_i].pdb\n";
		$domain_i = $chain_i + 1;
		print OUTPUT "load bilayer_domain_$domain_i.pdb\n";
		print OUTPUT "align prep_$cid_list[$chain_i], bilayer_domain_$domain_i\n";
		print OUTPUT "save aligned_prep_$cid_list[$chain_i].pdb,prep_$cid_list[$chain_i]\n";
	}
	print OUTPUT "quit\n";
	close OUTPUT;
	system("pymol -qc 15align_prep_prot.pml");
	$start_domain_lig = $domain_i + 1;


	########################################################################################prepare CONECT for each domain
	$C_N_dist_max = 1.6;  #this is 1.3 measured by pymol. 
	open OUTPUT1, ">14ssinfo_match.txt" or die "can not create!\n";

	$domain_count = 0;
	$itv_a = 0;
	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		$domain_count++;

		open INPUT, "aligned_prep_$cid_list[$chain_i].pdb" or die "can not open!\n";
		open OUTPUT, ">aligned_prep_ss_$cid_list[$chain_i].pdb" or die "can not open!\n";
		@atomre_lines = ();
		while(chomp($line=<INPUT>)){
			push @atomre_lines, $line;
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
			        $atom_mark .= "$gezis[$gezis_i]";
			}
			if (($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM"))  {
				print OUTPUT "$line\n";	
			}
		}
		close INPUT;
		print OUTPUT "TER\n";

		#$line_ini = 0;
		@resi_num_list_2 = (); 
		@atom_num_list_2 = (); 
		@atom_name_list_2 = ();
		@chain_ID_list_2 = (); 
		@pos_x_list_2 = ();
		@pos_y_list_2 = ();
		@pos_z_list_2 = ();
		$crt_resi_num = -99999;
		$resi_count = 0;
		#for ($line_i = $line_ini; $line_i < @atomre_lines; $line_i++){
		for ($line_i = 0; $line_i < @atomre_lines; $line_i++){
			
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
			#if ($atom_mark eq "TER   "){
				#$domain_count++;
				#$line_ini = $line_i+1;
				#last;
			#}
		}
		print OUTPUT1 "domain $domain_count\nresi_count: $resi_count\nlabel: $chain_$cid_list[$chain_i]\n";
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
		open INPUT, "./chain_$cid_list[$chain_i]/chain_$cid_list[$chain_i]_mm.aln" or die "can not open 2!\n";
		while(chomp($line=<INPUT>)){
			@re_item = split /\s+/, $line;
			if ($re_item[0] eq "chain_$cid_list[$chain_i]_em"){
				$em_seq.=$re_item[1];
			}
			if ($re_item[0] eq "chain_$cid_list[$chain_i]_md"){
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
        
		#open INPUT, "./chain_$cid_list[$chain_i]/chain_$cid_list[$chain_i]_itv.txt" or die "can not open 5!\n";
		#while(chomp($line=<INPUT>)){
		#	@itv_input = split /\s+/, $line;
		#	$itv_a = $itv_input[0];
		#}
		#close INPUT;

		push @conect_list, $fix_region[0] + $itv_a;
		push @conect_label_list, $cid_list[$chain_i];
		for ($fix_i = 1; $fix_i < @fix_region; $fix_i++){
			if ($fix_region[$fix_i] > $fix_region[$fix_i-1] + 1){
				push @conect_list, $fix_region[$fix_i-1] + $itv_a;
				push @conect_label_list, $cid_list[$chain_i];
				push @conect_list, $fix_region[$fix_i] + $itv_a;
				push @conect_label_list, $cid_list[$chain_i];
			}
		}
		push @conect_list, $fix_region[-1] + $itv_a;
		push @conect_label_list, $cid_list[$chain_i];
		print OUTPUT1 "fix region itv : @conect_list\n";
		
		$input_cnc_file = "./chain_$cid_list[$chain_i]/chain_$cid_list[$chain_i]"."_resi_CN_cnc.txt";
		open INPUT, "$input_cnc_file" or die "can not open 4!\n" ;
		$line_c = 0;
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 1){
				@cncs = split /\s+/, $line;
				push @conect_list_2, $cncs[0] + $itv_a;
				push @conect_label_list_2, $cid_list[$chain_i];
				push @conect_list_2, $cncs[1] + $itv_a;
				push @conect_label_list_2, $cid_list[$chain_i];
			}
		}
		close INPUT;
		print OUTPUT1 "CN cnc : @conect_list_2\n";


		open INPUT, "./chain_$cid_list[$chain_i]/ssbond_filter_model_$cid_list[$chain_i].txt" or die "can not open!\n";
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

		close OUTPUT;
	}
	close OUTPUT1;
	########################################################################################prepare CONECT for each domain
	
	open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not open!\n";
	open OUTPUT, ">nonprotein_forbuild.pdb" or die "can not create!\n";
	$protein_mark = 1;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$seg_mark = undef;
		for ($gezis_i = 71; $gezis_i <= 74 ; $gezis_i++){
			$seg_mark .= $gezis[$gezis_i];
		}
		if ($seg_mark eq "MEMB"){
			$protein_mark = 0;
		}
		if ($protein_mark == 0){
			print OUTPUT "$line\n";
		}
	}
	close OUTPUT;
	close INPUT;

	system("cp nonprotein_forbuild.pdb nonprotein_forbuild.pdb_bak");
	open INPUT, "nonprotein_forbuild.pdb_bak" or die "can not open!\n";
	open OUTPUT, ">nonprotein_forbuild.pdb" or die "can not create!\n";
	$crt_resi_name = "YYY";
	$crt_resi_num = 0;
	$crt_chain_id = "XXX";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
                $atom_mark = undef;
                for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
                	$atom_mark .= "$gezis[$gezis_i]";
		}
		$seg_mark = undef;
		for ($gezis_i = 71; $gezis_i <= 74 ; $gezis_i++){
			$seg_mark .= $gezis[$gezis_i];
		}
		if ($seg_mark eq "MEMB"){
			$resi_name = undef;
			for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
				$resi_name .= $gezis[$gezis_i];
			}
			$part_A = undef;
			for ($gezis_i = 0; $gezis_i <= 21 ; $gezis_i++){
				$part_A .= $gezis[$gezis_i];
			}
			$part_B = undef;
			for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
				$part_B .= $gezis[$gezis_i];
			}
			if ($crt_chain_id ne $gezis[21]){
				$crt_resi_num=1;
				$crt_resi_name = $resi_name;
				$crt_chain_id = $gezis[21];
				printf OUTPUT "$part_A%4g$part_B\n", $crt_resi_num;
				next;
			}
			elsif ($resi_name ne $crt_resi_name){
				$crt_resi_name = $resi_name;
				$crt_resi_num++;
				printf OUTPUT "$part_A%4g$part_B\n", $crt_resi_num;
				next;
			}
			else{
				printf OUTPUT "$part_A%4g$part_B\n", $crt_resi_num;
				next;
			}
		}
		elsif (($atom_mark eq "TER   ") ||  ($atom_mark eq "TER") || ($line=~"^TER.*")) {
			print OUTPUT "TER   \n";
		}
		else{
			print OUTPUT "$line\n";
		}
	}
	close OUTPUT;
	close INPUT;



	#@lig_pep_pdb_list = ();
	#open INPUT, "./4lig_pep_list.txt" or die "can not open!\n";
	#while(chomp($line=<INPUT>)){
	#	push @lig_pep_pdb_list, $line;
	#}
	#close INPUT;


	open OUTPUT1, ">15align_lig_pep.pml" or die "can not create!\n";
	for ($id_i = 0; $id_i < @lig_pep_pdb_list; $id_i++){
		print OUTPUT1 "set retain_order, 1\n";
		print OUTPUT1 "load bilayer_domain_$start_domain_lig.pdb\n";
		print OUTPUT1 "load $lig_pep_pdb_list[$id_i].pdb\n";
		print OUTPUT1 "align $lig_pep_pdb_list[$id_i],  bilayer_domain_$start_domain_lig\n";
		print OUTPUT1 "save aligned_$lig_pep_pdb_list[$id_i].pdb, $lig_pep_pdb_list[$id_i]\n";
		$start_domain_lig++;
	}
	print OUTPUT1 "quit\n";
	close OUTPUT1;

	system("pymol -qc 15align_lig_pep.pml");


	for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){

		open INPUT, "$lig_pep_pdb_list[$lig_i].pdb" or die "can not open!\n";
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



		open INPUT, "aligned_$lig_pep_pdb_list[$lig_i].pdb" or die "can not open!\n";
		open OUTPUT, ">aligned_ss_$lig_pep_pdb_list[$lig_i].pdb" or die "can not open!\n";
		@resi_num_list = ();
		@resi_name_list = ();
		@atom_num_list = ();
		@atom_name_list = ();
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
				print OUTPUT "$line\n";
			}
                }
		close INPUT;
		print OUTPUT "TER\n";

		@cnb_a = ();
		@cnb_b = ();
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

		for ($cnc_i = 0; $cnc_i < @cnb_b; $cnc_i++){
			printf OUTPUT "CONECT%5d%5d\n", $cnb_a[$cnc_i],  $cnb_b[$cnc_i];
		}

		close OUTPUT;

	}


	#@ssbond_list_w = @ssbond_list;
	#open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not open!\n";
	#open OUTPUT, ">bilayer_complex_prep_hs_edit.pdb" or die "can not create!\n";
	#while(chomp($line=<INPUT>)){
	#	@gezis = split //, $line;
	#	$atom_mark = undef;
	#	for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
	#		$atom_mark .= "$gezis[$gezis_i]";
	#	}
	#	$resi_name = undef;
	#	for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
	#		$resi_name .= "$gezis[$gezis_i]";
	#	}
	#	$resi_num = undef;
	#	for ($gezis_i = 22; $gezis_i <= 25; $gezis_i++){
	#		$resi_num .= "$gezis[$gezis_i]";
	#	}
	#	$atom_name = undef;
	#	for ($gezis_i = 12; $gezis_i <= 15; $gezis_i++){
	#		$atom_name .= "$gezis[$gezis_i]";
	#	}
	#	if ($atom_mark eq "ATOM  "){
	#		if (($resi_name eq "CYS") && ($resi_num == $ssbond_list_w[0])){
	#			$part_A = undef;
	#			for ($gezis_i = 0; $gezis_i <= 16; $gezis_i++){
	#				$part_A .= "$gezis[$gezis_i]";
	#			}
	#			$part_B = undef;
	#			for ($gezis_i = 20; $gezis_i < @gezis; $gezis_i++){
	#				$part_B .= "$gezis[$gezis_i]";
	#			}
	#			printf OUTPUT "$part_A%3s$part_B\n", "CYX";
	#			shift @ssbond_list_w;
	#		}
	#		elsif (($resi_name eq "HIE") && ($atom_name eq " HD1" )){
	#			next;
	#		}
	#		elsif (($resi_name eq "ASP") && ($atom_name eq " HD2" )){
	#			next;
	#		}
	#		elsif (($resi_name eq "GLU") && ($atom_name eq " HE2" )){
	#			next;
	#		}
	#		elsif (($resi_name eq "GLH") && ($atom_name eq " HCA" )){
	#			next;
	#		}
	#		elsif ($atom_name eq " HXT" ){
	#			next;
	#		}
	#		else{
	#			print OUTPUT "$line\n";
	#		}
	#	
	#	}
	#	else{
	#		print OUTPUT "$line\n";
	#	}
	#}
	#close INPUT;
        #
	#if (@ssbond_list_w){
	#	print "ssbond match error! ignore\n";
		#die;
		#}

	
	open OUTPUT, ">box_size.txt" or die "can not create!\n";
	open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not create!\n";
	$line_c = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
			$atom_mark .= "$gezis[$gezis_i]";
		}
		if ($atom_mark eq "ATOM  "){
			$line_c++;
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
			if ($line_c == 1){
				$up_x = $pos_x;
				$down_x = $pos_x;
				$up_y = $pos_y;
				$down_y = $pos_y;
				$up_z = $pos_z;
				$down_z = $pos_z;
			}
			else{
				if ($up_x < $pos_x){
					$up_x = $pos_x;
				}
				if ($down_x > $pos_x){
					$down_x = $pos_x;
				}
				if ($up_y < $pos_y){
					$up_y = $pos_y;
				}
				if ($down_y > $pos_y){
					$down_y = $pos_y;
				}
				if ($up_z < $pos_z){
					$up_z = $pos_z;
				}
				if ($down_z > $pos_z){
					$down_z = $pos_z;
				}
			}
		}
	
	}
	close INPUT;
	$len_x = $up_x - $down_x;
	$len_y = $up_y - $down_y;
	$len_z = $up_z - $down_z;

	$cen_x = ($up_x + $down_x)/2;
	$cen_y = ($up_y + $down_y)/2;
	$cen_z = ($up_z + $down_z)/2;

	print OUTPUT "bilayer_complex_prep_hs\t$len_x($up_x - $down_x)\t$len_y($up_y - $down_y)\t$len_z($up_z - $down_z) center: $cen_x\t$cen_y\t$cen_z\n" ; 
	close OUTPUT;


	@nature_aa = ("ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX", "HIE", "HID", "HIP", "GLH", "ASH" );

	@nst_lig = ();
	for ($lig_i = 0; $lig_i < @lig_name_list; $lig_i++){
		open INPUT, "aligned_ss_$lig_name_list[$lig_i].pdb" or die "can not open!\n";
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
				$atom_mark .= "$gezis[$gezis_i]";
			}
			if ($atom_mark eq "ATOM  "){
				$resi_name = undef;
				for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
					$resi_name .= "$gezis[$gezis_i]";
				}
				$match_naa = 0;
				for ($aa_i = 0; $aa_i < @nature_aa; $aa_i++){
					if ($nature_aa[$aa_i] eq $resi_name){
						$match_naa = 1;
						last;
					}
				}
				if ($match_naa == 0){
					$match_nst = 0;
					for ($aa_i = 0; $aa_i < @nst_lig; $aa_i++){
						if ($nst_lig[$aa_i] eq $resi_name){
							$match_nst = 1;
							last;
						}
					}
					if ($match_nst == 0){
						push @nst_lig, $resi_name;
					}
					
				}
			}
		}
		close INPUT;
	}

	@nst_pep = ();
	for ($pep_i = 0; $pep_i < @pep_name_list; $pep_i++){
		open INPUT, "apepned_ss_$pep_name_list[$pep_i].pdb" or die "can not open!\n";
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
				$atom_mark .= "$gezis[$gezis_i]";
			}
			if ($atom_mark eq "ATOM  "){
				$resi_name = undef;
				for ($gezis_i = 17; $gezis_i <= 19; $gezis_i++){
					$resi_name .= "$gezis[$gezis_i]";
				}
				$match_naa = 0;
				for ($aa_i = 0; $aa_i < @nature_aa; $aa_i++){
					if ($nature_aa[$aa_i] eq $resi_name){
						$match_naa = 1;
						last;
					}
				}
				if ($match_naa == 0){
					$match_nst = 0;
					for ($aa_i = 0; $aa_i < @nst_pep; $aa_i++){
						if ($nst_pep[$aa_i] eq $resi_name){
							$match_nst = 1;
							last;
						}
					}
					if ($match_nst == 0){
						push @nst_pep, $resi_name;
					}
					
				}
			}
		}
		close INPUT;
	}


	open OUTPUT, ">build_leap_temp.in" or die "can not create!\n";
	print OUTPUT "source leaprc.gaff
source leaprc.protein.ff19SB
source leaprc.water.tip3p
source leaprc.lipid21
loadamberparams frcmod.ionsjc_tip3p\n";
	for ($aa_i = 0; $aa_i < @nst_lig; $aa_i++){
		print OUTPUT "loadamberprep ../databases/nonaa/$nst_lig[$aa_i].prepin\n";
		print OUTPUT "loadamberparams ../databases/nonaa/$nst_lig[$aa_i]_gaff.frcmod\n";
		print OUTPUT "loadoff ../databases/nonaa/$nst_lig[$aa_i].lib\n";
	}
	for ($aa_i = 0; $aa_i < @nst_pep; $aa_i++){
		print OUTPUT "loadamberprep ../databases/nonaa/$nst_pep[$aa_i].prepin\n";
		print OUTPUT "loadamberparams ../databases/nonaa/$nst_pep[$aa_i]_ff14SB.frcmod\n";
		print OUTPUT "loadoff ../databases/nonaa/$nst_pep[$aa_i].lib\n";
	}
	print OUTPUT "npn = loadpdb nonprotein_forbuild.pdb\n";
	
	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		print OUTPUT "PR$chain_i = loadpdb aligned_prep_ss_$cid_list[$chain_i].pdb\n";
	}
	for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){
		$lig_name = "LG$lig_i";
		print OUTPUT "LG$lig_i = loadpdb aligned_ss_$lig_pep_pdb_list[$lig_i].pdb\n" ;
	}

	print OUTPUT "com = combine {";
	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		print OUTPUT "PR$chain_i ";
	}
	for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){
		print OUTPUT "LG$lig_i ";
	}
	print OUTPUT "npn}\n";

	#for ($ss_i = 0; $ss_i < @ssbond_list; $ss_i+=2){
	#	print OUTPUT "bond com.$ssbond_list[$ss_i].SG  com.$ssbond_list[$ss_i+1].SG\n";
	#}

	print OUTPUT "charge com
quit\n";

	close OUTPUT;
	system("tleap -s -f build_leap_temp.in > build_leap_temp.log");
	
	open INPUT, "build_leap_temp.log" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		#if ($line=~qw/^Total perturbed charge.*/){
		if ($line=~"^Total perturbed charge.*"){
			@items_charge = split /\s+/, $line;
			$charge_n = $items_charge[3];
		}
	}
	close INPUT;


	open OUTPUT, ">build_leap.in" or die "can not create!\n";
	print OUTPUT "source leaprc.gaff
source leaprc.protein.ff19SB
source leaprc.water.tip3p
source leaprc.lipid21
loadamberparams frcmod.ionsjc_tip3p\n";
	for ($aa_i = 0; $aa_i < @nst_lig; $aa_i++){
		print OUTPUT "loadamberprep ../databases/nonaa/$nst_lig[$aa_i].prepin\n";
		print OUTPUT "loadamberparams ../databases/nonaa/$nst_lig[$aa_i]_gaff.frcmod\n";
		print OUTPUT "loadoff ../databases/nonaa/$nst_lig[$aa_i].lib\n";
	}
	for ($aa_i = 0; $aa_i < @nst_pep; $aa_i++){
		print OUTPUT "loadamberprep ../databases/nonaa/$nst_pep[$aa_i].prepin\n";
		print OUTPUT "loadamberparams ../databases/nonaa/$nst_pep[$aa_i]_ff14SB.frcmod\n";
		print OUTPUT "loadoff ../databases/nonaa/$nst_pep[$aa_i].lib\n";
	}
	print OUTPUT "npn = loadpdb nonprotein_forbuild.pdb\n";

	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		print OUTPUT "PR$chain_i = loadpdb aligned_prep_ss_$cid_list[$chain_i].pdb\n";
	}
	for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){
		$lig_name = "LG$lig_i";
		print OUTPUT "LG$lig_i = loadpdb aligned_ss_$lig_pep_pdb_list[$lig_i].pdb\n" ;
	}
	if ($wt_inh){
		print OUTPUT "wtn = loadpdb wats_inhole_del_2.pdb\n";
		print OUTPUT "com = combine {";
		for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
			print OUTPUT "PR$chain_i ";
		}
		for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){
			print OUTPUT "LG$lig_i ";
		}
		print OUTPUT "npn wtn}\n";
	}
	else{
		print OUTPUT "com = combine {";
		for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
			print OUTPUT "PR$chain_i ";
		}
		for ($lig_i = 0; $lig_i < @lig_pep_pdb_list; $lig_i++){
			print OUTPUT "LG$lig_i ";
		}
		print OUTPUT "npn}\n";
	}



	print OUTPUT "set com box {$len_x $len_y $len_z }\n";
	if (&near_int($charge_n) < 0){
		print OUTPUT "addIonsRand com Na+ 0\n";
	}
	if (&near_int($charge_n) > 0){
		print OUTPUT "addIonsRand com Cl- 0\n";
	}

	print OUTPUT "saveAmberParm com system_start.prmtop system_start.inpcrd
savepdb com system_start.pdb
quit\n";

	close OUTPUT;
	system("tleap -s -f build_leap.in");
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

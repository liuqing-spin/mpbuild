	#input info: origin pdb, chain ID , ref pdb ID
	
	$complex_pdb = $ARGV[0];
	@crt_cID = ();
	push @crt_cID, $ARGV[1];

	open INPUT, "chain_$crt_cID[0]_raw_vs_em_resi_num.txt" or die "can not open 4!\n" ;
	$line_c=0;
	@em_resi_num_list = ();
	@raw_resi_num_list = ();
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c>1){
			@em_raw_vlist = split /\s+/, $line;
			push @em_resi_num_list, $em_raw_vlist[1];
			push @raw_resi_num_list, $em_raw_vlist[2];
		}
	}
	close INPUT;


 %aa_name = (
        "ARG" => "R",
        "HIS" => "H",
        "LYS" => "K",
        "ASP" => "D",
        "GLU" => "E",
        "SER" => "S",
        "THR" => "T",
        "ASN" => "N",
        "GLN" => "Q",
        "CYS" => "C",
        "GLY" => "G",
        "PRO" => "P",
        "ALA" => "A",
        "VAL" => "V",
        "ILE" => "I",
        "LEU" => "L",
        "MET" => "M",
        "PHE" => "F",
        "TYR" => "Y",
        "TRP" => "W",
);
	@template_pdb_list = ();
	@template_ID_list = ();
	open INPUT, "6rank_template_output_forss.txt" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;
		push @template_pdb_list, $items[0];
		push @template_ID_list, $items[1];
	}
	close INPUT;

	@tmpdb_forss = ();
	@tmpdb_forss_name = ();
	for ($pdb_i = 0; $pdb_i < @template_pdb_list; $pdb_i++){
		open INPUT, "$template_pdb_list[$pdb_i].pdb" or die "can not open!\n";
		$out_pdb = "$template_pdb_list[$pdb_i]_$template_ID_list[$pdb_i].pdb";
		$out_pdb_name = "$template_pdb_list[$pdb_i]_$template_ID_list[$pdb_i]";
		push @tmpdb_forss, $out_pdb;
		push @tmpdb_forss_name, $out_pdb_name;
		open OUTPUT, ">$out_pdb" or die "can not create!\n";
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			if ($line=~"^SSBOND.*"){
				if (($gezis[15] eq $template_ID_list[$pdb_i]) && ($gezis[29] eq $template_ID_list[$pdb_i]))  {
					print OUTPUT "$line\n";
				}
			}
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$atom_mark.= $gezis[$gezis_i];
			}
			if (($atom_mark eq "ATOM  ") && ($gezis[21] eq $template_ID_list[$pdb_i])){
				print OUTPUT "$line\n";
			}
		}
		close INPUT;
		close OUTPUT;
	}
	
	
	for ($pdb_i = 0; $pdb_i < @tmpdb_forss; $pdb_i++){
		system("../pdb2fasta $tmpdb_forss[$pdb_i] > $tmpdb_forss_name[$pdb_i].fasta");
		open INPUT, "$tmpdb_forss_name[$pdb_i].fasta" or die "can not open!\n";
		open OUTPUT, ">forss_$tmpdb_forss_name[$pdb_i].fasta" or die "can not open!\n";
		$line_c = 0;
		print OUTPUT ">$tmpdb_forss_name[$pdb_i]\n";
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 1){
				print OUTPUT "$line\n";
			}
		}
		close INPUT;
		close OUTPUT;
        
		system("cat forss_$tmpdb_forss_name[$pdb_i].fasta chain_$crt_cID[0]_em_2.fasta > forssw_$tmpdb_forss_name[$pdb_i].fasta");
		system("../clustalw2 forssw_$tmpdb_forss_name[$pdb_i].fasta");
	
	}


	#@origin_ssbond = ();
	open INPUT, "$complex_pdb" or die "can not open!\n";
	open OUTPUT, ">ssbond_chain_$crt_cID[0].txt" or die "can not create!\n";
	print OUTPUT "SSBOND XXX CYS $crt_cID[0] XXXX    CYS $crt_cID[0] XXXX\n";
	while(chomp($line=<INPUT>)){
		if ($line=~"^SSBOND.*"){
			@gezis = split //, $line;
			$cys_a = undef;
			for ($gezis_i = 17; $gezis_i <= 20 ; $gezis_i++){
				$cys_a .= $gezis[$gezis_i];
			}
			$cys_b = undef;
			for ($gezis_i = 31; $gezis_i <= 34 ; $gezis_i++){
				$cys_b .= $gezis[$gezis_i];
			}
			if (($gezis[15] eq $crt_cID[0]) && ($gezis[29] eq $crt_cID[0]))  {
				for ($num_i = 0; $num_i<@em_resi_num_list; $num_i++){
					if ($cys_a == $em_resi_num_list[$num_i]){
						$cys_a_raw = $raw_resi_num_list[$num_i];
					}
					if ($cys_b == $em_resi_num_list[$num_i]){
						$cys_b_raw = $raw_resi_num_list[$num_i];
					}
				}
				printf OUTPUT "SSBOND XXX CYS $crt_cID[0] %4g    CYS $crt_cID[0] %4g\n", $cys_a_raw, $cys_b_raw;
			}
		}
	}
	close INPUT;



	@target_seq_ss = ();
	@target_seq_ss_num = ();
	$target_seq_ss_index = -999999999;
	$target_seq_ss_resn = "XXX";
	open INPUT, "chain_$crt_cID[0].pdb" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$chain_ID = $gezis[21];
		if (($atom_mark eq "ATOM  ") && ($chain_ID eq $crt_cID[0])){
			$resi_name = undef;
			for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
				$resi_name .= $gezis[$gezis_i];
			}
			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			if (($resi_num != $target_seq_ss_index) || ($resi_name ne $target_seq_ss_resn)){
				$target_seq_ss_index = $resi_num;
				$target_seq_ss_resn = $resi_name;
				push @target_seq_ss, $aa_name{$resi_name};
				push @target_seq_ss_num, $resi_num;

			}
		}
	}
	close INPUT;

	for ($pdb_i = 0; $pdb_i < @tmpdb_forss; $pdb_i++){
		$ref_seq = undef;
		$target_seq = undef;
		open INPUT, "forssw_$tmpdb_forss_name[$pdb_i].aln" or die "can not open!\n";
		while(chomp($line=<INPUT>)){
			@re_item = split /\s+/, $line;
			if ($re_item[0] eq "$tmpdb_forss_name[$pdb_i]"){
				$ref_seq.=$re_item[1];
				@gezis = split //, $re_item[0];
				$ref_seq_ID = $gezis[-1];
			}
			if ($re_item[0] eq "chain_$crt_cID[0]_em"){
				$target_seq.=$re_item[1];
			}
		}
		close INPUT;

		@ref_seq_ss = ();
		@ref_seq_ss_num = ();
		$ref_seq_ss_index = -999999999;
		$ref_seq_ss_resn = "XXX";
		open INPUT, "$tmpdb_forss[$pdb_i]" or die "can not open!\n";
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$atom_mark .= $gezis[$gezis_i];
			}
			$chain_ID = $gezis[21];
			if (($atom_mark eq "ATOM  ") && ($chain_ID eq $template_ID_list[$pdb_i])){
				$resi_name = undef;
				for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
					$resi_name .= $gezis[$gezis_i];
				}
				$resi_num = undef;
				for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
					$resi_num .= $gezis[$gezis_i];
				}
				if (($resi_num != $ref_seq_ss_index) || ($resi_name ne $ref_seq_ss_resn)){
					$ref_seq_ss_index = $resi_num;
					$ref_seq_ss_resn = $resi_name;
					push @ref_seq_ss, $aa_name{$resi_name};
					push @ref_seq_ss_num, $resi_num;
        
				}
			}
		}
		close INPUT;
	
		@target_seq_pap = split //, $target_seq;
		$match_mark_target = 1;
		$ss_i = 0;
		print "@target_seq_pap \n @target_seq_ss\n";
		for ($seq_i = 0; $seq_i < @target_seq_pap; $seq_i++){
			if ($target_seq_pap[$seq_i] ne "-"){
				if ($target_seq_pap[$seq_i] ne $target_seq_ss[$ss_i]){
					$match_mark_target = 0;
					die "chain_$crt_cID[0] pdb sequecne is not same as pap sequence, please check!\n";
				}
				$ss_i++;
			}
		}
		if ($match_mark_target){
			print "OK! chain_$crt_cID[0] pdb sequecne is same as pap sequence\n";
		}
        
		@ref_seq_pap = split //, $ref_seq;
		$match_mark_ref = 1;
		$ss_i = 0;
		print "@ref_seq_pap \n @ref_seq_ss\n";
		for ($seq_i = 0; $seq_i < @ref_seq_pap; $seq_i++){
			if ($ref_seq_pap[$seq_i] ne "-"){
				if ($ref_seq_pap[$seq_i] ne $ref_seq_ss[$ss_i]){
					$match_mark_ref = 0;
					die "$tmpdb_forss_name[$pdb_i] pdb sequecne is not same as pap sequence, please check!\n";
				}
				$ss_i++;
			}
		}
		if ($match_mark_ref){
			print "OK! $tmpdb_forss_name[$pdb_i] pdb sequecne is same as pap sequence\n";
		}
	
		open INPUT, "$tmpdb_forss[$pdb_i]" or die "can not open!\n";
		@pote_ss_pair = ();
		while(chomp($line=<INPUT>)){
			if ($line=~"^SSBOND.*"){
				@gezis = split //, $line;
				$cys_a = undef;
				for ($gezis_i = 17; $gezis_i <= 20 ; $gezis_i++){
					$cys_a .= $gezis[$gezis_i];
				}
				$cys_b = undef;
				for ($gezis_i = 31; $gezis_i <= 34 ; $gezis_i++){
					$cys_b .= $gezis[$gezis_i];
				}
				if (($gezis[15] eq $template_ID_list[$pdb_i]) && ($gezis[29] eq  $template_ID_list[$pdb_i]))  {
					push @pote_ss_pair, $cys_a;
					push @pote_ss_pair, $cys_b;
				}
			}
		}
		close INPUT;

		for ($pote_i = 0; $pote_i < @pote_ss_pair; $pote_i+=2){
			for ($seq_i = 0; $seq_i < @ref_seq_ss; $seq_i++){
				if ($ref_seq_ss_num[$seq_i] == $pote_ss_pair[$pote_i]){
					$ss_ref_indexA = $seq_i;
				}
				if ($ref_seq_ss_num[$seq_i] == $pote_ss_pair[$pote_i+1]){
					$ss_ref_indexB = $seq_i;
				}
			}
			$pap_ref_c = -1;
			for ($seq_i2 = 0; $seq_i2 < @ref_seq_pap; $seq_i2++){
				if ($ref_seq_pap[$seq_i2] ne "-"){
					$pap_ref_c++;
				}
				if ($pap_ref_c == $ss_ref_indexA){
					$pap_target_indexA = $seq_i2;
				}
				if ($pap_ref_c == $ss_ref_indexB){
					$pap_target_indexB = $seq_i2;
				}
			}
			$pap_target_c = -1;
			$extr_ss_mark = 0;
			for ($seq_i3 = 0; $seq_i3 < @target_seq_pap; $seq_i3++){
				if ($target_seq_pap[$seq_i3] ne "-"){
					$pap_target_c++;
				}
				if (($seq_i3 == $pap_target_indexA) && ($target_seq_pap[$seq_i3] eq "C")) {
					$ss_target_indexA = $pap_target_c;
					$extr_ss_mark++;
				}
				if (($seq_i3 == $pap_target_indexB) && ($target_seq_pap[$seq_i3] eq "C")) {
					$ss_target_indexB = $pap_target_c;
					$extr_ss_mark++;
				}
			}
			if ($extr_ss_mark == 2){
				printf OUTPUT "SSBOND XXX CYS $crt_cID[0] %4g    CYS $crt_cID[0] %4g\n", $target_seq_ss_num[$ss_target_indexA], $target_seq_ss_num[$ss_target_indexB];
			}
		}
	}

	close OUTPUT;

	
	open INPUT, "ssbond_chain_$crt_cID[0].txt" or die "can not create!\n";
	#open OUTPUT1, ">ssbond_filter_$crt_cID[0].txt" or die "can not create!\n";
	open OUTPUT2, ">ssbond_filter_model_$crt_cID[0].txt" or die "can not create!\n";
	$line_c = 0;
	@ssbond_tmp_list_A = ();
	@ssbond_tmp_list_B = ();
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c > 1){
			@items = split /\s+/, $line;
			$ss_match = 0;
			for ($ss_i = 0; $ss_i < @ssbond_tmp_list_A; $ss_i++){
				if (($items[4] == $ssbond_tmp_list_A[$ss_i]) && ($items[7] == $ssbond_tmp_list_B[$ss_i])){
					$ss_match = 1;
				}
				elsif (($items[4] == $ssbond_tmp_list_B[$ss_i]) && ($items[7] == $ssbond_tmp_list_A[$ss_i])){
					$ss_match = 1;
				}
			}
			if ($ss_match == 0){
				push @ssbond_tmp_list_A , $items[4];
				push @ssbond_tmp_list_B , $items[7];
				#$model_indexA = $items[4] - ($target_seq_ss_num[0] - 1);
				#$model_indexB = $items[7] - ($target_seq_ss_num[0] - 1);
				printf OUTPUT2 "SSBOND XXX CYS $crt_cID[0] %4g    CYS $crt_cID[0] %4g\n", $items[4], $items[7];
				#printf OUTPUT1 "SSBOND XXX CYS $crt_cID[0] %4g    CYS $crt_cID[0] %4g\n", $items[4], $items[7];
				#printf OUTPUT2 "SSBOND XXX CYS $crt_cID[0] %4g    CYS $crt_cID[0] %4g\n", $model_indexA, $model_indexB;
			}
			#$ss_match_mark = 0;
			#for ($ssb_i = 0; $ssb_i < @origin_ssbond; $ssb_i+=2){
			#	if (($model_indexA == $origin_ssbond[$ssb_i]) && ($model_indexB == $origin_ssbond[$ssb_i+1])){
			#		$ss_match_mark = 1;
			#	}
			#	if (($model_indexB == $origin_ssbond[$ssb_i]) && ($model_indexA == $origin_ssbond[$ssb_i+1])){
			#		$ss_match_mark = 1;
			#	}
			#}
		}
	}
	close INPUT;
	close OUTPUT1;
	close OUTPUT2;


	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

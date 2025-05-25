	
	open INPUT, "$ARGV[0]_hstate.txt" or die "can not open!\n";

	$line_c = -1;
	@resi_h = ();
	@resn_h = ();
	@chan_h = ();
	@rest_h = ();
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c>=1){
			@items = split /\s+/, $line;
			push @resi_h, $items[0];
			push @resn_h, $items[1];
			push @chan_h, $items[2];
			push @rest_h, $items[3];
		}
	}
	close INPUT;

	@ss_a = ();
	@ss_b = ();
	@ss_label = ();
	for ($ag_i = 1; $ag_i < @ARGV; $ag_i++){
		open INPUT, "./chain_$ARGV[$ag_i]/ssbond_filter_model_$ARGV[$ag_i].txt" or die "can not open!\n";
		$line_c = 0;
		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c > 0){
				@ss_info = split /\s+/, $line;
				push @ss_a, $ss_info[4];
				push @ss_b, $ss_info[7];
				push @ss_label, $ss_info[3];
			}
			
		}
		close INPUT;
	}
	print "@ss_a\n@ss_b\n@ss_label\n";

	open INPUT, "$ARGV[0]_prep.pdb" or die "can not open!\n";
	#@name_split = split /\./, $ARGV[1];
	#$output_file = "$ARGV[0]"."_hs.pdb";
	$output_file = "complex_prep_hs.pdb";
	open OUTPUT, ">$output_file" or die "can not create!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0;$gezis_i <= 5; $gezis_i++){
			$atom_mark.=$gezis[$gezis_i];
		}
		
		$protein_mark = 0;
		$chain_ID = $gezis[21];
		for ($ag_i = 1; $ag_i < @ARGV; $ag_i++){
			if ($chain_ID eq "$ARGV[$ag_i]"){
				$protein_mark = 1;
				last;
			}
		}
		if ($protein_mark == 0){
			print OUTPUT "$line\n";
			next;
		}

		if ($atom_mark eq "ATOM  "){
			$resi_name = undef;
			for ($gezis_i = 17;$gezis_i <= 19; $gezis_i++){
				$resi_name.=$gezis[$gezis_i];
			}
			$resi_num = undef;
			for ($gezis_i = 22;$gezis_i <= 25; $gezis_i++){
				$resi_num.=$gezis[$gezis_i];
			}
			#$chain_ID = undef;
			#for ($gezis_i = 72; $gezis_i<=75; $gezis_i++){
			#	$chain_ID.= $gezis[$gezis_i];
			#}


			if (($resi_name eq "HIS") || ($resi_name eq "HSD") || ($resi_name eq "HSE") || ($resi_name eq "HSP")){
				$partA = undef;
				for ($gezis_i = 0;$gezis_i <= 16; $gezis_i++){
					$partA.=$gezis[$gezis_i];
				}
				$partB = undef;
				for ($gezis_i = 20;$gezis_i < @gezis; $gezis_i++){
					$partB.=$gezis[$gezis_i];
				}
				for ($hs_i = 0; $hs_i < @resi_h; $hs_i++){
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
					if (("HID" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
						$output_line = $partA."HID".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
					if (("HIE" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
						$output_line = $partA."HIE".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 2)){
					if (("HIP" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 2)){
						$output_line = $partA."HIP".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
				}
			}
			elsif (($resi_name eq "LYS") || ($resi_name eq "LYN") ){
				$partA = undef;
				for ($gezis_i = 0;$gezis_i <= 16; $gezis_i++){
					$partA.=$gezis[$gezis_i];
				}
				$partB = undef;
				for ($gezis_i = 20;$gezis_i < @gezis; $gezis_i++){
					$partB.=$gezis[$gezis_i];
				}
				for ($hs_i = 0; $hs_i < @resi_h; $hs_i++){
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
					if (("LYN" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
						print "$resi_name\t$resn_h[$hs_i]\t$chain_h[$hs_i]\t$rest_h[$hs_i]\n";
						$output_line = $partA."LYN".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
					if (("LYS" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
						$output_line = $partA."LYS".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
				}
			}
			elsif (($resi_name eq "ASP") || ($resi_name eq "ASH") ){
				$partA = undef;
				for ($gezis_i = 0;$gezis_i <= 16; $gezis_i++){
					$partA.=$gezis[$gezis_i];
				}
				$partB = undef;
				for ($gezis_i = 20;$gezis_i < @gezis; $gezis_i++){
					$partB.=$gezis[$gezis_i];
				}
				for ($hs_i = 0; $hs_i < @resi_h; $hs_i++){
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
					if (("ASP" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
						$output_line = $partA."ASP".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
					if (("ASH" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
						$output_line = $partA."ASH".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
				}
			}
			elsif (($resi_name eq "GLU") || ($resi_name eq "GLH") ){
				$partA = undef;
				for ($gezis_i = 0;$gezis_i <= 16; $gezis_i++){
					$partA.=$gezis[$gezis_i];
				}
				$partB = undef;
				for ($gezis_i = 20;$gezis_i < @gezis; $gezis_i++){
					$partB.=$gezis[$gezis_i];
				}
				for ($hs_i = 0; $hs_i < @resi_h; $hs_i++){
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
					if (("GLU" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 0)){
						$output_line = $partA."GLU".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
					#if (($resi_name eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
					if (("GLH" eq $resi_h[$hs_i]) && ($resi_num == $resn_h[$hs_i]) && ($chan_h[$hs_i] eq $chain_ID) && ($rest_h[$hs_i] == 1)){
						$output_line = $partA."GLH".$partB;
						print OUTPUT "$output_line\n";
						last;
					}
				}
			}
			elsif ($resi_name eq "CYS"){
				$ss_mark_1 = 0;
				$ss_mark_2 = 0;
				$partA = undef;
				for ($gezis_i = 0;$gezis_i <= 16; $gezis_i++){
					$partA.=$gezis[$gezis_i];
				}
				$partB = undef;
				for ($gezis_i = 20;$gezis_i < @gezis; $gezis_i++){
					$partB.=$gezis[$gezis_i];
				}
				for ($ss_i = 0; $ss_i < @ss_a; $ss_i++){
					if (($resi_num == $ss_a[$ss_i]) && ($chain_ID eq $ss_label[$ss_i])) {
						$output_line = $partA."CYX".$partB;
						print OUTPUT "$output_line\n";
						$ss_mark_1 = 1;
						last;
					}
				}
				for ($ss_i = 0; $ss_i < @ss_b; $ss_i++){
					if (($resi_num == $ss_b[$ss_i]) && ($chain_ID eq $ss_label[$ss_i])) {
						$output_line = $partA."CYX".$partB;
						print OUTPUT "$output_line\n";
						$ss_mark_2 = 1;
						last;
					}
				}
				if (($ss_mark_1 == 0) && ($ss_mark_2 == 0)){
					print OUTPUT "$line\n";
				}
			}
			else{
				print OUTPUT "$line\n";
			}
		}
		elsif (($atom_mark eq "TER   ")  || ($line=~"^TER.*"))   {
			print OUTPUT "$line\n";
		}
	}
	close INPUT;
	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

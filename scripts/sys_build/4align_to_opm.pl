	
	$com_pdb = undef;
	@tmm_list = ();
	@cid_list = ();
	@lig_list = ();
	@pep_list = ();
	for ($argv_i = 0; $argv_i < @ARGV ; $argv_i++){
		if ($ARGV[$argv_i] eq "-h"){
			print "	-p	target complex pdb file. 
	-tm	transmembrane chain ID.
	-ci	the list of chain IDs.
	-lg	the list of ligand pdbs.
	-lp	the list of peptide pdbs.
       	";
		die "	-h	print above information\n";
		}
		if ($ARGV[$argv_i] eq "-p"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$com_pdb = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-tm"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @tmm_list, $ARGV[$argv_i2];
				}
			}
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
	}


	open INPUT, "./tmm_id_list.txt" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@tm_id_list = split /\s+/, $line;
		$tm_ref_pdbs{$tm_id_list[0]}++;
		$tm_model_id{$tm_id_list[0]} = $tm_id_list[1];
	}
	close INPUT;

	$tmm_max_count = 0;
	$tmm_max_id = undef;
	$tm_ID = undef;
	foreach $tmpdb (keys %tm_ref_pdbs) {
		if ($tm_ref_pdbs{$tmpdb} > $tmm_max_count){
			$tmm_max_count = $tm_ref_pdbs{$tmmpdb};
			$tmm_max_id = $tmpdb;
			$tm_ID = $tm_model_id{$tmpdb};
		}
	}

	$opm_name = $tmm_max_id."_opm";
	print "opm_name: $opm_name\n";

	#open INPUT, "./chain_$tm_ID/6rank_template_output_forss.txt";
	#while(chomp($line=<INPUT>)){
	#	@items = split /\s+/, $line;
	#	if ($items[0] =~ "^ZT.*"){
	#		next;
	#	}
	#	else{
	#		$opm_name = $items[0]."_opm";
	#		last;
	#	}
	#}
	#close INPUT;

	@up_let = ('A'..'Z');
	@cid_temp_list = ();
	open INPUT, "$com_pdb" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
		        $atom_mark .= "$gezis[$gezis_i]";
		}
		if (($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM"))  {
			$cid_mark = 0;
			for ($cid_i = 0; $cid_i < @cid_temp_list; $cid_i++){
				if ($cid_temp_list[$cid_i] eq $gezis[21]){
					$cid_mark = 1;
				}
			}
			if ($cid_mark == 0){
				push @cid_temp_list, $gezis[21];
			}
		}
	}
	close INPUT;
	print "cid_temp_list: @cid_temp_list\n";

	@lig_id_list = ();
	if (@lig_list){
		for ($lig_i = 0; $lig_i < @lig_list; $lig_i++){
			for ($up_i = 0; $up_i < @up_let; $up_i++){
				$id_match = 0;
				for ($cid_i = 0; $cid_i < @cid_temp_list; $cid_i++){
					if ($cid_temp_list[$cid_i] eq $up_let[$up_i]){
						$id_match = 1;
						last;
					}
				}
				if ($id_match == 0){
					push @cid_temp_list, $up_let[$up_i];
					push @lig_id_list, $up_let[$up_i];
					$crt_cid = $up_let[$up_i];
					last;
				}
			}

			system("mv $lig_list[$lig_i]  bak_$lig_list[$lig_i]");
			open INPUT, "bak_$lig_list[$lig_i]" or die "can not open!\n";	
			open OUTPUT, ">$lig_list[$lig_i]" or die "can not open!\n";	
			while(chomp($line=<INPUT>)){
				@gezis = split //, $line;
				$atom_mark = undef;
				for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
					$atom_mark .= $gezis[$gezis_i];
				}
				if (($atom_mark eq "ATOM  ") || ($atom_mark eq "HETATM")){
					$part_A = undef;
					for ($gezis_i = 0; $gezis_i <= 20 ; $gezis_i++){
						$part_A .= $gezis[$gezis_i];
					}
					$part_B = undef;
					for ($gezis_i = 22; $gezis_i < @gezis ; $gezis_i++){
						$part_B .= $gezis[$gezis_i];
					}
					printf OUTPUT "$part_A%1s$part_B\n", $crt_cid;
				}
				if (($atom_mark eq "TER   ") || ($atom_mark eq "TER") || ($line=~"^TER.*")){
					print OUTPUT "$line\n";
				}
				if ($atom_mark eq "CONECT"){
					print OUTPUT "$line\n";
				}
			}
			close INPUT;
			close OUTPUT;

		}
	}

	@pep_id_list = ();
	if (@pep_list){
		for ($pep_i = 0; $pep_i < @pep_list; $pep_i++){
			for ($up_i = 0; $up_i < @up_let; $up_i++){
				$id_match = 0;
				for ($cid_i = 0; $cid_i < @cid_temp_list; $cid_i++){
					if ($cid_temp_list[$cid_i] eq $up_let[$up_i]){
						$id_match = 1;
						last;
					}
				}
				if ($id_match == 0){
					push @cid_temp_list, $up_let[$up_i];
					push @pep_id_list, $up_let[$up_i];
					$crt_cid = $up_let[$up_i];
					last;
				}
			}

			system("mv $pep_list[$pep_i]  bak_$pep_list[$pep_i]");
			open INPUT, "bak_$pep_list[$pep_i]" or die "can not open!\n";	
			open OUTPUT, ">$pep_list[$pep_i]" or die "can not open!\n";	
			while(chomp($line=<INPUT>)){
				@gezis = split //, $line;
				$atom_mark = undef;
				for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
					$atom_mark .= $gezis[$gezis_i];
				}
				if (($atom_mark eq "ATOM  ") || ($atom_mark eq "HETATM")){
					$part_A = undef;
					for ($gezis_i = 0; $gezis_i <= 20 ; $gezis_i++){
						$part_A .= $gezis[$gezis_i];
					}
					$part_B = undef;
					for ($gezis_i = 22; $gezis_i < @gezis ; $gezis_i++){
						$part_B .= $gezis[$gezis_i];
					}
					printf OUTPUT "$part_A%1s$part_B\n", $crt_cid;
				}
				if (($atom_mark eq "TER   ") || ($atom_mark eq "TER") || ($line=~"^TER.*")){
					print OUTPUT "$line\n";
				}
				if ($atom_mark eq "CONECT"){
					print OUTPUT "$line\n";
				}
			}
			close INPUT;
			close OUTPUT;

		}
	}

	print "cid_temp_list: @cid_temp_list\n";

	open OUTPUT, ">4align_to_opm.pml" or die "can not create!\n";
	open OUTPUT1, ">4lig_pep_list.txt" or die "can not create!\n";
	@orig_pdb = split /\./, $com_pdb;
	print OUTPUT "load $orig_pdb[0].pdb\n";
	print OUTPUT "load $opm_name.pdb\n";

	if (@lig_list){
		print OUTPUT "set retain_order, 1\n";
		for ($lig_i = 0; $lig_i < @lig_list; $lig_i++){
			print OUTPUT "load $lig_list[$lig_i]\n";
		}
		print OUTPUT  "create lig_cb, $orig_pdb[0]";
		for ($lig_i = 0; $lig_i < @lig_list; $lig_i++){
			@lig_name = split /\./, $lig_list[$lig_i];
			print OUTPUT " or $lig_name[0]";
		}
		print OUTPUT "\n";
		print OUTPUT "align lig_cb and chain $tm_ID, $opm_name\n";
		for ($lig_i = 0; $lig_i < @lig_list; $lig_i++){
			@lig_name = split /\./, $lig_list[$lig_i];
			print OUTPUT "save $lig_name[0]_align.pdb, lig_cb and chain $lig_id_list[$lig_i]\n";
			print OUTPUT1 "$lig_name[0]\n";
		}
	}

	if (@pep_list){
		print OUTPUT "set retain_order, 1\n";
		for ($pep_i = 0; $pep_i < @pep_list; $pep_i++){
			print OUTPUT "load $pep_list[$pep_i]\n";
		}
		print OUTPUT  "create pep_cb, $orig_pdb[0]";
		for ($pep_i = 0; $pep_i < @pep_list; $pep_i++){
			@pep_name = split /\./, $pep_list[$pep_i];
			print OUTPUT " or $pep_name[0]";
		}
		print OUTPUT "\n";
		print OUTPUT "align pep_cb and chain $tm_ID, $opm_name\n";
		for ($pep_i = 0; $pep_i < @pep_list; $pep_i++){
			@pep_name = split /\./, $pep_list[$pep_i];
			print OUTPUT "save $pep_name[0]_align.pdb, pep_cb and chain $pep_id_list[$pep_i]\n";
			print OUTPUT1 "$pep_name[0]\n";
		}
	}



	print OUTPUT "align $orig_pdb[0] and chain $tm_ID, $opm_name\n";

	for ($chain_i = 0; $chain_i < @cid_list ; $chain_i++){
		print OUTPUT "load ./chain_$cid_list[$chain_i]/chain_$cid_list[$chain_i]_fe_model_align_fix_refine.pdb\n";
		print OUTPUT "align chain_$cid_list[$chain_i]_fe_model_align_fix_refine, $orig_pdb[0] and chain $cid_list[$chain_i]\n";
		print OUTPUT "save chain_$cid_list[$chain_i]_fe_model_align.pdb,  chain_$cid_list[$chain_i]_fe_model_align_fix_refine\n";
		#print OUTPUT "save hetatm_$cid_list[$chain_i]_align.pdb,   $orig_pdb[0] and chain $cid_list[$chain_i] and hetatm and not resn HOH\n";
	}

	print OUTPUT "bg white\nhide line\nshow cartoon\nzoom\n";
	close OUTPUT;
	close OUTPUT1;

	system("pymol -qc 4align_to_opm.pml");
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

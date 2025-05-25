	#input:  resi name, mol2 file,  atom list for del
	#pdb file: HETATM replace by ATOM, residue name has to be edited to non-standard name, same charge atoms have to be edit names(4 char name contain same 3 chrs and 1 diff chr in last!)
	#mol2 file: atom names have to be same with pdb file.
	
	%atom_weights = (
		"C" => 6,
		"O" => 8,
		"N" => 7,
		"H" => 1,
		"P" => 15,
		"S" => 16,
		"F" => 9,
		"Br" => 35,
		"Cl" => 17,
		"Si" => 14,
	);


	open INPUT, "$ARGV[0].prepin" or die "can not open!\n";
	open OUTPUT, ">$ARGV[0]_cut.prepin" or die "can not open!\n";
	$resi_name_mark = 0;
	$dumm_start_mark = 0;
	$dumm_end_mark = 0;
	$loop_mark = 0;
	$improper_mark = 0;
	$atom_count = 3;
	@atom_name_forlib = ();
	@atom_type_forlib = ();
	@pos_x_forlib = ();
	@pos_y_forlib = ();
	@pos_z_forlib = ();
	@charge_forlib = ();
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;
		if ($items[0] eq $ARGV[0]){
			$resi_name_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		elsif (($items[2] eq "DUMM") && ($dumm_end_mark == 0) && ($dumm_end_mark == 0)){
			$dumm_start_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		elsif ((@items == 0) && ($dumm_end_mark == 0) && ($dumm_start_mark == 1)){
			$dumm_end_mark = 1;
			$dumm_start_mark = 0;
			print OUTPUT "$line\n";
			next;
		}
		elsif (($dumm_end_mark == 1) && ($items[0] eq "LOOP") && ($loop_mark == 0)){
			$loop_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		elsif (($loop_mark == 1) && (@items == 0)){
			$loop_mark = 0;
			print OUTPUT "$line\n";
			next;
		}
		elsif (($improper_mark == 0) && ($items[0] eq "IMPROPER")){
			$improper_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		elsif (($improper_mark == 1) && (@items == 0)){
			$improper_mark = 0;
			print OUTPUT "$line\n";
			next;
		}
		elsif ((@items > 8) && ($resi_name_mark == 1) && ($dumm_start_mark == 1) && ($dumm_end_mark == 0)){
			$del_match = 0;
			for ($del_i = 2; $del_i < @ARGV; $del_i++){
				if ($items[2] eq "$ARGV[$del_i]"){
					$del_match = 1;
					last;
				}
			}
			if ($del_match == 1){
				next;
			}
			push @atom_name_forlib, $items[2];
			push @atom_type_forlib, $items[3];
			push @pos_x_forlib, $items[8];
			push @pos_y_forlib, $items[9];
			push @pos_z_forlib, $items[10];
			push @charge_forlib, $items[11];

			$atom_count++;
			@gezis = split //, $line;
			$part_B = undef;
			for ($gezis_i = 4; $gezis_i < @gezis; $gezis_i++){
				$part_B.=$gezis[$gezis_i];
			}
			printf OUTPUT "%4d$part_B\n", $atom_count;
			next;
			
		}
		elsif ($loop_mark == 1){
			$del_match = 0;
			for ($del_i = 2; $del_i < @ARGV; $del_i++){
				if ($items[2] eq "$ARGV[$del_i]"){
					$del_match = 1;
					last;
				}
			}
			if ($del_match == 1){
				next;
			}
			print OUTPUT "$line\n";
		}
		elsif ($improper_mark == 1){
			$del_match = 0;
			for ($del_i = 2; $del_i < @ARGV; $del_i++){
				if ($items[2] eq "$ARGV[$del_i]"){
					$del_match = 1;
					last;
				}
			}
			if ($del_match == 1){
				next;
			}
			print OUTPUT "$line\n";
		}
		else{
			print OUTPUT "$line\n";
		}
	}
	close INPUT;
	close OUTPUT;

	open OUTPUT, ">$ARGV[0]_cut.lib" or die "can not create!\n";
	print OUTPUT "!!index array str\n \"$ARGV[0]\"\n";
	print OUTPUT "!entry.$ARGV[0].unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i + 1;
		@atom_type_forlib_gezi = split //, $atom_type_forlib[$atom_i];
		print OUTPUT " \"$atom_name_forlib[$atom_i]\" \"$atom_type_forlib[$atom_i]\" 0 1 131072 $atom_count $atom_weights{$atom_type_forlib_gezi[0]} $charge_forlib[$atom_i]\n";
		if ($atom_name_forlib[$atom_i] eq "N"){
			$cnct_start = $atom_count;
		}
		if ($atom_name_forlib[$atom_i] eq "C"){
			$cnct_end = $atom_count;
		}
	}
	print OUTPUT "!entry.$ARGV[0].unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i + 1;
		print OUTPUT " \"$atom_name_forlib[$atom_i]\" \"$atom_type_forlib[$atom_i]\" 0 -1 0.0\n";
	}
	print OUTPUT "!entry.$ARGV[0].unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.$ARGV[0].unit.childsequence single int
 2
!entry.$ARGV[0].unit.connect array int
 $cnct_start\n $cnct_end\n";
 	
	open INPUT, "$ARGV[1]" or die "can not mol2!\n";
	$bond_mark = 0;
	$atom_mark = 0;
	@atom_name_mol2 = ();
	@bond_A = ();
	@bond_B = ();
	push @atom_name_mol2, "XXXX";
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;
		if ($items[0] eq "@<TRIPOS>ATOM"){
			$atom_mark = 1;
			next;
		}
		if (($items[0] eq "@<TRIPOS>BOND") && ($atom_mark == 1)) {
			$bond_mark  = 1;
			$atom_mark = 0;
			next;
		}
		if (($items[0] =~ "^@<TRIPOS>.*") && ($bond_mark == 1)){
			$bond_mark = 0;
			last;
		}
		if ($atom_mark == 1){
			push @atom_name_mol2, $items[2];
		}
		if ($bond_mark == 1){
			push @bond_A, $items[2];
			push @bond_B, $items[3];
		}

	}
	close INPUT;

	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		for ($atom_i2 = 1; $atom_i2 < @atom_name_mol2; $atom_i2++){
			if ($atom_name_forlib[$atom_i] eq $atom_name_mol2[$atom_i2]){
				$mol2prep{$atom_i2}=$atom_i+1;
				last;
			}
		}
	}

	@bond_A_2 = ();
	@bond_B_2 = ();
	for ($bond_i = 0;$bond_i < @bond_A; $bond_i++){
		if (!$mol2prep{$bond_A[$bond_i]}){
			next;
		}
		if (!$mol2prep{$bond_B[$bond_i]}){
			next;
		}
		push @bond_A_2, $mol2prep{$bond_A[$bond_i]};
		push @bond_B_2, $mol2prep{$bond_B[$bond_i]};
	}

	@bond_A_3 = @bond_A_2;
	@bond_B_3 = @bond_B_2;
	@bond_A_4 = ();
	@bond_B_4 = ();
	@del_index = ();
	for ($cc_i = 1; $cc_i <= @bond_A_2; $cc_i++){
		$min_atom_num = 999999;
		for ($bond_i=0;$bond_i<@bond_A_3;$bond_i++){
			$del_match = 0;
			for ($del_i = 0; $del_i < @del_index; $del_i++){
				if ($del_index[$del_i] == $bond_i){
					$del_match = 1;
				}
			}
			if (($bond_A_3[$bond_i]<=$min_atom_num) && ($del_match == 0)){
				$min_index=$bond_i;
				$min_atom_num = $bond_A_3[$bond_i];
			}
		}
		push @bond_A_4, $min_atom_num;
		push @bond_B_4, $bond_B_3[$min_index];
		push @del_index, $min_index;
	}

	$crt_atom_num = $bond_A_4[0];
	@itv_list = ();
	$itv_list[0]=0;
	for ($cc_i = 1; $cc_i < @bond_A_4; $cc_i++){
		if ($bond_A_4[$cc_i] > $crt_atom_num){
			push @itv_list, $cc_i-1;
			push @itv_list, $cc_i;
			$crt_atom_num = $bond_A_4[$cc_i];
		}
	}
	push @itv_list, $cc_i-1;

	@bond_A_5 = ();
	@bond_B_5 = ();
	for ($itv_i = 0; $itv_i < @itv_list; $itv_i+=2){
	
		@del_index = ();
		for ($cc_i = $itv_list[$itv_i]; $cc_i <= $itv_list[$itv_i+1]; $cc_i++){
			$min_atom_num = 999999;
			for ($bond_i=$itv_list[$itv_i];$bond_i<=$itv_list[$itv_i+1];$bond_i++){
				$del_match = 0;
				for ($del_i = 0; $del_i < @del_index; $del_i++){
					if ($del_index[$del_i] == $bond_i){
						$del_match = 1;
					}
				}
				if (($bond_B_4[$bond_i]<=$min_atom_num) && ($del_match == 0)){
					$min_index=$bond_i;
					$min_atom_num = $bond_B_4[$bond_i];
				}
			}
			push @bond_A_5, $bond_A_4[$min_index];
			push @bond_B_5, $min_atom_num;
			push @del_index, $min_index;
		}

	}

	print OUTPUT "!entry.$ARGV[0].unit.connectivity table  int atom1x  int atom2x  int flags\n";
	for ($bond_i = 0; $bond_i < @bond_A_5; $bond_i++){
		print OUTPUT " $bond_A_5[$bond_i] $bond_B_5[$bond_i] 1\n";
	}

	print OUTPUT "!entry.$ARGV[0].unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx\n \"U\" 0 \"R\" 1\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i+1;
		print OUTPUT " \"R\" 1 \"A\" $atom_count\n";
	}

	print OUTPUT "!entry.$ARGV[0].unit.name single str
 \"$ARGV[0]\"\n";
 	print OUTPUT "!entry.$ARGV[0].unit.positions table  dbl x  dbl y  dbl z\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		print OUTPUT " $pos_x_forlib[$atom_i] $pos_y_forlib[$atom_i] $pos_z_forlib[$atom_i]\n";
	}
	print OUTPUT "!entry.$ARGV[0].unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x\n";
	print OUTPUT " $cnct_start $cnct_end 0 0 0 0\n";
	print OUTPUT "!entry.$ARGV[0].unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx\n";
	$atom_count_more = @atom_name_forlib + 1;
	print OUTPUT " \"$ARGV[0]\" 1 $atom_count_more 1 \"p\" 0\n";
	print OUTPUT "!entry.$ARGV[0].unit.residuesPdbSequenceNumber array int
 0
!entry.$ARGV[0].unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.$ARGV[0].unit.velocities table  dbl x  dbl y  dbl z\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		print OUTPUT " 0.0 0.0 0.0\n";
	}
	close OUTPUT;

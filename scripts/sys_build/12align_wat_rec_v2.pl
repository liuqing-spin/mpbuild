	

	if ($ARGV[-1] == 0){
		die "not add waters into the pocket!";
	}
	elsif ($ARGV[-1] == 1){
		print "generate the waters into pocket!";
	}
	else {
		system("cp $ARGV[-1] wats_inhole_del_2.pdb ");
		die "water pdb is provided!\n";
	}

	@tm_ids = ();
	for ($ag_i = 0; $ag_i < @ARGV-1; $ag_i++){
		$tm_ids[0] .= $ARGV[$ag_i];
	}

	open INPUT, "bilayer_complex_prep_hs.pdb" or die "can not open!\n";
	open OUTPUT, ">bilayer_complex_prep_hs_dry.pdb" or die "can not create!\n";
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;

		if ($items[-1] eq "MEMB"){
			last;
		}
		else{
			print OUTPUT "$line\n";
		}
		
	}
	close INPUT;
	close OUTPUT;



	open OUTPUT, ">align_wat_rec.pml" or die "can not create!\n";
	print OUTPUT "load chain_$tm_ids[0]_for_rism_dry.pdb\n";
	print OUTPUT "load bilayer_complex_prep_hs_dry.pdb\n";
	print OUTPUT "load wats_inhole_del.pdb\n";
	print OUTPUT "create chain_$tm_ids[0]_wats, chain_$tm_ids[0]_for_rism_dry or wats_inhole_del\n";
	print OUTPUT "align chain_$tm_ids[0]_wats, bilayer_complex_prep_hs_dry\n";
	print OUTPUT "save wats_inhole_del_align.pdb, chain_$tm_ids[0]_wats and resn WAT\nquit\n";
	close OUTPUT;

	system("pymol -qc align_wat_rec.pml");
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

	
	system("cp $ARGV[0] $ARGV[0]_bak");
	open INPUT, "$ARGV[0]" or die "can not open!\n";
	@prep_lines = ();
	while(chomp($line=<INPUT>)){
		push @prep_lines, $line;
	}
	close INPUT;

	open OUTPUT, ">$ARGV[0]" or die "can not create!\n";
	for ($line_i=0; $line_i < @prep_lines; $line_i++){
		print OUTPUT "$prep_lines[$line_i]\n";
	}


	@lig_pep_list = ();
	open INPUT, "./4lig_pep_list.txt" or die "can not open!\n";
	while(chomp($line=<INPUT>)){
		push @lig_pep_list, $line;
	}
	close INPUT;

	$atom_c = @prep_lines;
	for ($lig_i = 0; $lig_i < @lig_pep_list; $lig_i++){
		open INPUT, "$lig_pep_list[$lig_i]_align.pdb" or die "can not open!\n";
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
			}
		}
		print OUTPUT "TER\n";
		close INPUT;
	}

	close OUTPUT;

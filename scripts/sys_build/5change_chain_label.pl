
	$output_file_2 = "complex_align_label.pdb";
	open OUTPUT2, ">$output_file_2" or die "can not open!\n" ;
	
	for ($chain_i = 0; $chain_i < @ARGV; $chain_i++){
		$input_file = "chain_$ARGV[$chain_i]_fe_model_align.pdb";
		$output_file = "chain_$ARGV[$chain_i]_fe_model_align_label.pdb";
		open INPUT, "$input_file" or die "can not open!\n" ;
		open OUTPUT, ">$output_file" or die "can not open!\n"; 
		while(chomp($line=<INPUT>)){
			@gezis = split //, $line;
			$atom_mark = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$atom_mark .= $gezis[$gezis_i];
			}
			if ($atom_mark eq "ATOM  "){
				$part_A = undef;
				for ($gezis_i = 0; $gezis_i <= 20 ; $gezis_i++){
					$part_A .= $gezis[$gezis_i];
				}
				$part_B = undef;
				for ($gezis_i = 22; $gezis_i < @gezis ; $gezis_i++){
					$part_B .= $gezis[$gezis_i];
				}
				$output_line = $part_A . $ARGV[$chain_i] . $part_B;
				print OUTPUT "$output_line\n";
				print OUTPUT2 "$output_line\n";
				
			}
		}
		print OUTPUT "TER\n";
		print OUTPUT2 "TER\n";
		close OUTPUT;
		close INPUT;
	}
	close OUTPUT2;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

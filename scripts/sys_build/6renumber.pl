
	$output_file_2 = "complex_align_label_renumber.pdb";
	open OUTPUT2, ">$output_file_2" or die "can not open!\n" ;
	

	$current_num = -99999999999999;
	for ($chain_i = 0; $chain_i < @ARGV; $chain_i++){
		$input_itv = "./chain_$ARGV[$chain_i]/chain_$ARGV[$chain_i]_itv.txt";
		open INPUT1, "$input_itv" or die "can not open!\n" ;
		while(chomp($line=<INPUT1>)){
			@items = split /\s+/, $line;
			$itv_cr = $items[0];
		}
		close INPUT1;

		$input_file = "chain_$ARGV[$chain_i]_fe_model_align_label.pdb";
		$output_file = "chain_$ARGV[$chain_i]_fe_model_align_label_renumber.pdb";
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
				for ($gezis_i = 0; $gezis_i <= 21 ; $gezis_i++){
					$part_A .= $gezis[$gezis_i];
				}
				$part_B = undef;
				for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
					$part_B .= $gezis[$gezis_i];
				}
				$resi_num = undef;
				for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
					$resi_num .= $gezis[$gezis_i];
				}
				if ($resi_num != $current_num){
					$current_num = $resi_num;
					$true_num = $resi_num + $itv_cr;
					printf OUTPUT "$part_A%4d$part_B\n", $true_num;
					printf OUTPUT2 "$part_A%4d$part_B\n", $true_num;
				}
				else{
					printf OUTPUT "$part_A%4d$part_B\n", $true_num;
					printf OUTPUT2 "$part_A%4d$part_B\n", $true_num;
					
				}
				
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

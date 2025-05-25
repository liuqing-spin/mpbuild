	open INPUT, "$ARGV[0]" or "can not open $ARGV[0]!\n";
	$line_c = 0;
	@full_gezis = ();
	while (chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c>1){
			@gezis = split //, $line;
			push @full_gezis, @gezis;
		}
		else{
			@gezis = split //, $line;
			$mol_name = undef;
			for ($gezis_i = 1; $gezis_i < @gezis; $gezis_i++){
				$mol_name .= "$gezis[$gezis_i]";
			}
		}
	}
	close INPUT;

	@input_file_name_list = split /\./, $ARGV[0];
	$out_file_name = $input_file_name_list[0];
	open OUTPUT, ">$out_file_name.ali" or die "can not create!\n";
	print OUTPUT ">P1;$mol_name\nsequence:$mol_name:\::::::0.00: 0.00\n";
	$gezis_c = 0;
	$end_mark = 0;
	for ($gezis_i = 0; $gezis_i < @full_gezis; $gezis_i++){
		$gezis_c++;
		if ($gezis_c == 75){
			if ($gezis_i == @full_gezis - 1){
				print OUTPUT "$full_gezis[$gezis_i]*\n";
				$end_mark = 1;
			}
			else{
				print OUTPUT "$full_gezis[$gezis_i]\n";
				$gezis_c = 0;
			}
		}
		else{
			print OUTPUT "$full_gezis[$gezis_i]";
		}
	}
	if ($end_mark == 0){
		print OUTPUT "*\n";
	}
	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

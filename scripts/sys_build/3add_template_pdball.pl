      open INPUT, "$ARGV[0]" or die "can not open!\n";
      @chain_ID_list = ();
      $file_count = 1;
      $output_file = "ZT0"."$file_count";
      open OUTPUT, ">$output_file.pdb" or die "can not create!\n";
      while(chomp($line=<INPUT>)){
	      @items = split /\s+/, $line;
	      @gezis = split //, $line;
	      if ($items[0] eq "ATOM") {
			print OUTPUT "$line\n";
			$chain_ID = $gezis[21];
		}
		if ($items[0] eq "TER"){
			print OUTPUT "$line\n";
			push @chain_ID_list, $chain_ID;
			close OUTPUT;
			system("./pdb2fasta $output_file.pdb > $output_file.fasta");
                        $file_count++ ;
                        $output_file = "ZT0"."$file_count";
      			open OUTPUT, ">$output_file.pdb" or die "can not create!\n";
		}
      }
      close INPUT;
      close OUTPUT;
	system("./pdb2fasta $output_file.pdb > $output_file.fasta");

	$command_add_pdball = "cat";
	for ($file_i = 1; $file_i < $file_count; $file_i++){
      		$input_file = "ZT0"."$file_i";
      		open INPUT, "$input_file.fasta" or die "can not open!\n";
      		open OUTPUT, ">for_$input_file.txt" or die "can not open!\n";
		print OUTPUT "C; Produced by MODELLER

>P1;$input_file$chain_ID_list[$file_i-1]
structureX:$input_file\:26:F:506:F:MOL_ID  1; MOLECULE  PREFUSION RSV F (DS-CAV1),ENVELOPE GLYCOPROTEIN; CHAIN  F; ENGINEERED  YES:MOL_ID  1; ORGANISM_SCIENTIFIC  HUMAN RESPIRATORY SYNCYTIAL VIRUS, HUM IMMUNODEFICIENCY VIRUS 1; ORGANISM_TAXID  11250, 11676; EXPRESSION_SYSTEM  HOMO SAPIENS; EXPRESSION_SYSTEM_TAXID  9606: 2.60: 0.19\n";
		$line_c=0;
      		while(chomp($line=<INPUT>)){
			$line_c++;
			if ($line_c>1){
				$resi_count = 0;
	      			@gezis = split //, $line;
				$output_line = undef;
				for ($gezis_i = 0; $gezis_i < @gezis; $gezis_i++){
					$resi_count++;
					$output_line.=$gezis[$gezis_i];
					if ($resi_count==75){
						print OUTPUT "$output_line\n";
						$resi_count=0;
						$output_line=undef;
					}
				}
				if (($resi_count>0) && ($resi_count<75)){
					print OUTPUT "$output_line";
				}
				print OUTPUT "\*\n";
			}
		}
		close INPUT;
	
		$command_add_pdball.=" for_$input_file.txt"
	}
	
	$command_add_pdball.=" >> pdball.pir";
	#$command_add_pdball.=" >> temp.pir";
	system("$command_add_pdball");


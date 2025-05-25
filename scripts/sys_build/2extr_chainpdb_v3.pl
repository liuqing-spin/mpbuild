 @chain_array = @ARGV;
 if (@chain_array == 1){
	 open INPUT, "chain_id_list.txt" or die "can not open chain_id_list.txt!\n";
	 while(chomp($line=<INPUT>)){
	 	@chain_list = split /\s+/, $line;
	 }
	 push @chain_array, @chain_list;
	 #print "chain id not find, so read from chain_id_list.txt\n";
 }

	system("mv $ARGV[0]  $ARGV[0]_bak");
	open INPUT, "$ARGV[0]_bak" or die "can not open!\n";
	open OUTPUT, ">$ARGV[0]" or die "can not open!\n";
       	while(chomp($line=<INPUT>)){
                @gezis = split //, $line;
                $atom_mark = undef;
                for ($gezis_i = 0; $gezis_i <= 5; $gezis_i++){
                	$atom_mark .= "$gezis[$gezis_i]";
		}
		if (($atom_mark eq "ATOM  ") ||  ($atom_mark eq "HETATM")) {
			$part_A = undef;
			for ($gezis_i = 0; $gezis_i <= 65 ; $gezis_i++){
				$part_A .= $gezis[$gezis_i];
			}
			$part_B = undef;
			for ($gezis_i = 77; $gezis_i < @gezis ; $gezis_i++){
				$part_B .= $gezis[$gezis_i];
			}
			print OUTPUT "$part_A           $part_B\n";
		}
		elsif (($atom_mark eq "TER   ") ||  ($atom_mark eq "TER")) {
			print OUTPUT "TER\n";
		}
		else{
			print OUTPUT "$line\n";
		}
	}
	close INPUT;
	close OUTPUT;

	
#@chain_array = ("A","B","C","D","E","F");
 for ($chain_c = 1; $chain_c < @chain_array; $chain_c++){
      open INPUT, "$ARGV[0]" or die "can not open!\n";
      open OUTPUT, ">chain_$chain_array[$chain_c].pdb" or die "can not create!\n";
      while(chomp($line=<INPUT>)){
	      @items = split /\s+/, $line;
	      @gezis = split //, $line;
	      if (($items[0] eq "ATOM") && ($gezis[21] eq $chain_array[$chain_c])){
			print OUTPUT "$line\n";
		}
      }
      system("mkdir -p chain_$chain_array[$chain_c]");
      system("cp chain_$chain_array[$chain_c].pdb ./chain_$chain_array[$chain_c]");
      close INPUT;
      close OUTPUT;
 }
 #print  "#################################################\n";
 #print  "#                  By Liu Qing                  #\n";
 #print  "# University of Science and Technology of China #\n";
 #print  "#################################################\n";
 	shift @chain_array;
	print "@chain_array\n";

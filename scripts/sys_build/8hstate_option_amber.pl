	
	while($ARGV[0]){
		system("sleep 10s");
		open INPUT, "$ARGV[0].log" or die "can not open!\n" ;
		while(chomp($line=<INPUT>)){
			$log_line = $line;
		}
		close INPUT;
		if ($log_line=~/^DONE.*/){
			print "schrodinger prep is completed!\n";
			last;
		}

	}

	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "#   Macau University of Science and Technology  #\n";
	print  "#################################################\n";
	
	open INPUT, "$ARGV[0]_prep.pdb" or die "can not open!\n";
	#@name_split = split /\./, $ARGV[0];
	$output_file = "$ARGV[0]"."_hstate.txt";
	open OUTPUT, ">$output_file" or die "can not create!\n";


	print OUTPUT "Residue  SeqIDd  ChainID  GromacsOption\n";
	@hstate_resi = ("HIS","HSD","HSE","HSP","HID","HIE","HIP","LYS","GLU","ASP");
	$current_resi_num = 0;
	$mark0 = 0;
	$mark1 = 0;
	@list_resi_name = ();
	@list_atom_name = ();
	@list_resi_num = ();
	@list_chain_ID = ();

	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i<=19; $gezis_i++){
			$resi_name.= $gezis[$gezis_i];
		}
		$atom_name = undef;
		for ($gezis_i = 12; $gezis_i<=15; $gezis_i++){
			$atom_name.= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i<=25; $gezis_i++){
			$resi_num.= $gezis[$gezis_i];
		}
		#$chain_ID = undef;
		#for ($gezis_i = 72; $gezis_i<=75; $gezis_i++){
		#	$chain_ID.= $gezis[$gezis_i];
		#}
		$chain_ID = $gezis[21];

		for ($hresi_i = 0; $hresi_i < @hstate_resi; $hresi_i++){
			if ($resi_name eq $hstate_resi[$hresi_i]){
				push @list_resi_name, $resi_name;
				push @list_atom_name, $atom_name;
				push @list_resi_num, $resi_num;
				push @list_chain_ID, $chain_ID;
			}
		}
	}
	close INPUT;

	@list_HIS = ();
	@list_LYS = ();
	@list_GLU = ();
	@list_ASP = ();
	$temp_resi_num = $list_resi_num[0];
	$temp_chain_ID = $list_chain_ID[0];
	for ($list_i = 0; $list_i < @list_resi_name; $list_i++){
		if (($list_resi_num[$list_i] > $temp_resi_num) || ($list_chain_ID[$list_i] ne $temp_chain_ID) || ($list_i == @list_resi_name -1) ){
			if (($list_resi_name[$list_i-1] eq "HIS") || ($list_resi_name[$list_i-1] eq "HID") || ($list_resi_name[$list_i-1] eq "HIE") || ($list_resi_name[$list_i-1] eq "HIP")){
				push @list_HIS, "$list_resi_num[$list_i-1]"."_$list_chain_ID[$list_i-1]";
			}
			if ($list_resi_name[$list_i-1] eq "LYS"){
				push @list_LYS, "$list_resi_num[$list_i-1]"."_$list_chain_ID[$list_i-1]";
			}
			if ($list_resi_name[$list_i-1] eq "GLU"){
				push @list_GLU, "$list_resi_num[$list_i-1]"."_$list_chain_ID[$list_i-1]";
			}
			if ($list_resi_name[$list_i-1] eq "ASP"){
				push @list_ASP, "$list_resi_num[$list_i-1]"."_$list_chain_ID[$list_i-1]";
			}
			$temp_resi_num = $list_resi_num[$list_i];
			$temp_chain_ID = $list_chain_ID[$list_i];
		}
	}

	for ($list_i = 0; $list_i < @list_HIS; $list_i++){
		@item = split /_/, $list_HIS[$list_i];
		for ($list_i2 = 0; $list_i2 < @list_resi_name; $list_i2++){
			if ((($list_resi_name[$list_i2] eq "HIS" ) || ($list_resi_name[$list_i2] eq "HID" ) || ($list_resi_name[$list_i2] eq "HIE" )  || ($list_resi_name[$list_i2] eq "HIP" )) && ($list_resi_num[$list_i2] == $item[0]) && ($list_chain_ID[$list_i2] eq $item[1])){
				if ($list_atom_name[$list_i2] eq " HD1"){
					$mark0 = 1;
				}
				elsif ($list_atom_name[$list_i2] eq " HE2"){
					$mark1 = 1;
				}
				
			}
		}
		if (($mark0 == 1) && ($mark1 == 0)) {
			print OUTPUT "HID @item 0\n";
		}
		elsif (($mark0 == 0) && ($mark1 == 1)) {
			print OUTPUT "HIE @item 1\n";
		}
		elsif (($mark0 == 1) && ($mark1 == 1)) {
			print OUTPUT "HIP @item 2\n";
		}
		$mark0 = 0;
		$mark1 = 0;
	}

	for ($list_i = 0; $list_i < @list_LYS; $list_i++){
		@item = split /_/, $list_LYS[$list_i];
		for ($list_i2 = 0; $list_i2 < @list_resi_name; $list_i2++){
			if (($list_resi_name[$list_i2] eq "LYS" ) && ($list_resi_num[$list_i2] == $item[0]) && ($list_chain_ID[$list_i2] eq $item[1])){
				if ($list_atom_name[$list_i2] eq " HZ1"){
					$mark1 = 1;
				}
			}
		}
		if($mark1 == 1){
			print OUTPUT "LYS @item 1\n";
		}
		else {
			print OUTPUT "LYN @item 0\n";
		}
		$mark1 = 0;
	}
        
	for ($list_i = 0; $list_i < @list_ASP; $list_i++){
		@item = split /_/, $list_ASP[$list_i];
		for ($list_i2 = 0; $list_i2 < @list_resi_name; $list_i2++){
			if (($list_resi_name[$list_i2] eq "ASP" ) && ($list_resi_num[$list_i2] == $item[0]) && ($list_chain_ID[$list_i2] eq $item[1])){
				if ($list_atom_name[$list_i2] eq " HD2"){
					$mark1 = 1;
				}
			}
		}
		if($mark1 == 1){
			print OUTPUT "ASH @item 1\n";
		}
		else {
			print OUTPUT "ASP @item 0\n";
		}
		$mark1 = 0;
	}
        
	for ($list_i = 0; $list_i < @list_GLU; $list_i++){
		@item = split /_/, $list_GLU[$list_i];
		for ($list_i2 = 0; $list_i2 < @list_resi_name; $list_i2++){
			if (($list_resi_name[$list_i2] eq "GLU" ) && ($list_resi_num[$list_i2] == $item[0]) && ($list_chain_ID[$list_i2] eq $item[1])){
				if ($list_atom_name[$list_i2] eq " HE2"){
					$mark1 = 1;
				}
			}
		}
		if($mark1 == 1){
			print OUTPUT "GLH @item 1\n";
		}
		else {
			print OUTPUT "GLU @item 0\n";
		}
		$mark1 = 0;
	}

	close OUTPUT;
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

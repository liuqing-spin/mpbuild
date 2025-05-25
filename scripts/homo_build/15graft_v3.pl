
	@items = split /_/, $ARGV[0];
	open OUTPUT, ">15align_model_to_raw.pml" or die "can not open 3!\n";
	print OUTPUT "load $ARGV[0].pdb\n";
	print OUTPUT "load $ARGV[0]_fe_model.pdb\n";
	print OUTPUT "align $ARGV[0]_fe_model, $ARGV[0]\n";
	print OUTPUT "save $ARGV[0]_fe_model_align.pdb, $ARGV[0]_fe_model\n";
	close OUTPUT;
	system("pymol -qc 15align_model_to_raw.pml");

	$em_seq = undef;
	$md_seq = undef;
	open INPUT, "chain_$items[1]_mm.aln" or die "can not open 2!\n";
	while(chomp($line=<INPUT>)){
		@re_item = split /\s+/, $line;
		if ($re_item[0] eq "chain_$items[1]_em"){
			$em_seq.=$re_item[1];
		}
		if ($re_item[0] eq "chain_$items[1]_md"){
			$md_seq.=$re_item[1];
		}
	}

	@md_seq_list = split //, $md_seq;
	@em_seq_list = split //, $em_seq;
	@fix_region = ();
	for ($em_i = 0; $em_i < @em_seq_list; $em_i++){
		if ($em_seq_list[$em_i] eq "-") {
			push @fix_region, $em_i + 1;
		}
	}
	close INPUT; 
	print "fix_region: @fix_region\n";
	open OUTPUT, ">chain_$items[1]_fix_region.txt" or die "can not open 3!\n";
	print OUTPUT "fix_region: @fix_region\n";
	close OUTPUT;

	open INPUT, "chain_$items[1]_itv.txt" or die "can not open 5!\n";
	while(chomp($line=<INPUT>)){
		@itv_input = split /\s+/, $line;
		$itv_a = $itv_input[0];
	}
	close INPUT;

	@ss_a = ();
	@ss_b = ();
	open INPUT, "ssbond_filter_model_$items[1].txt" or die "can not open 6!\n" ;
	while(chomp($line=<INPUT>)){
		@ss_items = split /\s+/, $line;
		push @ss_a, $ss_items[4];
		push @ss_b, $ss_items[7];
	}
	close INPUT;

	open INPUT, "$ARGV[0]_fe_model_align.pdb" or die "can not open 4!\n" ;
	@fix_lines = ();
	@fix_resi_lines = ();
	@ss_lines_A = ();
	@ss_lines_B = ();
	@ss_resi_lines = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		if (($gezis[12] eq "H") || ($gezis[13] eq "H")){
			next;
		}
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
			$resi_num .= $gezis[$gezis_i];
		}
		$resi_name = undef;
		for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
			$resi_name .= $gezis[$gezis_i];
		}

		$fix_match = 0;
		for ($fix_i = 0; $fix_i < @fix_region; $fix_i++){
			if ($resi_num == $fix_region[$fix_i]){
				$fix_match = 1;
			}
		}

		if (($atom_mark eq "ATOM  ") && ($fix_match == 1)) {
			push @fix_lines, $line;
			push @fix_resi_lines, $resi_num;

		}

		$ss_mark = 0;
		for ($ss_i = 0; $ss_i < @ss_a; $ss_i++){
			if ($ss_a[$ss_i] == $resi_num){
				$ss_mark = 1;
				last
			}
		}
		for ($ss_i = 0; $ss_i < @ss_b; $ss_i++){
			if ($ss_b[$ss_i] == $resi_num){
				$ss_mark = 1;
				last
			}
		}
		if (($resi_name eq "CYS") && ($ss_mark == 1)){
			$part_A = undef;
			for ($gezis_i = 0; $gezis_i <= 16 ; $gezis_i++){
				$part_A .= $gezis[$gezis_i];
			}
			$part_B = undef;
			for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
				$part_B .= $gezis[$gezis_i];
			}
			#$ss_md_line = printf $part_A."CYX $items[1]%4d".$part_B, $resi_num+$itv_a;
			push @ss_lines_A, $part_A;
			push @ss_lines_B, $part_B;
			push @ss_resi_lines, $resi_num+$itv_a;
		}
	}
	close INPUT;


	open INPUT, "$ARGV[0].pdb" or die "can not open 10!\n" ;
	open OUTPUT, ">$ARGV[0]_fe_model_align_graft.pdb" or die "can not create!\n";
	$fix_ini = 0;
	$fix_end_mark = 0;
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		if (($gezis[12] eq "H") || ($gezis[13] eq "H")){
			next;
		}
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
			$resi_num .= $gezis[$gezis_i];
		}
		$resi_num_rn = $resi_num - $itv_a;
		$part_A = undef;
		for ($gezis_i = 0; $gezis_i <= 20 ; $gezis_i++){
			$part_A .= $gezis[$gezis_i];
		}
		$part_A.="A";
		$part_B = undef;
		for ($gezis_i = 26; $gezis_i < @gezis ; $gezis_i++){
			$part_B .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			if ($fix_end_mark == 1){
				printf OUTPUT "$part_A%4d$part_B\n", $resi_num_rn;
				next;
			}
			for ($fix_i = $fix_ini; $fix_i < @fix_resi_lines; $fix_i++){
				if ($fix_resi_lines[$fix_i] > $resi_num_rn){
					printf OUTPUT "$part_A%4d$part_B\n", $resi_num_rn;
					$fix_ini = $fix_i;
					last;
				}
				elsif ($fix_resi_lines[$fix_i] < $resi_num_rn){
					print OUTPUT "$fix_lines[$fix_i]\n";

				}
				else{
					print "fix region wrong!\n";
					die;
				}
				if ($fix_i == @fix_resi_lines - 1){
					$fix_end_mark = 1;
					printf OUTPUT "$part_A%4d$part_B\n", $resi_num_rn;

				}
			}
			if (@fix_resi_lines == 0){
				printf OUTPUT "$part_A%4d$part_B\n", $resi_num_rn;
			}
		}
	}
	close INPUT;
	if (($fix_ini < @fix_resi_lines) && ($fix_end_mark == 0)) {
		for ($fix_i = $fix_ini; $fix_i < @fix_resi_lines; $fix_i++){
			print OUTPUT "$fix_lines[$fix_i]\n";
		}
	}

	close OUTPUT;


	open OUTPUT, ">15save_graft.pml" or die "can not open 3!\n";
	print OUTPUT "load $ARGV[0]_fe_model_align_graft.pdb\n";
	print OUTPUT "save $ARGV[0]_graft.pdb, $ARGV[0]_fe_model_align_graft\n";
	close OUTPUT;
	system("pymol -qc 15save_graft.pml");

	@conect_list = ();
	push @conect_list, $fix_region[0];
	for ($fix_i = 1; $fix_i < @fix_region; $fix_i++){
		if ($fix_region[$fix_i] > $fix_region[$fix_i-1] + 1){
			push @conect_list, $fix_region[$fix_i-1];
			push @conect_list, $fix_region[$fix_i];
		}
	}
	push @conect_list, $fix_region[-1];
	open OUTPUT, ">>chain_$items[1]_fix_region.txt" or die "can not open 3!\n";
	print OUTPUT "fix region itv : @conect_list\n";
	close OUTPUT;
	
	@conect_list_2 = ();
	$input_cnc_file = "chain_$items[1]"."_resi_CN_cnc.txt";
	open INPUT, "$input_cnc_file" or die "can not open 4!\n" ;
	$line_c = 0;
	while(chomp($line=<INPUT>)){
		$line_c++;
		if ($line_c > 1){
			@cncs = split /\s+/, $line;
			push @conect_list_2, $cncs[0];
			push @conect_list_2, $cncs[1];
		}
	}
	close INPUT;


	@resi_num_list_2 = ();
	@atom_num_list_2 = ();
	@atom_name_list_2 = ();
	open INPUT, "$ARGV[0]_graft.pdb" or die "can not open 10!\n" ;
	open OUTPUT, ">$ARGV[0]_graft_conect.pdb" or die "can not create!\n";
	open OUTPUT1, ">graft_cnc_match.txt" or die "can not create!\n";
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		$resi_num = undef;
		for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
			$resi_num .= $gezis[$gezis_i];
		}
		$atom_num = undef;
		for ($gezis_i = 6; $gezis_i <= 10 ; $gezis_i++){
			$atom_num .= $gezis[$gezis_i];
		}
		$atom_name = undef;
		for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
			$atom_name .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			print OUTPUT "$line\n";
			push @resi_num_list_2, $resi_num;
			push @atom_num_list_2, $atom_num;
			push @atom_name_list_2, $atom_name;
		}
	}
	print OUTPUT "TER\n";

	for ($cne_i = 0; $cne_i <= 1; $cne_i+=2){
		@conect_pair_b = ();
		$mark_3 = 0;
		$mark_4 = 0;
		for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
			if (($conect_list[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_3 = 1;
			}
			if (($conect_list[$cne_i + 1] + 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_4 = 1;
			}
		}
		if (($mark_3 == 1)  && ($mark_4 == 1) ){
			print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
			}
		else{
			print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
		}
	}
	for ($cne_i = @conect_list-2; $cne_i < @conect_list; $cne_i+=2){
		@conect_pair_a = ();
		$mark_1 = 0;
		$mark_2 = 0;
		for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
			if (($conect_list[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")){
				push @conect_pair_a, $atom_num_list_2[$ret_i];
				$mark_1 = 1;
			}
			if (($conect_list[$cne_i] - 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ")){
				push @conect_pair_a, $atom_num_list_2[$ret_i];
				$mark_2 = 1;
			}
		}
		if (($mark_1 == 1) && ($mark_2 == 1) ){
			print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			printf OUTPUT "CONECT%5d%5d\n", $conect_pair_a[0],$conect_pair_a[1];
			}
		else{
			print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
		}
	}
	for ($cne_i = 2; $cne_i < @conect_list-2; $cne_i+=2){
		@conect_pair_a = ();
		@conect_pair_b = ();
		$mark_1 = 0;
		$mark_2 = 0;
		$mark_3 = 0;
		$mark_4 = 0;
		for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
			if (($conect_list[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")){
				push @conect_pair_a, $atom_num_list_2[$ret_i];
				$mark_1 = 1;
			}
			if (($conect_list[$cne_i] - 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ")){
				push @conect_pair_a, $atom_num_list_2[$ret_i];
				$mark_2 = 1;
			}
			if (($conect_list[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_3 = 1;
			}
			if (($conect_list[$cne_i + 1] + 1 == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_4 = 1;
			}
		}
		if (($mark_1 == 1) && ($mark_2 == 1)  && ($mark_3 == 1)  && ($mark_4 == 1) ){
			print "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			printf OUTPUT "CONECT%5d%5d\n", $conect_pair_a[0],$conect_pair_a[1];
			printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
			}
		else{
			print "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
			print OUTPUT1 "NO loop terminals match: $conect_list[$cne_i]  $conect_list[$cne_i+1]!\n";
		}
	}
	for ($cne_i = 0; $cne_i < @conect_list_2; $cne_i+=2){
		@conect_pair_b = ();
		$mark_3 = 0;
		$mark_4 = 0;
		for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
			if (($conect_list_2[$cne_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " C  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_3 = 1;
			}
			if (($conect_list_2[$cne_i + 1] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " N  ")){
				push @conect_pair_b, $atom_num_list_2[$ret_i];
				$mark_4 = 1;
			}
		}
		if (($mark_3 == 1)  && ($mark_4 == 1) ){
			print "loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
			print OUTPUT1 "loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
			printf OUTPUT "CONECT%5d%5d\n", $conect_pair_b[0],$conect_pair_b[1];
			}
		else{
			print "NO loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
			print OUTPUT1 "NO loop terminals match2: $conect_list_2[$cne_i]  $conect_list_2[$cne_i+1]!\n";
		}
	}
	for ($ss_i = 0; $ss_i < @ss_a; $ss_i++){
		@ss_pair = ();
		$mark_1 = 0;
		$mark_2 = 0;
		for ($ret_i = 0; $ret_i < @resi_num_list_2; $ret_i++){
			if (($ss_a[$ss_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " SG ")){
				push @ss_pair, $atom_num_list_2[$ret_i];
				$mark_1 = 1;
			}
			if (($ss_b[$ss_i] == $resi_num_list_2[$ret_i]) && ($atom_name_list_2[$ret_i] eq " SG ")){
				push @ss_pair, $atom_num_list_2[$ret_i];
				$mark_2 = 1;
			}
		}
		if (($mark_1 == 1) && ($mark_2 == 1) ){
			print "s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
			print OUTPUT1 "s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
			printf OUTPUT "CONECT%5d%5d\n", $ss_pair[0],$ss_pair[1];
		}
		else{
			print "NO s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
			print OUTPUT1 "NO s-s bond match: $ss_a[$ss_i]  $ss_b[$ss_i]!\n";
		}
		
	}

	close INPUT;
	close OUTPUT;

	#system("\$SCHRODINGER/utilities/prepwizard -rehtreat -disulfides -propka_pH 7.0 -rmsd 5 -fillsidechains $ARGV[0]_graft_conect.pdb $ARGV[0]_graft_conect_prep.pdb");
	system("$ARGV[1]/utilities/prepwizard -rehtreat -disulfides -propka_pH 7.0 -rmsd 20   -fillsidechains $ARGV[0]_graft_conect.pdb $ARGV[0]_graft_conect_prep.pdb");

	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

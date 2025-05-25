	#input: pdb file name, resi name, nproc , gaussian exe, mol2 file anme, predefine atom group name
	#pdb file: HETATM replace by ATOM, residue name has to be edited to non-standard name, same charge atoms have to be edit names(4 char name contain same 3 chrs and 1 diff chr in last!)
	#mol2 file: atom names have to be same with pdb file.
	#predefine atom group name file were converted from prepin files, in which the columns of  atom name and charge were maintained, and the first file column has to be whitespace; file name show "groupname_charge.txt", the groupname has to be "3-letter" format.
	

	$nonaa_pdb = undef;
	$nonaa_mol2 = undef;
	@nonaa_name = ();
	$gau_name = undef;
	$thr_num = 1;
	$nonaa_charge = 0;
	@pre_group = ();
	for ($argv_i = 0; $argv_i < @ARGV ; $argv_i++){
		if ($ARGV[$argv_i] eq "-h"){
			print "		-pdb	the pdb file of nonaa with caps.
	-mol2	the mol2 file of nonaa without caps, and with the same atom names as the pdb file.
	-name	the nonaa name.
	-gau	the execute name of Gaussian.
	-np	the threads number used.
	-ac	the charge of nonaa.
	-pf	the atom group name with predefined atom charges. The charge infos have to be provided in the txt file like 'XXX'_charge.txt.
       	";
		die "	-h	print above information\n";
		}
		if ($ARGV[$argv_i] eq "-pdb"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$nonaa_pdb = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-mol2"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$nonaa_mol2 = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-name"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @nonaa_name,  $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-gau"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$gau_name = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-np"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$thr_num = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-ac"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					$nonaa_charge = $ARGV[$argv_i2];
				}
			}
		}
		if ($ARGV[$argv_i] eq "-pf"){
			for ($argv_i2 = $argv_i+1; $argv_i2 < @ARGV; $argv_i2++){
				if ($ARGV[$argv_i2] =~ "^-"){
					last;
				}
				else{
					push @pre_group,  $ARGV[$argv_i2];
				}
			}
		}
	}



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

	@gr_atom_name_list = ();
	@gr_atom_nametop_list = ();
	@gr_name_list = ();
	@gr_charge_list = ();
	for ($grpu_i = 5; $grpu_i < @pre_group; $grpu_i++){
		open INPUT, "$pre_group[$grpu_i]_charge.txt" or die "can not open!\n";
		while(chomp($line=<INPUT>)){
			@items = split /\s+/, $line;

			@gezi_temp = split //, $items[1];
			if (@gezi_temp == 1){
				push @gr_atom_name_list, " $items[1]  ";
				push @gr_atom_nametop_list, " $items[1] ";
			}
			elsif (@gezi_temp == 2){
				push @gr_atom_name_list, " $items[1] ";
				push @gr_atom_nametop_list, " $items[1]";
			}
			elsif (@gezi_temp == 3){
				push @gr_atom_name_list, " $items[1]";
				$atom_name_top = " ";
				for ($gezi_i = 0; $gezi_i < @gezi_temp-1; $gezi_i++){
					$atom_name_top.=$gezi_temp[$gezi_i];
				}
				push @gr_atom_nametop_list, $atom_name_top;
			}
			elsif (@gezi_temp == 4){
				push @gr_atom_name_list, "$items[1]";
				$atom_name_top = undef;
				for ($gezi_i = 0; $gezi_i < @gezi_temp-1; $gezi_i++){
					$atom_name_top.=$gezi_temp[$gezi_i];
				}
				push @gr_atom_nametop_list, $atom_name_top;
			}

			push @gr_name_list, $pre_group[$grpu_i];
			push @gr_charge_list, $items[2];

		}
		close INPUT;
	}


	open OUTPUT, ">$nonaa_name[0]_opt.com" or die "can not create!\n";
	print OUTPUT "%nproc=$thr_num
%mem=512MB
%chk=$nonaa_name[0]_opt.chk
#P b3lyp/6-31g* opt

$nonaa_name[0]

$nonaa_charge 1
";
	open INPUT, "$nonaa_pdb" or die "can not open!\n";
	@atom_type_list = ();
	@pos_x_list = ();
	@pos_y_list = ();
	@pos_z_list = ();
	$atom_count = 0;
	$atom_name_list = ();
	$atom_nametop_list = ();
	$resi_name_list = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		@items = split /\s+/, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if (($atom_mark eq "ATOM  ") || ($atom_mark eq "HETATM")) {
			$atom_count++;
			$atom_name = undef;
			$atom_type_mark = 0;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				$atom_name .= $gezis[$gezis_i];
				if (($gezis[$gezis_i] ne " ") && ($atom_type_mark == 0)){
					$atom_type = $gezis[$gezis_i];
					$atom_type_mark = 1;
					if ($atom_type ne "$items[-1]"){
						print "##############error? atom type wrong!\n";
					}
				}
			}
			push @atom_type_list, $items[-1];

			$atom_name = undef;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				$atom_name .= $gezis[$gezis_i];
			}
			push @atom_name_list, $atom_name;

			$atom_nametop = undef;
			for ($gezis_i = 12; $gezis_i <= 14 ; $gezis_i++){
				$atom_nametop .= $gezis[$gezis_i];
			}
			push @atom_nametop_list, $atom_nametop;

			$resi_name = undef;
			for ($gezis_i = 17; $gezis_i <= 19 ; $gezis_i++){
				$resi_name .= $gezis[$gezis_i];
			}
			push @resi_name_list, $resi_name;


			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
				$pos_x .= $gezis[$gezis_i];
			}
			push @pos_x_list, $pos_x;

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
				$pos_y .= $gezis[$gezis_i];
			}
			push @pos_y_list, $pos_y;

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
				$pos_z .= $gezis[$gezis_i];
			}
			push @pos_z_list, $pos_z;
		}
	}
	close INPUT;


	for ($atom_i = 0; $atom_i < @atom_type_list; $atom_i++){
		printf OUTPUT " $atom_type_list[$atom_i]            %12.8f  %12.8f  %12.8f\n", $pos_x_list[$atom_i], $pos_y_list[$atom_i], $pos_z_list[$atom_i];
	}
	print OUTPUT "\n";
	close OUTPUT;



	system("$gau_name $nonaa_name[0]_opt.com");


	open OUTPUT, ">$nonaa_name[0]_hf.com" or die "can not create!\n";
	print OUTPUT "%nprocshared=$thr_num
%mem=512MB
%chk=$nonaa_name[0]_hf.chk
#P HF/6-31G*  SCF=Tight Pop=MK IOp(6/33=2)

$nonaa_name[0]

$nonaa_charge 1
";
	for ($atom_i = 0; $atom_i < @atom_type_list; $atom_i++){
		printf OUTPUT " $atom_type_list[$atom_i]            %12.8f  %12.8f  %12.8f\n", $pos_x_list[$atom_i], $pos_y_list[$atom_i], $pos_z_list[$atom_i];
	}
	print OUTPUT "\n";
	close OUTPUT;
	system("$gau_name $nonaa_name[0]_hf.com");

	system("espgen -i $nonaa_name[0]_hf.log -o esp.dat");


	open OUTPUT, ">resp.in" or die "can not create!\n";
	print OUTPUT "capped-resp run
 &cntrl
 nmol=1,
 ihfree=1,
 qwt=0.0005,
 iqopt=2,
 /\n";
 	print OUTPUT "    1\n$nonaa_name[0]\n";
	printf OUTPUT " %4d %4d\n", $nonaa_charge,$atom_count;

	@charge_indi = ();
	for ($atom_i = 0; $atom_i < $atom_count; $atom_i++){
		$charge_indi[$atom_i] = 0;
	}
	for ($atom_i = 0; $atom_i < $atom_count; $atom_i++){
		#if (($atom_i < 6) || ($atom_i >= $atom_count - 6)){
		if ($atom_i < 6) {
			$charge_indi[$atom_i] = -1;
		}
		else{
			$gr_match_2 = 0;
			for ($grpu_i = 0; $grpu_i < @gr_name_list; $grpu_i++){
				if ($resi_name_list[$atom_i] eq $gr_name_list[$grpu_i]){
					$charge_indi[$atom_i] = -1;
					$gr_match_2 = 1;
					last;
				}
			}
			if ($gr_match_2 == 0){
				for ($atom_i2 = $atom_i + 1; $atom_i2 < $atom_count; $atom_i2++){
					@crt_atom_name_rry = split //, $atom_name_list[$atom_i];
					@fpr_atom_name_rry = split //, $atom_name_list[$atom_i2];
					$diff_mark = 0;
					for ($chr_i = 0; $chr_i < @crt_atom_name_rry - 1; $chr_i++){
						if ($crt_atom_name_rry[$chr_i] ne $fpr_atom_name_rry[$chr_i]){
							$diff_mark = 1;
						}
					}
					if ($diff_mark == 0){
						$charge_indi[$atom_i2] = $atom_i + 1;
						
					}
				}
			}
		}
	}

	for ($atom_i = 0; $atom_i < $atom_count; $atom_i++){
		printf OUTPUT " %4d %4d\n", $atom_weights{$atom_type_list[$atom_i]},$charge_indi[$atom_i];
	}

	print OUTPUT "\n";
	close OUTPUT;


	open OUTPUT, ">resp.qin" or die "can not create!\n";
	for ($atom_i = 0; $atom_i < $atom_count; $atom_i++){
		if ($atom_i < 6){
			if ($atom_name_list[$atom_i] eq " C  "){
				printf OUTPUT "%10.6f", 0.5972;
			}
			elsif ($atom_name_list[$atom_i] eq " O  "){
				printf OUTPUT "%10.6f", -0.5679;
			}
			elsif ($atom_name_list[$atom_i] eq " CH3"){
				printf OUTPUT "%10.6f", -0.3662;
			}
			else{
				printf OUTPUT "%10.6f", 0.1123;
			}
		}
		#elsif ( $atom_i >= $atom_count - 6){
		#	if ($atom_name_list[$atom_i] eq " N  "){
		#		printf OUTPUT "%10.6f", -0.4157;
		#	}
		#	elsif ($atom_name_list[$atom_i] eq " H  "){
		#		printf OUTPUT "%10.6f", 0.2719;
		#	}
		#	elsif ($atom_name_list[$atom_i] eq " CA "){
		#		printf OUTPUT "%10.6f", -0.1490;
		#	}
		#	else{
		#		printf OUTPUT "%10.6f", 0.0976;
		#	}
		#}
		#elsif (( $atom_i < $atom_count - 6) && ($atom_i >=6)){
		elsif ($atom_i >=6){
			$gr_match = 0;
			for ($grpu_i = 0; $grpu_i < @gr_atom_name_list;$grpu_i++){
				if (($resi_name_list[$atom_i] eq $gr_name_list[$grpu_i]) && ($atom_nametop_list[$atom_i] eq $gr_atom_nametop_list[$grpu_i])){
					printf OUTPUT "%10.6f", $gr_charge_list[$grpu_i];
					$gr_match = 1;
					last;
				}
			}
			if ($gr_match == 0) {
				printf OUTPUT "%10.6f", 0;
			}
		}
		if (($atom_i + 1) % 8 == 0 ){
			print OUTPUT "\n";
		}
	}
	if ($atom_count % 8 != 0 ){
		print OUTPUT "\n";
	}
	close OUTPUT;

	system("resp -O -i resp.in -o resp.out -p resp.pch -t resp.chg  -q resp.qin -e esp.dat");
	system("antechamber -fi gout -i $nonaa_name[0]_hf.log -bk $nonaa_name[0] -fo ac -o $nonaa_name[0].ac -c rc -cf resp.chg -at amber");

	for ($atom_i = 0; $atom_i < $atom_count; $atom_i++){
		if (($atom_name_list[$atom_i] eq " CA ") && ($resi_name_list[$atom_i] eq $nonaa_name[0]))   {
			$main_chain_indi = $atom_i + 1;
		}
		if (($atom_name_list[$atom_i] eq " N  ") && ($resi_name_list[$atom_i] eq $nonaa_name[0])) {
			$head_name_indi = $atom_i + 1;
		}
		if (($atom_name_list[$atom_i] eq " C  ") && ($resi_name_list[$atom_i] eq $nonaa_name[0])) {
			$tail_name_indi = $atom_i + 1;
		}
	}

	open INPUT, "$nonaa_name[0].ac" or die "can not open!\n";
	$atom_count_2 = 0;
	@atom_name_list_2 = ();
	@omit_name_list_head = ();
	@omit_name_list_tail = ();
	while(chomp($line=<INPUT>)){
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if ($atom_mark eq "ATOM  "){
			$atom_count_2++;
			$atom_name = undef;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				if ($gezis[$gezis_i] ne " "){
					$atom_name .= $gezis[$gezis_i];
				}
			}
			push @atom_name_list_2, $atom_name;
			if ($atom_count_2 == $main_chain_indi){
				$main_chain_name = $atom_name;
			}
			if ($atom_count_2 == $head_name_indi){
				$head_name = $atom_name;
			}
			if ($atom_count_2 == $tail_name_indi){
				$tail_name = $atom_name;
			}
			if ($atom_count_2 <= 6){
				push @omit_name_list_head, $atom_name;
			}
			#if ($atom_count_2 > $atom_count - 6){
			#	push @omit_name_list_tail, $atom_name;
			#}
		
		}
	}
	close INPUT;

	if ($atom_count != $atom_count_2){
		die "pdb and ac file not match!\n";
	}

	open OUTPUT, ">$nonaa_name[0].mc" or die "can not open!\n";
	print OUTPUT "HEAD_NAME $head_name\n";
	#print OUTPUT "TAIL_NAME $tail_name\n";
	print OUTPUT "MAIN_CHAIN $main_chain_name\n";
	for ($name_i = 0; $name_i < @omit_name_list_head; $name_i++){
		print OUTPUT "OMIT_NAME $omit_name_list_head[$name_i]\n";
	}
	#for ($name_i = 0; $name_i < @omit_name_list_tail; $name_i++){
	#	print OUTPUT "OMIT_NAME $omit_name_list_tail[$name_i]\n";
	#}
	#print OUTPUT "PRE_HEAD_TYPE C\nPOST_TAIL_TYPE N\nCHARGE $nonaa_charge\n";
	print OUTPUT "PRE_HEAD_TYPE C\nCHARGE $nonaa_charge\n";
	close OUTPUT;

	system("prepgen -i $nonaa_name[0].ac -o $nonaa_name[0].prepin -m $nonaa_name[0].mc -rn $nonaa_name[0]");
	system("parmchk2 -i $nonaa_name[0].ac -f ac -o $nonaa_name[0]_gaff.frcmod");
	system("parmchk2 -i $nonaa_name[0].ac -f ac -o $nonaa_name[0]_ff14SB.frcmod -a Y -p \$AMBERHOME/dat/leap/parm/parm10.dat");


	system("mv $nonaa_name[0].prepin $nonaa_name[0].prepin_bak");
	open INPUT, "$nonaa_name[0].prepin_bak" or die "can not open!\n";
	open OUTPUT, ">$nonaa_name[0].prepin" or die "can not open!\n";

	for ($atom_i = 0; $atom_i < @atom_name_list; $atom_i++){
		@pdb_atom_name = split //, $atom_name_list[$atom_i];
		$pdb_atom_clear = undef;
		for ($name_i = 0; $name_i < @pdb_atom_name; $name_i++){
			if ($pdb_atom_name[$name_i] ne " "){
				$pdb_atom_clear .= $pdb_atom_name[$name_i];
			}
		}
		push @pdb_atom_clear_list, $pdb_atom_clear;
	}


	$resi_name_mark = 0;
	$dumm_start_mark = 0;
	$dumm_end_mark = 0;
	$impr_start_mark = 0;
	$impr_end_mark = 0;
	$loop_start_mark = 0;
	$loop_end_mark = 0;
	@atom_type_order = ();
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;
		if ($items[0] eq $nonaa_name[0]){
			$resi_name_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if (($items[2] eq "DUMM") && ($dumm_end_mark == 0) && ($dumm_end_mark == 0)){
			$dumm_start_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items == 0) && ($dumm_end_mark == 0) && ($dumm_start_mark == 1)){
			$dumm_end_mark = 1;
			$dumm_start_mark = 0;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items > 8) && ($resi_name_mark == 1) && ($dumm_start_mark == 1) && ($dumm_end_mark == 0)){
			@gezis = split //, $line;
			$atom_name = undef;
			for ($gezis_i = 6; $gezis_i <= 9 ; $gezis_i++){
				if($gezis[$gezis_i] ne " "){
					$atom_name .= $gezis[$gezis_i];
				}
			}
			$part_A = undef;
			for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
				$part_A .= $gezis[$gezis_i];
			}
			$part_B = undef;
			for ($gezis_i = 10; $gezis_i < @gezis ; $gezis_i++){
				$part_B .= $gezis[$gezis_i];
			}
			for ($atom_i = 0; $atom_i < @atom_name_list_2; $atom_i++){
				if ($atom_name eq $atom_name_list_2[$atom_i]){
					printf OUTPUT "$part_A%-4s$part_B\n", $pdb_atom_clear_list[$atom_i];
					push @atom_type_order, $atom_type_list[$atom_i];
				}	
			}
			next;
			
		}
		if (($items[0] eq "IMPROPER") && ($resi_name_mark == 1) && ($impr_start_mark == 0) && ($impr_end_mark == 0)){
			$impr_start_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
			
		if ((@items == 0) && ($resi_name_mark == 1) && ($impr_start_mark == 1) && ($impr_end_mark == 0)){
			$impr_start_mark = 0;
			$impr_end_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items != 0) && ($resi_name_mark == 1) && ($impr_start_mark == 1) && ($impr_end_mark == 0)){
			for ($items_i = 1; $items_i < @items; $items_i++){
				@item_atom_name = split //, $items[$items_i];
				$item_atom_clear = undef;
				for ($name_i = 0; $name_i < @item_atom_name; $name_i++){
					if ($item_atom_name[$name_i] ne " "){
						$item_atom_clear .= $item_atom_name[$name_i];
					}
				}
				$item_match = 0;
				#for ($atom_i = 6; $atom_i < @atom_name_list_2 -6 ; $atom_i++){
				for ($atom_i = 6; $atom_i < @atom_name_list_2 ; $atom_i++){
					if ($item_atom_clear eq $atom_name_list_2[$atom_i]){
						printf OUTPUT " %4s",  $pdb_atom_clear_list[$atom_i];
						$item_match = 1;
						last;
					}
				}
				if ($item_match == 0){
						printf OUTPUT " %4s",  $items[$items_i];
				}
			}
			print OUTPUT "\n";
			next;
		}

		if (($items[0] eq "LOOP") && ($resi_name_mark == 1) && ($loop_start_mark == 0) && ($loop_end_mark == 0)){
			$loop_start_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
			
		if ((@items == 0) && ($resi_name_mark == 1) && ($loop_start_mark == 1) && ($loop_end_mark == 0)){
			$loop_start_mark = 0;
			$loop_end_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items != 0) && ($resi_name_mark == 1) && ($loop_start_mark == 1) && ($loop_end_mark == 0)){
			for ($items_i = 1; $items_i < @items; $items_i++){
				@item_atom_name = split //, $items[$items_i];
				$item_atom_clear = undef;
				for ($name_i = 0; $name_i < @item_atom_name; $name_i++){
					if ($item_atom_name[$name_i] ne " "){
						$item_atom_clear .= $item_atom_name[$name_i];
					}
				}
				$item_match = 0;
				#for ($atom_i = 6; $atom_i < @atom_name_list_2 -6 ; $atom_i++){
				for ($atom_i = 6; $atom_i < @atom_name_list_2  ; $atom_i++){
					if ($item_atom_clear eq $atom_name_list_2[$atom_i]){
						printf OUTPUT " %4s",  $pdb_atom_clear_list[$atom_i];
						$item_match = 1;
						last;
					}
				}
				if ($item_match == 0){
						printf OUTPUT " %4s",  $items[$items_i];
				}
			}
			print OUTPUT "\n";
			next;
		}
		else{
			print OUTPUT "$line\n";
		}
	}
	close OUTPUT;
	close INPUT;

	open INPUT, "$nonaa_name[0].prepin" or die "can not open!\n";
	$resi_name_mark = 0;
	$dumm_start_mark = 0;
	$dumm_end_mark = 0;
	@atom_name_forlib = ();
	@atom_type_forlib = ();
	@pos_x_forlib = ();
	@pos_y_forlib = ();
	@pos_z_forlib = ();
	@charge_forlib = ();
	while(chomp($line=<INPUT>)){
		@items = split /\s+/, $line;
		if ($items[0] eq $nonaa_name[0]){
			$resi_name_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if (($items[2] eq "DUMM") && ($dumm_end_mark == 0) && ($dumm_end_mark == 0)){
			$dumm_start_mark = 1;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items == 0) && ($dumm_end_mark == 0) && ($dumm_start_mark == 1)){
			$dumm_end_mark = 1;
			$dumm_start_mark = 0;
			print OUTPUT "$line\n";
			next;
		}
		if ((@items > 8) && ($resi_name_mark == 1) && ($dumm_start_mark == 1) && ($dumm_end_mark == 0)){
			#@gezis = split //, $line;
			#$atom_name = undef;
			#for ($gezis_i = 6; $gezis_i <= 9 ; $gezis_i++){
			#	if($gezis[$gezis_i] ne " "){
			#		$atom_name .= $gezis[$gezis_i];
			#	}
			#}
			push @atom_name_forlib, $items[2];
			push @atom_type_forlib, $items[3];
			push @pos_x_forlib, $items[8];
			push @pos_y_forlib, $items[9];
			push @pos_z_forlib, $items[10];
			push @charge_forlib, $items[11];
			next;
			
		}
	}
	close INPUT;

	open OUTPUT, ">$nonaa_name[0].lib" or die "can not create!\n";
	print OUTPUT "!!index array str\n \"$nonaa_name[0]\"\n";
	print OUTPUT "!entry.$nonaa_name[0].unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i + 1;
		print OUTPUT " \"$atom_name_forlib[$atom_i]\" \"$atom_type_forlib[$atom_i]\" 0 1 131072 $atom_count $atom_weights{$atom_type_order[$atom_i]} $charge_forlib[$atom_i]\n";
		if ($atom_name_forlib[$atom_i] eq "N"){
			$cnct_start = $atom_count;
		}
		#if ($atom_name_forlib[$atom_i] eq "C"){
		#	$cnct_end = $atom_count;
		#}
	}
	print OUTPUT "!entry.$nonaa_name[0].unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i + 1;
		print OUTPUT " \"$atom_name_forlib[$atom_i]\" \"$atom_type_forlib[$atom_i]\" 0 -1 0.0\n";
	}
	print OUTPUT "!entry.$nonaa_name[0].unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.$nonaa_name[0].unit.childsequence single int
 2
!entry.$nonaa_name[0].unit.connect array int
 $cnct_start\n 0\n";
 	
	open INPUT, "$nonaa_mol2" or die "can not mol2!\n";
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

	print OUTPUT "!entry.$nonaa_name[0].unit.connectivity table  int atom1x  int atom2x  int flags\n";
	for ($bond_i = 0; $bond_i < @bond_A_5; $bond_i++){
		print OUTPUT " $bond_A_5[$bond_i] $bond_B_5[$bond_i] 1\n";
	}

	print OUTPUT "!entry.$nonaa_name[0].unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx\n \"U\" 0 \"R\" 1\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		$atom_count = $atom_i+1;
		print OUTPUT " \"R\" 1 \"A\" $atom_count\n";
	}

	print OUTPUT "!entry.$nonaa_name[0].unit.name single str
 \"$nonaa_name[0]\"\n";
 	print OUTPUT "!entry.$nonaa_name[0].unit.positions table  dbl x  dbl y  dbl z\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		print OUTPUT " $pos_x_forlib[$atom_i] $pos_y_forlib[$atom_i] $pos_z_forlib[$atom_i]\n";
	}
	print OUTPUT "!entry.$nonaa_name[0].unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x\n";
	#print OUTPUT " $cnct_start $cnct_end 0 0 0 0\n";
	print OUTPUT " $cnct_start 0 0 0 0 0\n";
	print OUTPUT "!entry.$nonaa_name[0].unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx\n";
	$atom_count_more = @atom_name_forlib + 1;
	print OUTPUT " \"$nonaa_name[0]\" 1 $atom_count_more 1 \"p\" 0\n";
	print OUTPUT "!entry.$nonaa_name[0].unit.residuesPdbSequenceNumber array int
 0
!entry.$nonaa_name[0].unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.$nonaa_name[0].unit.velocities table  dbl x  dbl y  dbl z\n";
	for ($atom_i = 0; $atom_i < @atom_name_forlib; $atom_i++){
		print OUTPUT " 0.0 0.0 0.0\n";
	}
	close OUTPUT;

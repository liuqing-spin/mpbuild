	
	opendir (DIR, '.') or die "can not open dir!\n";
	open OUTPUT, ">align_select_model.pml" or die "can not create!\n";
	@file_list_cur = readdir DIR;
	@items_a = split /_/, $ARGV[1];
	print OUTPUT "load chain_$items_a[1].pdb, EM_struc\n";
	for ($file_i = 0; $file_i < @file_list_cur; $file_i++){
		if ($file_list_cur[$file_i] =~ /^$ARGV[1]\.B99/){
			@file_name = split /\./, $file_list_cur[$file_i];
			print OUTPUT "load $file_list_cur[$file_i]\n";
			print OUTPUT "align $file_name[0].$file_name[1] , EM_struc\n";
			print OUTPUT "save $file_name[0]_$file_name[1]_align.pdb,  $file_name[0].$file_name[1]\n";
			#print "$file_list_cur[$file_i]\n";
		}
	}
	closedir DIR;
	close OUTPUT;
	system("pymol -qc align_select_model.pml > align_select_model.log");

	@model_list = ();
	@rmsd_list = ();
	open INPUT, "align_select_model.log" or die "can not open!\n";
	$model_mark = 0;
	while(chomp($line=<INPUT>)){
		if ($line=~/^PyMOL>align.*/){
			#print "$line\n";
			@items = split /\s+/, $line;
			@model_name = split /\./, $items[1];
			$align_name = $model_name[0]."_$model_name[1]"."_align.pdb";
			push @model_list, $align_name;
			next;
		}
		if ($line=~/^ ExecutiveAlign.*/){
			$model_mark = 1;
			next;
		}
		if (($line=~/^ Executive: RMS.*/) && ($model_mark == 1)){
			$model_mark = 0;
			@items = split /\s+/, $line;
			print "rmsd: $items[4]\n";
			push @rmsd_list , $items[4];
		}
	}
	close INPUT;

	$em_seq = undef;
	$md_seq = undef;
	@items = split /_/, $ARGV[1];
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

	open INPUT, "$ARGV[0]" or die "can not open 3!\n" ;
	@line_list = ();
	@pos_x_list = ();
	@pos_y_list = ();
	@pos_z_list = ();
	@resi_num_list = ();
	while(chomp($line=<INPUT>)){
		push @line_list, $line;
		@gezis = split //, $line;
		$atom_mark = undef;
		for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
			$atom_mark .= $gezis[$gezis_i];
		}
		if (($atom_mark eq "ATOM  ") && ($gezis[21] ne $items[1])) {
			$atom_name = undef;
			for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
				$atom_name .= $gezis[$gezis_i];
			}

			$pos_x = undef;
			for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
				$pos_x .= $gezis[$gezis_i];
			}

			$pos_y = undef;
			for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
				$pos_y .= $gezis[$gezis_i];
			}

			$pos_z = undef;
			for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
				$pos_z .= $gezis[$gezis_i];
			}

			$resi_num = undef;
			for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
				$resi_num .= $gezis[$gezis_i];
			}
			push @pos_x_list, $pos_x;
			push @pos_y_list, $pos_y;
			push @pos_z_list, $pos_z;
			push @resi_num_list, $resi_num;
			#print "get EM struture!: $pos_x $pos_y $pos_z $resi_num\n";
		}
	}
	close INPUT;
                
	$rmsd_cutoff = 2;
	@min_atom_dist_list = ();
	for ($list_i = 0; $list_i < @rmsd_list; $list_i++){
		if ($rmsd_list[$list_i] < $rmsd_cutoff){
			open INPUT, "$model_list[$list_i]" or die "can not open 4!\n" ;
			@line_list_2 = ();
			@pos_x_list_2 = ();
			@pos_y_list_2 = ();
			@pos_z_list_2 = ();
			@resi_num_list_2 = ();
			while(chomp($line=<INPUT>)){
				push @line_list_2, $line;
				@gezis = split //, $line;
				$atom_mark = undef;
				for ($gezis_i = 0; $gezis_i <= 5 ; $gezis_i++){
					$atom_mark .= $gezis[$gezis_i];
				}
				$resi_num = undef;
				for ($gezis_i = 22; $gezis_i <= 25 ; $gezis_i++){
					$resi_num .= $gezis[$gezis_i];
				}

				$fix_match = 0;
				for ($fix_i = 0; $fix_i < @fix_region; $fix_i++){
					if ($resi_num == $fix_region[$fix_i]){
						$fix_match = 1;
						last;
					}
				}

				if (($atom_mark eq "ATOM  ") && ($fix_match == 1)) {
					$atom_name = undef;
					for ($gezis_i = 12; $gezis_i <= 15 ; $gezis_i++){
						$atom_name .= $gezis[$gezis_i];
					}
                
					$pos_x = undef;
					for ($gezis_i = 30; $gezis_i <= 37 ; $gezis_i++){
						$pos_x .= $gezis[$gezis_i];
					}
                
					$pos_y = undef;
					for ($gezis_i = 38; $gezis_i <= 45 ; $gezis_i++){
						$pos_y .= $gezis[$gezis_i];
					}
                
					$pos_z = undef;
					for ($gezis_i = 46; $gezis_i <= 53 ; $gezis_i++){
						$pos_z .= $gezis[$gezis_i];
					}
                
					push @pos_x_list_2, $pos_x;
					push @pos_y_list_2, $pos_y;
					push @pos_z_list_2, $pos_z;
					push @resi_num_list_2, $resi_num;
					#print "get model: $pos_x $pos_y $pos_z $resi_num\n";
                
				}
			}
			close INPUT;
			$min_atom_dist = 9999999;
			for ($pos_i = 0; $pos_i < @pos_x_list; $pos_i++){
				for ($pos_i_2 = 0; $pos_i_2 < @pos_x_list_2; $pos_i_2++){
						$atom_dist = sqrt(($pos_x_list[$pos_i] - $pos_x_list_2[$pos_i_2])**2 + ($pos_y_list[$pos_i] - $pos_y_list_2[$pos_i_2])**2 + ($pos_z_list[$pos_i] - $pos_z_list_2[$pos_i_2])**2 );
						#print "atom_dist:$atom_dist\n";
						if ($atom_dist < $min_atom_dist){
							$min_atom_dist = $atom_dist;
						}
				}
				
			}
			push @min_atom_dist_list , $min_atom_dist;
                
		}
	}
	
	print "min_atom_dist_list:@min_atom_dist_list\n";
	$max_model_dist = 0;
	for ($min_i = 0; $min_i < @min_atom_dist_list; $min_i++){
		if ($min_atom_dist_list[$min_i] > $max_model_dist){
			$max_model_dist = $min_atom_dist_list[$min_i] ;
			$model_index = $min_i;
		}
	}

	print "max_atom_dist: $max_model_dist in $model_list[$model_index]\n";
	system("cp $model_list[$model_index] $ARGV[1]_model.pdb");

	#open OUTPUT, ">check_spst.pml" or die "can not create pml!\n";
	#print OUTPUT "load $ARGV[1]_model.pdb\n";
	#print OUTPUT "load $ARGV[0], EM_struct\n";
	#print OUTPUT "align $ARGV[1]_model, EM_struct\n";
	#print OUTPUT "bg white\nhide line\nshow cartoon\n";
	#print OUTPUT "color red, $ARGV[1]_model\n";
	#print OUTPUT "show stick, resn CYS\nlabel resn CYS, resi\nzoom\n";
	#close OUTPUT; 

	#system ("pymol check_spst.pml");

	
	#open OUTPUT, ">manual_select_model.pml" or die "can not create pml!\n";
	#print OUTPUT "load $ARGV[0], EM_struct\n";
	#print OUTPUT "load $ARGV[1].B99990001.pdb\n";
	#print OUTPUT "load $ARGV[1].B99990002.pdb\n";
	#print OUTPUT "load $ARGV[1].B99990003.pdb\n";
	#print OUTPUT "load $ARGV[1].B99990004.pdb\n";
	#print OUTPUT "load $ARGV[1].B99990005.pdb\n";
	#print OUTPUT "align $ARGV[1].B99990001, EM_struct\n";
	#print OUTPUT "align $ARGV[1].B99990002, EM_struct\n";
	#print OUTPUT "align $ARGV[1].B99990003, EM_struct\n";
	#print OUTPUT "align $ARGV[1].B99990004, EM_struct\n";
	#print OUTPUT "align $ARGV[1].B99990005, EM_struct\n";
	#print OUTPUT "bg white\nhide line\nshow cartoon\n";
	#print OUTPUT "show stick, resn CYS\nlabel resn CYS, resi\nzoom\n";
	#close OUTPUT;

	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

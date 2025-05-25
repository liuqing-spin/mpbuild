
	


 open INPUT, "$ARGV[0]" or die "can not open!\n";
 @miss_resi = ();
 @miss_chain = ();
 @miss_num = ();
 @stru_resi = ();
 @stru_chain = ();
 @stru_num = ();
 while(chomp($line=<INPUT>)){
   @items = split /\s+/, $line;
   if (($items[0] eq "REMARK") && ($items[1] == 465 )){
        push @miss_resi, $items[2];
	push @miss_chain, $items[3];
	push @miss_num, $items[4];
   }
   if ($items[0] eq "ATOM"){
	if (!$stru_num[0]){
		@gezi = split //, $line;
		$resi_num = "$gezi[22]"."$gezi[23]"."$gezi[24]"."$gezi[25]";
		push @stru_num, $resi_num;
		$resi_name = "$gezi[17]"."$gezi[18]"."$gezi[19]";
                push @stru_resi, $resi_name;
		$chain_name = "$gezi[21]";
		push @stru_chain, $chain_name;
        }
	else{
		@gezi = split //, $line;
		$resi_num = "$gezi[22]"."$gezi[23]"."$gezi[24]"."$gezi[25]";
		if ($resi_num == $stru_num[@stru_num-1]){
	                next;
	        }
		else{
			push @stru_num, $resi_num;
			$resi_name = "$gezi[17]"."$gezi[18]"."$gezi[19]";
                        push @stru_resi, $resi_name;
			$chain_name = "$gezi[21]";
			push @stru_chain, $chain_name;
		}
	}
   }
 }

 @chain_type = ();
 push @chain_type, $stru_chain[0];
 for ($resi_c = 1; $resi_c < @stru_chain; $resi_c++){
 	$chain_match = 0;
	for ($chain_s = 0; $chain_s < @chain_type; $chain_s++){
	 	if ($stru_chain[$resi_c] eq $chain_type[$chain_s]){
			$chain_match = 1;
			last;
		}
	 }
	if ($chain_match == 0){
		push @chain_type, $stru_chain[$resi_c];
	}
 }
 print "@chain_type\n";
 close INPUT;

open OUTPUT, ">chain_id_list.txt" or die "can not create!\n";
print OUTPUT "@chain_type\n";
close OUTPUT;


 %aa_name = (
 	"ARG" => "R",
	"HIS" => "H",
	"LYS" => "K",
	"ASP" => "D",
	"GLU" => "E",
	"SER" => "S",
	"THR" => "T",
	"ASN" => "N",
	"GLN" => "Q",
	"CYS" => "C",
	"GLY" => "G",
	"PRO" => "P",
	"ALA" => "A",
	"VAL" => "V",
	"ILE" => "I",
	"LEU" => "L",
	"MET" => "M",
	"PHE" => "F",
	"TYR" => "Y",
	"TRP" => "W",
);
  
#foreach $key (keys %aa_name){
	 #print "$key $aa_name{$key}\n";
#}

 for ($chain_s = 0; $chain_s < @chain_type; $chain_s++){
	$file_out = "chain_$chain_type[$chain_s].txt";
	$file_out1 = "chain_$chain_type[$chain_s].fasta";
	@stru2_num = ();
	@stru2_resi = ();
	@stru2_chain = ();
	open OUTPUT, ">$file_out" or die "can not create!\n";
	open OUTPUT1, ">$file_out1" or die "can not create!\n";
	print OUTPUT1 ">chain$chain_type[$chain_s]\n";
 	for ($resi_c = 0; $resi_c < @stru_num; $resi_c++){
		if ($stru_chain[$resi_c] eq $chain_type[$chain_s]){
		    push @stru2_num, $stru_num[$resi_c];
		    push @stru2_resi, $stru_resi[$resi_c];
		    push @stru2_chain, $stru_chain[$resi_c];
		}
	}
 	for ($miss_c = 0; $miss_c < @miss_num; $miss_c++){
		if ($miss_chain[$miss_c] eq $chain_type[$chain_s]){
		    push @miss2_num, $miss_num[$miss_c];
		    push @miss2_resi, $miss_resi[$miss_c];
		    push @miss2_chain, $miss_chain[$miss_c];
		}
	}
	$miss2_mark = 0;
	$fasta_c = 0;
        for ($resi2_c = 0; $resi2_c < @stru2_num; $resi2_c++){
                for ($miss2_c = $miss2_mark; $miss2_c < @miss2_num; $miss2_c++){
              		if (($stru2_num[$resi2_c] > $miss2_num[$miss2_c]) && ($stru2_chain[$resi2_c] eq $miss2_chain[$miss2_c])){
              			print OUTPUT "$miss2_num[$miss2_c] $miss2_chain[$miss2_c] $miss2_resi[$miss2_c]\n";
              			$miss2_mark = $miss2_c+1;
				if ($fasta_c<80){
				    print OUTPUT1 "$aa_name{$miss2_resi[$miss2_c]}";
				    $fasta_c++;
				}
				else{
				    print OUTPUT1 "$aa_name{$miss2_resi[$miss2_c]}\n";
				    $fasta_c = 0;
				}

              		}	
             	 }
              	print OUTPUT "$stru2_num[$resi2_c] $stru2_chain[$resi2_c] $stru2_resi[$resi2_c]\n";
		if ($fasta_c<80){
		    print OUTPUT1 "$aa_name{$stru2_resi[$resi2_c]}";
		    $fasta_c++;
		}
		else{
		    print OUTPUT1 "$aa_name{$stru2_resi[$resi2_c]}\n";
		    $fasta_c = 0;
		}
             	if ($resi2_c == @stru2_num - 1) {
                     	for ($miss2_c = $miss2_mark; $miss2_c < @miss2_num; $miss2_c++){
                     	 	if (($stru2_num[$resi2_c] < $miss2_num[$miss2_c]) && ($stru2_chain[$resi2_c] eq $miss2_chain[$miss2_c])){
                             		print OUTPUT "$miss2_num[$miss2_c] $miss2_chain[$miss2_c] $miss2_resi[$miss2_c]\n";
					if ($fasta_c<80){
					    print OUTPUT1 "$aa_name{$miss2_resi[$miss2_c]}";
					    $fasta_c++;
					}
					else{
					    print OUTPUT1 "$aa_name{$miss2_resi[$miss2_c]}\n";
					    $fasta_c = 0;
					}
               		        }	
                    	}
		}
        }
       	print "TER\n";
	print OUTPUT1 "\n";
	close OUTPUT;
	close OUTPUT1;
 }
 # $miss_mark = 0;
 # for ($resi_c = 0; $resi_c < @stru_num; $resi_c++){
 #	for ($miss_c = $miss_mark; $miss_c < @miss_num; $miss_c++){
 #		if (($stru_num[$resi_c] > $miss_num[$miss_c]) && ($stru_chain[$resi_c] eq $miss_chain[$miss_c])){
 ##			print "$miss_num[$miss_c] $miss_chain[$miss_c] $miss_resi[$miss_c]\n";
 #       		$miss_mark = $miss_c+1;
 #       	}	
 #       }
 ##       if (($resi_c > 0) && ($stru_chain[$resi_c] eq $stru_chain[$resi_c-1])){
 #       	print "$stru_num[$resi_c] $stru_chain[$resi_c] $stru_resi[$resi_c]\n";
 #       }
 #       else{
 #               for ($miss_c = $miss_mark; $miss_c < @miss_num; $miss_c++){
 #               	if (($stru_num[$resi_c-1] < $miss_num[$miss_c]) && ($stru_chain[$resi_c-1] eq $miss_chain[$miss_c])){
 #                      		print "$miss_num[$miss_c] $miss_chain[$miss_c] $miss_resi[$miss_c]\n";
 #       			$miss_mark = $miss_c+1;
 #                      }	
 ##               }
 #       	print "TER\n";
 #
 #       }
 #
 #}
 #for ($miss_c = $miss_mark; $miss_c < @miss_num; $miss_c++){
 #	if ($stru_num[$resi_c] < $miss_num[$miss_c]){
 #       	print "$miss_num[$miss_c] $miss_chain[$miss_c] $miss_resi[$miss_c]\n";
 #       }	
 #}
 #
 #
 #
	print  "#################################################\n";
	print  "#                  By Liu Qing                  #\n";
	print  "# University of Science and Technology of China #\n";
	print  "#################################################\n";

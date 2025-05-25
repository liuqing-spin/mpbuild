#make internal non-standard residue
perl make_nonaa_inr_v5.pl -pdb CYA_min.pdb -mol2 CYA_nocap.mol2   -name CYA -gau g16 -np 10 -ac 0 
#make c-termianl non-standard residue
perl make_nonaa_ct_v5.pl -pdb CYA_min.pdb -mol2 CYA_nocap.mol2   -name CYA -gau g16 -np 10 -ac 0 
#make n-termianl non-standard residue
perl make_nonaa_nt_v5.pl -pdb CYA_min.pdb -mol2 CYA_nocap.mol2   -name CYA -gau g16 -np 10 -ac 0 

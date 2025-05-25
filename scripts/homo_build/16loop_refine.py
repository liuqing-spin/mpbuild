# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
import sys
import re 

#inpstr = input("input EM structure, target ID and start and end residue numbers of the loop:\n")
#inplst = inpstr.split()

input_model = sys.argv[2] + "_model_relax_align.pdb"
input_seq = sys.argv[2]
loop_indexA = sys.argv[3] + ":A"
loop_indexB = sys.argv[4] + ":A"

input_ssbond = sys.argv[5]
ssbond_list = []
with open(input_ssbond) as f:
    while True:
        line = f.readline()
        if not line: break
        line = line.rstrip('\r\n')
        spl_result = re.split(r'\s+',line)
        ssbond_list.append(spl_result[4])
        ssbond_list.append(spl_result[3])
        ssbond_list.append(spl_result[7])
        ssbond_list.append(spl_result[6])

#input_model = inplst[1] + "_model.pdb"
#input_seq = inplst[1]
#loop_indexA = inplst[2] + ":A"
#loop_indexB = inplst[3] + ":A"

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Create a new class based on 'LoopModel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(LoopModel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # 10 residue insertion 
        return Selection(self.residue_range(loop_indexA, loop_indexB))
    def special_patches(self, aln):
        for idx in range(0,len(ssbond_list),4):
            #disuA = str(ssbond_list[idx]) + ':' + ssbond_list[idx+1]
            disuA = str(ssbond_list[idx]) + ':A'
            #disuB = str(ssbond_list[idx+2]) + ':' + ssbond_list[idx+3]
            disuB = str(ssbond_list[idx+2]) + ':A'
            self.patch(residue_type='DISU', residues=(self.residues[disuA],self.residues[disuB]))


m = MyLoop(env,
           inimodel=input_model, # initial model of the target
           sequence=input_seq)          # code of the target

m.loop.starting_model= 1           # index of the first loop model 
m.loop.ending_model  = 20         # index of the last loop model
m.loop.md_level = refine.slow # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

m.make()

loop_pml = sys.argv[2] + '_loop_check.pml'
loop_pml_fo = open(loop_pml,'w')

#output_line = 'load ' + '../' + sys.argv[1] + '.pdb\n'
output_line = 'load ' + './' + sys.argv[1] + '.pdb\n'
loop_pml_fo.write(output_line)

for x in range(m.loop.starting_model,m.loop.ending_model+1):
    swr_list = []
    for i in list(str(x)):
        swr_list.append(i)
    zero_len = 4 - len(swr_list)
    zero_need = '0' * zero_len
    output_line = 'load '+ sys.argv[2] + '.BL' + zero_need + str(x) + '0001.pdb\n'
    loop_pml_fo.write(output_line)
    output_line = 'align ' + sys.argv[2] + '.BL' + zero_need + str(x) + '0001, ' + sys.argv[1] + '\n'
    loop_pml_fo.write(output_line)
    output_line = 'save loop_refine_align_' + str(x) + '.pdb, ' + sys.argv[2] + '.BL' + zero_need + str(x) + '0001' + '\n'
    loop_pml_fo.write(output_line)

loop_pml_fo.write('show cartoon\nhide lines\nzoom\nquit\n')

#list(str(sys.argv[0]))
print(loop_pml_fo.name)
print ( "#################################################\n")
print ( "#             Edited By Liu Qing                #\n")
print ( "# University of Science and Technology of China #\n")
print ( "#################################################\n")

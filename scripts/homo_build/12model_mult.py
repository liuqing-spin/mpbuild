from modeller import *
from modeller.automodel import *

import re
import sys
#inpstr = input("input target_sequence ID and ssbond file, seperated by whitespace:\n")
#inplst = inpstr.split()
#target_ID = inplst[0]
target_ID = sys.argv[1]
input_ali = target_ID + "-mult.ali"

template_list = []
with open(input_ali) as f:
    while True:
        line = f.readline()
        if not line: break
        line = line.rstrip('\r\n')
        if re.match('^>P1.*', line): 
            spl_result = line.split(';')
            #print(spl_result[1])
            template_list.append(spl_result[1])
#print(template_list)             

#input_ssbond = inplst[1]
input_ssbond = sys.argv[2]
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

#print(inplst)

log.verbose()

class MyModel(AutoModel):
    def special_patches(self, aln):
        for idx in range(0,len(ssbond_list),4):
            #disuA = str(ssbond_list[idx]) + ':' + ssbond_list[idx+1]
            disuA = str(ssbond_list[idx]) + ':A'
            #disuB = str(ssbond_list[idx+2]) + ':' + ssbond_list[idx+3]
            disuB = str(ssbond_list[idx+2]) + ':A'
            self.patch(residue_type='DISU', residues=(self.residues[disuA],self.residues[disuB]))



env = Environ()
a = MyModel(env, alnfile=input_ali,
              #knowns=(template_list[0],template_list[1],template_list[2]), sequence=template_list[3])
              knowns=(template_list[0],template_list[1]), sequence=template_list[2],
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
#if int(inplst[2]) > 5:
#    a.ending_model = int(inplst[2])
if int(sys.argv[3]) > 5:
    a.ending_model = int(sys.argv[3])
else:
    a.ending_model = 5

a.make()

print ( "#################################################\n")
print ( "#             Edited By Liu Qing                #\n")
print ( "# University of Science and Technology of China #\n")
print ( "#################################################\n")

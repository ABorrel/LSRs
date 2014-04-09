
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline


#result_handle = NCBIWWW.qblast("blastp", "pdb", "/home/borrel/Yue_project/dataset/AMP/4EPM/4EPM.fasta")



blastp_cline = NcbiblastpCommandline(query="/home/borrel/Yue_project/dataset/AMP/4EPM/4EPM.fasta", db="pdbaa", outfmt=5, out="tt.xml")


print(blastp_cline)

stdout, stderr = blastp_cline()

result_handle = open("tt.xml")


blast_records = NCBIXML.read(result_handle)
# blast_records = list(blast_records)
print vars (blast_records)
print vars(blast_records.alignments[0])
print vars(blast_records.alignments[0].hsps[0])
#print result_handle

# save_file = open("my_blast.xml", "w")
# save_file.write(result_handle.read())
# save_file.close()
# blast_records = NCBIXML.parse(result_handle)
# blast_record = NCBIXML.read(result_handle)
# 
# E_VALUE_THRESH = 0.04
# 
# for alignment in blast_record.alignments:
#     for hsp in alignment.hsps:
# #         if hsp.expect < E_VALUE_THRESH:
#             print('****Alignment****')
#             print('sequence:', alignment.title)
#             print('length:', alignment.length)
#             print('e value:', hsp.expect)
#             print(hsp.query[0:75] + '...')
#             print(hsp.match[0:75] + '...')
#             print(hsp.sbjct[0:75] + '...')
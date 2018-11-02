import mbuild as mb
import scripts.xml_utils as xml_utils 

dopc = mb.load('dopc.gro')
dopc.name = 'DOPC'
chl = mb.load('new_chl1.gro')
chl.name = 'CHL1'

new_dopc = xml_utils.align_cmpd(dopc, [0, 87])
new_dopc.name = dopc.name
new_chl = xml_utils.align_cmpd(chl, [1,64])
new_chl.name = chl.name

new_dopc.save('aligned_dopc.gro', residues='DOPC')
new_chl.save('aligned_chl.gro', residues='CHL1')

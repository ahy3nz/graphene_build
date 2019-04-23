#import mvp
import mbuild as mb
from mbuild.lib.moieties.ch3 import CH3
import pdb

# A monomer of PVP
class mvp(mb.Compound):
    def __init__(self):
        super().__init__(name='mvp')

        n = mb.Particle(name='N')
        c1 = mb.Particle(name='C')
        c2 = mb.Particle(name='C')
        c3 = mb.Particle(name='C')
        c4 = mb.Particle(name='C')
        o41 = mb.Particle(name='O')
        c5 = mb.Particle(name='C')
        c6 = mb.Particle(name='C')


        self.add([n, c1, c2,  c3,  c4,
            o41, c5, c6 ])

        self.add_bond([n, c1])
        self.add_bond([c1, c2])
        self.add_bond([c3, c2])
        self.add_bond([c3, c4])
        self.add_bond([c4, o41])
        self.add_bond([c4, n])
        self.add_bond([n, c5])
        self.add_bond([c5, c6])

        self.add(mb.Port(anchor=c5), 'up')
        self.add(mb.Port(anchor=c6), 'down')

# Implicit hydrogens
# XYZ doesn't matter, we just need bonds and particle names
# So we can upload to charmm-gui
cap_with = 'H'
pvp = mb.Polymer(mvp(), n=2)

if cap_with=='C':
    c_top = mb.Particle(name='C')
    cap_top = mb.Compound()
    cap_top.add(c_top)
    cap_top.add(mb.Port(anchor=c_top), 'down')
    mb.force_overlap(cap_top, cap_top['down'], pvp['up'])

    c_bot = mb.Particle(name='C')
    cap_bot = mb.Compound()
    cap_bot.add(c_bot)
    cap_bot.add(mb.Port(anchor=c_bot), 'up')
    mb.force_overlap(cap_bot, cap_bot['up'], pvp['down'])
    pvp.add(cap_top)
    pvp.add(cap_bot)

elif cap_with =='H':
    pass
    
# Override the parmed structure so we can actually get double bonds
struc = pvp.to_parmed()
for bond in struc.bonds:
    if 'O'==bond.atom1.name or 'O' == bond.atom2.name:
        bond.order = 2.0
struc.save('pvp.pdb',overwrite=True)
struc.save('pvp.mol2',overwrite=True)

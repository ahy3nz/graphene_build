import subprocess

import itertools

### this is a [bad] wrapper around the load_0.bash just so i can use itertools
### to make this load_2, load_4, or load_6, modify the _exec_bash my_command

def main():
    angles = ['15', '30', '45']
    #angles = ['0', '15', '30', '45']
    forces = ['125']
    #forces = ['250', '500', '1000']

    trials = ['a', 'b', 'c']
    for angle, force in itertools.product(angles,forces):
        for trial in trials:
            output_dir = 'alt6_dna/k{0}_{1}_{2}'.format(force, angle, trial)
            correct_top = '6_dna.top'
            GF = 'generated-files'
            _exec_bash(GF, angle, force, output_dir, correct_top)


def _exec_bash(GF, angle, force, output_dir, correct_top):
    my_command = """
export GF={GF}
export ANGLE={angle}
export FORCE={force}
export output_dir={output_dir}

export correct_top={correct_top}
mkdir -p ${{output_dir}}
mkdir -p ${{GF}}

# make 2nd sheet for top
gmx editconf -f inputs/dna.gro -rotate 0 0 180 -o ${{GF}}/dna-top2.gro
gmx editconf -f ${{GF}}/dna-top2.gro -translate 6 6 0 -o ${{GF}}/dna-top2.gro
#
gmx editconf -f inputs/dna.gro -rotate 0 0 90 -o ${{GF}}/dna-top3.gro
gmx editconf -f ${{GF}}/dna-top3.gro -translate 6 0 0 -o ${{GF}}/dna-top3.gro
#
#
## make 1st sheet for bottom
gmx editconf -f inputs/dna.gro -rotate 180 0 0 -o ${{GF}}/dna-bottom1.gro
gmx editconf -f ${{GF}}/dna-bottom1.gro -translate 0 6 6 -o ${{GF}}/dna-bottom1.gro
#
## make 2nd sheet for bottom
gmx editconf -f ${{GF}}/dna-bottom1.gro -rotate 0 0 180 -o ${{GF}}/dna-bottom2.gro
gmx editconf -f ${{GF}}/dna-bottom2.gro -translate 6 6 0 -o ${{GF}}/dna-bottom2.gro
#
gmx editconf -f ${{GF}}/dna-bottom1.gro -rotate 0 0 90 -o ${{GF}}/dna-bottom3.gro
gmx editconf -f ${{GF}}/dna-bottom3.gro -translate 6 0 0 -o ${{GF}}/dna-bottom3.gro

#
## combine the gro files
export GIM='gmx insert-molecules'
export OUT='generated-files/graphene-dna.gro'
export ROT='-rot none'
export DR='-dr 0.0'
export IP='-ip inputs/positions.dat'
export MOL='-nmol 1'
export SC='-scale 0.00001'
$GIM -f inputs/sheet.gro -ci inputs/dna.gro -o $OUT $MOL $ROT $DR $IP $SC
#for input in ${GF}/dna-top2.gro ${GF}/dna-bottom1.gro ${GF}/dna-bottom2.gro ${GF}/dna-bottom3.gro
for input in ${{GF}}/dna-top2.gro ${{GF}}/dna-top3.gro ${{GF}}/dna-bottom1.gro ${{GF}}/dna-bottom2.gro ${{GF}}/dna-bottom3.gro
do
$GIM -f $OUT -ci $input -o $OUT $MOL $ROT $DR $IP $SC
done

# at this point we have the gnf-dna complex, now we need to add it to the bilayer
# but first, we need to change the box size of the gnf-dna complex to match the bilayer
bilayer_box=`tail -n 1 inputs/bilayer-no-water.gro | awk '{{print $1, $2, $3'}}`
gmx editconf -f ${{GF}}/graphene-dna.gro -box ${{bilayer_box}} -o ${{GF}}/graphene-dna-resized.gro
gmx editconf -f ${{GF}}/graphene-dna-resized.gro -rotate 90 45 45 \
    -o ${{GF}}/graphene-dna-resized.gro -c
gmx editconf -f ${{GF}}/graphene-dna-resized.gro -translate 0 0 2.5 \
    -o ${{GF}}/graphene-dna-resized.gro

cd inputs
python graphene_rotations.py --angle $ANGLE --force $FORCE --gro ../${{GF}}/graphene-dna-resized.gro --out ../${{GF}}/spun.gro
cd ..
# now move the sheet as necessary and insert into bilayer system
$GIM -f inputs/bilayer-no-water.gro \
    -ci ${{GF}}/spun.gro \
    -o ${{GF}}/bilayer-graphene-dna.gro \
    $ROT $MOL $DR $IP $SC

# add water; copy the topol file incase we need to run again
#cp inputs/topol-bilayer-graphene-dna.top topol.top
cp inputs/${{correct_top}} topol.top
gmx solvate -cp ${{GF}}/bilayer-graphene-dna.gro \
    -cs inputs/water.gro \
    -o ${{GF}}/bilayer-graphene-dna-water.gro \
    -p topol.top \
    -scale 0.57
    #-scale 1.2

# now need to add ions, but to do that you need a tpr file, so we need to run grompp
export MDP='inputs/grompp.mdp'
export GRO=${{GF}}/bilayer-graphene-dna-water.gro
export TOP='topol.top'
export TPR=${{GF}}/'no-ions.tpr'
# ignore warning about generating velocities
gmx grompp -f ${{MDP}} -c ${{GRO}} -p ${{TOP}} -o ${{TPR}} -maxwarn 2

# finally, add ions and re-grompp
export TPR_IN='-s generated-files/no-ions.tpr'
export TOP='-p topol.top'
export OUT='-o generated-files/bilayer-graphene-dna-water-ions.gro'
export NDX='-n index.ndx'
export PNAME='-pname NA'
export NEUT='-neutral'
export MDP='inputs/em.mdp'
export TPR='em.tpr'
echo q | gmx make_ndx -f $GRO
echo SOL | gmx genion $TPR_IN $TOP $OUT $NDX $PNAME $NEUT
echo q | gmx make_ndx -f ${{GF}}/bilayer-graphene-dna-water-ions.gro

# re-grompp, ignoring warning about generating velocities
gmx grompp -f ${{MDP}} -c ${{GF}}/bilayer-graphene-dna-water-ions.gro ${{TOP}} -o ${{TPR}} -maxwarn 1

# visualize
#vmd bilayer-graphene-dna-water.gro

# remove old backup files
rm \#*
rm ${{GF}}/*
cp inputs/angled_insertion.mdp .
cp *.* ${{output_dir}}
""".format(**locals())
    with open('make_all.log','w') as f:
        p = subprocess.Popen(my_command, shell=True, stdout=f, stderr=f)
        p.wait()
if __name__ == "__main__":
    main()

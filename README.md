# ultra_enumerate
enumerate a large amount of rdkit valid SMILES.

Sometime, the doRandom=True flag of Chem.MolToSmiles in rdkit just doesn't cut it. For example, the topical nasal decongestant `c1ccccc1CC(C)NC` has merely 40 unique random SMILES. Making use of the "use ring operators as normal bonds" trick, we can do better than that and generate a large amount of ugly but rdkit-valid SMILES! Because using this trick makes it possible to use any atom ordering, this means a total of `N!` SMILES of this type are possible, with `N` the amount of atoms and `!` the factorial (I'm not shouting, I'm calm). In case there's symmetry the number will be a bit lower. 

## reqs

should work in any python env with reasonably up to date rdkit.

## usage
`python ultra_enumerate.py --smiles "c1ccccc1CC(C)NC" -o output.smi -N 100000`

with `-N` the amount of SMILES you want, `-o` the output (if left empty it will print to stdout) and `--smiles` the SMILES you want to have variants of.

Note that they will look like this.
```
c:1:2.C345.c:6:7.C4.c:8:1.C93.N5%10.c:%11:29.C%10.c:6:%11.c:7:8
C1.C2.c:3:4.c:5:6.c:6:7.C819.c:3:%10.c:4:5.c:%10:7%11.N92.C%118
C123.c:4:5.c:6:7.C2.C8.c:9:%10%11.N38.c:4:9.c:5:6.C%111.c:7:%10
c:1:2.c:3:4.C567.c:8:49.C6.c:2:%10.c:1:8.C95.C%11.N7%11.c:%10:3
c:1:2.c:3:4.N56.c:7:1.c:4:28.c:9:7.C6.c:3:9.C%10.C8%11.C%11%105
c:1:2.c:3:4.c:4:56.c:7:5.c:2:7.C89%10.N%10%11.C9.c:3:1.C68.C%11
C123.c:4:5.C2.c:6:7.C8.c:5:9%10.c:4:6.N38.c:7:%11.c:%11:9.C%101
c:1:2.C3.c:4:5.c:1:6.c:6:4.C789.C%107.c:2:%11%10.C8.N93.c:5:%11
c:1:23.c:4:5.c:6:7.c:7:2.c:4:1.C8.c:5:6.C39.N%10%11.C%11.C98%10
c:1:23.N45.c:6:7.C894.c:7:2.c:%10:%11.C5.C38.c:%10:1.c:%11:6.C9
```

The speed is relatively practical, for example the output 100000 nasal decongestant random SMILES take about 2 seconds to generate.

## limitations
For the moment stereochemistry is ignored.

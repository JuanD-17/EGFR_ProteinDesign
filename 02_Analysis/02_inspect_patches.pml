# Visual check of the top patches. Numbering here is PDB.
load ../01_Target/structure/6ARU.pdb, egfr
hide everything
show cartoon
color grey80
color palecyan, chain A
select fab, not chain A
color salmon, fab

select patch1, chain A and resi 342+403+405+406+407+408+409+410+411+434
show surface, patch1
color orange, patch1

select patch2, chain A and resi 407+408+409+410+411+434+436+459+463+464+465
show surface, patch2
color marine, patch2

select patch3, chain A and resi 384+406+407+408+409+410+411+436
show surface, patch3
color forest, patch3

select patch4, chain A and resi 440+449+450+463+464+465+466+469+472
show surface, patch4
color purple, patch4

select patch5, chain A and resi 342+400+403+405+406+407+431+433+434+458+459
show surface, patch5
color yellow, patch5

set transparency, 0.2
bg_color white
orient

# UniProt position = PDB position + 24

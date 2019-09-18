;Energy minimize initial configuration

gmx grompp -f em.mdp -p topol.top -c solvated.gro -o em
gmx mdrun -v -deffnm em

;Run short equilibration with berendsen pressure control:

gmx grompp -f berendsen.mdp -p topol.top -c em.gro -o berendsen
gmx mdrun -v -deffnm berendsen

;Run production simulation switching to Parrinello-Rahman pressure control

gmx grompp -f npt.mdp -p topol.top -c berendsen.gro -o npt
gmx mdrun -v -deffnm npt

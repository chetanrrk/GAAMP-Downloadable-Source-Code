mol.prm          - GAFF parameters file (equilibrium bond length and angle may be adjusted)
mol-tor.prm      - Parameter file with fitted torsion parameters
mol_final.prm    - Parameter file with fitted parameters and CHARMM-GUI compatible atom types

mol-esp.rtf      - GAFF rtf file with fitted charges from ESP fitting
mol-esp-wat.rtf  - Rtf with fitted charges from ESP and water-compound interactions fitting
mol-tor.rtf      - Rtf with fitted charges and updated after torsion fitting
mol_final.rtf    - Rtf with fitted charges and updated after torsion fitting and CHARMM-GUI compatible atom types

mol-opt.crd      - Crd file with optimized geometry from given structure at HF/6-31G* level

report-esp.txt   - The result of ESP charge fitting
report-esp-wat.txt - The result of charge fitting with ESP and water-compound interactions
result-1D.html   - The result of 1D energy profiles after 1D torsion fitting or 1D torsion and rotamer energies fitting (if there are more than one soft dihedrals)
fitting-1d-*.dat - The 1D energy profile for a soft dihedral after 1D torsion fitting. Format: phi, E_QM, E_MM
1d-qm-mm-*.dat   - The 1D energy profile for a soft dihedral after 1D torsion and rotamer energies fitting (if there are more than one soft dihedrals). Format: phi, E_QM, E_MM
rotamer-E.dat    - The rotamer energies after 1D torsion and rotamer energies fitting. Format: index, E_QM, E_MM
fit-mol.conf     - The configuration file used for charge fitting

qm               - The directory containing all QM data used in parameter fitting
qm/qm-mol-opt.out - The Gaussian output of QM geometry optimization
qm/qm-mol-esp.out - The Gaussian output of QM electrostatic potential used for non-polarizable model
qm/qm-mol-wat-donor*.gjf - The Gaussian input of the QM calculation to determine Emin and Rmin for a H donor interacting with a water molecule
qm/qm-mol-wat-donor*.out - The Gaussian output of the QM calculation to determine Emin and Rmin for a H donor interacting with a water molecule
qm/qm-mol-wat-donor*.pdb - The snapshot in which the pose of a water is optimized in QM for a H donor interacting with a water molecule
qm/qm-mol-wat-acceptor*.gjf - The Gaussian input of the QM calculation to determine Emin and Rmin for a H acceptor interacting with a water molecule
qm/qm-mol-wat-acceptor*.out - The Gaussian output of the QM calculation to determine Emin and Rmin for a H acceptor interacting with a water molecule
qm/qm-mol-wat-acceptor*.pdb - The snapshot in which the pose of a water is optimized in QM for a H acceptor interacting with a water molecule
qm/qm-scan-*-*.gjf  - The Gaussian input of QM 1D torsion scan for a specific soft dihedral
qm/qm-1d-phi-*.out  - The Gaussian output of QM 1D torsion scan for a specific soft dihedral
qm/qm-rotamer-*.gjf - The Gaussian input of QM geometry optimization for a specific rotamer
qm/qm-rotamer-*.out - The Gaussian output of QM geometry optimization for a specific rotamer

mm-acceptor-*.pdb   - The snapshot in which the pose of a water is optimized in MM for a H acceptor interacting with a water molecule
mm-donor-*.pdb      - The snapshot in which the pose of a water is optimized in MM for a H donor interacting with a water molecule


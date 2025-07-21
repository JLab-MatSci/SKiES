In this example basic features of SKiES code are presented.

The examples for two materials Ag and Pd are given. While Ag has a linear dependence of
electrical resistivity, Pd demostrates a non-linear dependence at high temperatures and 
this behavior can only be described with the use of general formulas of Allen's method
(this is connected with the differences in electronic DOS of Ag and Pd near the Fermi level,
see the figures Ag_resistivity.png, Pd_resistivity.png, DOS_cmp.png in this folder).

The subfolders with the "_not_ready" suffix do contain all the files necessary for
launching preliminary Quantum ESPRESSO/EPW calculations. See the corresponding README.md
files for the details.

The other folders (with the "_ready" suffix) hold the files 'crystal.fmt', 'epwdata.fmt'
and 'vmedata.fmt' which are the byproducts of the preliminary Quantum ESPRESSO/EPW
calculations which are described above. The only file that is missing is the binary
heavy file with the suffix '.epmatwp' containing e-ph matrix elements. One has to
download it from the remote repository as described in the corresponding README.md files.
After this file is downloaded one is able to run the SKiES test calculations in the subfolders
'skies_bands', 'skies_phonons', 'skies_transport_lowT' and 'skies_transport_general'
as described in details in the bash scripts given and in the manual.pdf in the root
folder of SKiES.

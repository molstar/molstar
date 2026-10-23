# Periodic trajectories

MD runs in a periodic box. Engines keep atoms (or molecule centers) inside the primary box, so a molecule that crosses
a face is written partly on one side and partly on the other. Mol\* draws the frame as given: bonds and cartoon segments
stretch across the whole box ("portalling"). This folder holds the analysis ([analysis.md](analysis.md)) and the plan
([plan.md](plan.md)).

Status: phases 0 to 4 done (readers, molecules whole, in the box or centered on a selection, no jump, the Prism command and bar menu, the ensemble structure); MVS and LAMMPS image flags are open questions. Open questions for the maintainers are at the end of
[plan.md](plan.md#open-questions).

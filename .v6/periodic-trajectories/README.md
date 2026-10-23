# Periodic trajectories

MD runs in a periodic box. Engines keep atoms (or molecule centers) inside the primary box, so a molecule that crosses
a face is written partly on one side and partly on the other. Mol\* draws the frame as given: bonds and cartoon segments
stretch across the whole box ("portalling"). This folder holds the analysis ([analysis.md](analysis.md)) and the plan
([plan.md](plan.md)).

Status: phases 1 and 2 done (molecules whole; in the box, or centered on a selection), and the GRO triclinic box; no jump, LAMMPS boxes and the Prism command are open. Open questions for the maintainers are at the end of
[plan.md](plan.md#open-questions).

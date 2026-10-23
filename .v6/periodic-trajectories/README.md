# Periodic trajectories

MD runs in a periodic box. Engines keep atoms (or molecule centers) inside the primary box, so a molecule that crosses
a face is written partly on one side and partly on the other. Mol\* draws the frame as given: bonds and cartoon segments
stretch across the whole box ("portalling"). This folder holds the analysis ([analysis.md](analysis.md)) and the plan
([plan.md](plan.md)).

Status: done (phases 0 to 5 and the distance links), decisions at the end of [plan.md](plan.md#decisions). Open:
periodic images, NAMD `.xsc` origins.

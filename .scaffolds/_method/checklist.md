# Checklist: <work id>

<!-- The template for a work folder's checklist, `execution/<work-id>/checklist.md` (§plan-and-checklist): one line per
     item of the plan, ticked with its commit as it lands (`- [x] ... (<commit>)`); the loop picks its next items here
     and Finalize starts from the unticked ones. Nothing else goes here: no notes, no dates, no findings (`state.md`,
     `findings.md`). An item cut by the author's word is ticked with "cut (the author, <date>)"; one carried to the next
     plan, with "carried". A lane is ticked with its merge commit, after its look and review (loop, Each step). The
     closed-work check fails a closed folder with an unticked line. Delete this comment once filled. -->

<!-- Stage 0: designs -->
- [ ] <draft>, its mock and its review
- [ ] GATE 0: the drafts agreed, with their mocks; the plan agreed (the author's words in the tick's commit)

<!-- Stage <n>: <name> -->
- [ ] **<LANE>** (`briefs/<lane>.md`, <model>) <what>
- [ ] **<LANE>** (`briefs/<lane>.md`, <model>) <what>
- [ ] GATE <n>
- [ ] **<fix lane>** (`briefs/<fix lane>.md`, fast) `fixes.md`'s open items after stage <n>

<!-- Last stage: use, measure, finalize -->
- [ ] **<the novice pass>** after the last lane that draws, its bugs reproduced before they are placed
- [ ] **<the measurements>** <the representatives at full size, the bundles against their budgets>
- [ ] GATE: last
- [ ] FINALIZE: the unticked items run, cut or carried; reconcile; `foldback.md` with the author's decision; the gate;
      the record; `state.md`'s closing line; the next plan

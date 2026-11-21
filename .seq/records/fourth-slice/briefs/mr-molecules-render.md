# Lane MR: molecules drawn and hit

**Task**: build the render half of [molecules](../../../design/molecules.md) (Decision: On an axis, A section of a complex, Behaviors; the matrix cells; Consequences `view/`, `render/`, `seqx-viewer`): a section of a complex (`classes`, `form`, chip labels), the `letters` mark and component code (`modified`, `code`; `letters.modifiedMark`), residue labels of a molecule and element symbols on atoms, the pairs layer's ends by row, a `ranges` layer over pairs, hover text with component, alternatives, atom name and residue, the hover record's identifiers, `links.hover` (off by default, on in examples), `block.form` in the context menu; ligand chips flowing in text (EC-99) and a block's label inline on its first text line (EC-95). Labels follow [layouts](../../../design/layouts.md) (Labels: inside if they fit, else beside; leaders never cross; a run of one label once). The Design page `Seq* (next)/Design/Molecules`.

**Output**: render and viewer changes with tests: drawing against the data (every modified position in view has its mark, every chip its label, every link with both ends in view its arc and hit box), a twin per consumer, a quick representative of the frame target with 200 merged copies and a glycan chip (EC-67); `_spec` Resolved view, Render, viewer HoverModel and Behaviors.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/mr-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.

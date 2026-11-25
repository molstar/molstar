# Lane MC: the mmCIF runtime

**Task**: build `runtimes/mmcif` (opt-in entry of `seqx-core`) as [molecules](../../../design/molecules.md) (Reading, Configuration: the allowlist, `assembly`, `water`, `name`) and [seq-data-ports](../../../design/seq-data-ports.md) (Seq io) say: CIF categories read through `@molstar/io`'s allowlist (add its entries; a hot-file change to list), converted to a `StructureDescription` and then `structureNodes` by DP's code (not a private converter), `values` by `field` (B-factor, pLDDT), an altloc rule, hydrogens refused. The conformance suite splits as the design says: the CIF reading against counts recorded from the file (every entity, chain, branch link, `struct_conn` row and author number; an entity without canonical letters; an altloc) on 1EHZ, 1FIN, 4HHB and a glycoprotein; `structureNodes` on descriptions. Fixtures re-recorded. A size row for the entry (under 30 kB, loaded only when composed). `pnpm build` and `check:workspace` after (`@molstar/io` changed).

**Output**: the runtime, its suite, the fixtures and their record script, the size row.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/mc-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.

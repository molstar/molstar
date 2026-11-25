# Lane SV: services in Storybook

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). The first fetch of each service's fixtures is allowed (the author, 2026-11-21: "download what you need for fixtures etc"): record each once into the repo, trimmed to what the stories need, with its source URL, date and licence beside it (UniProt CC-BY 4.0, PDBe SIFTS EMBL-EBI terms, AlphaFold DB CC-BY 4.0); tests never reach the network. Checks before the report: as `fx5-author-look.md` lists them; the look on port 6010.

**Task**: compose UniProt, PDBe SIFTS and AlphaFold DB as application-layer runtimes in the Storybook ([examples](../../../design/examples.md), Showcases; Gate 0, 6 and 12) from recorded fixtures with a live switch through the asset policy's allowed hosts (tests read fixtures only, never the network); the predictions-onto-a-structure and addressing showcases on ST3's fake viewer; a user's own file dropped onto a view for the formats built (`seq`, FASTA, Stockholm, the records, mmCIF through MC; no GenBank). SIFTS builds `structure` maps as MK does. Each with a headless twin on a fake clock.

**Output**: the runtimes, fixtures, stories, twins, a drop target.

**Model**: fast (sonnet).

**Where**: A worktree on a branch `seq/sv-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.

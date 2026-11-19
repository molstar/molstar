# Lane SR: finalize's reconcile, prepared (read-only)

**Goal** (the author, 2026-11-19: "Spec reconcile prep"): most of `loop.md`'s Finalize step 1 found before the slice's
last lanes land, so finalize applies a list instead of searching. **Read only**: you edit nothing; your report is the
list. Lanes never edit `.seq/`, and RK and CS are changing `packages/seq` beside you.

Compare what is built on `packages-seq` (the code, its tests, its stories) with what the documents say, and list each
disagreement:

1. **The `_spec`s' Configuration tables** (`packages/seq/core/_spec/README.md`, `packages/seq/viewer/_spec/README.md`):
   every config item, spec param, `RenderConfig`, `SpecBounds`, `ViewerConfig` field, port and call option the code
   has, against each design's Configuration table (`.seq/design/*.md`, "## Configuration"): missing rows, rows with a
   wrong default or "whose", rows for what was removed.
2. **`.seq/design/features.md`**: each row's status (`done`, `building: …`, `planned`, `deferred`) against the code: a
   row marked planned that is built, a `building` row whose remainder is now done or is not.
3. **`.seq/status.md`**: what it says is built or open, against the code and `.seq/execution/third-slice/state.md`.
4. **`.seq/design/edge-cases.md`**: each EC this slice claims (the third slice's plan, `plan.md`) has the outcome and the
   test it names; a named test file that does not exist.
5. **`drift/`** (one note, `themes-spellings.md`): resolved by the code or not.

For each finding: the file and line, what it says, what the code does (file and line), and the edit you propose.
Mark what RK (the keyboard cursor's rows) or CS (p53's chromosome) will change, from their briefs in this folder, as
"after RK" or "after CS" rather than proposing an edit now.

**Model**: fast (sonnet). **Where**: the main checkout, read only (no worktree, no commands that write: no installs,
builds or formatters; `git log`, `grep`, reading files and `pnpm test:unit` with no `--update` are fine).

**Stop when** the five lists are done, or after about an hour's work: report what you have, by list.

**Report**: the five lists; a count by list; anything that looks like a design question rather than a stale document
(those become drift notes).

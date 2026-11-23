# Lane <ID>: <what it does, in a few words>

<!-- The template for a lane's brief, `execution/<work-id>/briefs/<id>-<slug>.md` (§briefs-name-checks). With the rule
     sections and designs it names, it is the lane's whole context: a lane reads those, not the whole scaffold. Short:
     a page. Delete the comments once filled. -->

**Model**: <strong | fast> (<why, by scope: a closed item list and named outputs run fast; a design's core runs
strong>). **Base**: `<integration branch>` at the orchestrator's head (`git reset --hard <integration branch>` first).
Branch `<prefix>/<lane>`. <What runs beside it, or "Nothing beside you".>

**Why**: <the plan's stage, the design, the author's words with the date; for application-level work (a composition
on the public surface, a showcase), say so and name the library parts it uses and what the library had to gain>.

**Read first**: <the designs by file and section, their mock and precedent where the lane draws, the `_spec` sections,
the rule sections by id>.

**Items**, a commit at each:

<!-- Numbered, each with what proves it landed: a test by name, an output, a section of a spec. The orchestrator ticks
     them against that proof at the merge; an item with no proof is not landed. Name the edge cases each covers, so the
     lane's tests name them. A fix's item says how its test is shown failing on the parent. -->

1. **<item>**: <what, where>. Proof: <the test, the output, the section>.
2. ...

**How**:
- **host neutral**: <what the host owns that the lane never takes>;
- **two consumers**: a shared piece is tested with <its two consumers>, live in one test;
- **minimal code**: reuse <the base library's helpers> before writing new ones; a file grown past 800 lines is a
  question in the report, past 1500 a failure (§module-size);
- **its performance representative**: <one quick case; a bound under the heavy lane for a time; an allocation case
  when it adds to a frame's path>;
- **lazy parts awaited**: a test over a view with lazy parts awaits them before the gesture it checks;
- **the mock**: <the design's mock, when the lane draws a new part>; the built view matches it gesture by gesture;
- **a rewrite keeps its identifiers**: <when the lane rewrites text across many files (citations, ids, links)> list
  the identifiers the files cite (`EC-<n>`, rule ids, spec sections) at the base and at your head, and name each one
  gone with why; <the scaffold's id check in its compare mode, against the base, where it has one>.

**Where**: a worktree on `<prefix>/<lane>`. Touches <files>; out of scope <files, parts>. Headers <the repository's
rule>. Hot files never staged: <the lockfile, manifests, configs, the work folder, the changelog folder, the TODO>;
declare what you need in the report. <What stops the lane: a change past the agreed surface of a base library.> Write
only in your worktree and your own temporary space. Never touch the network.

**Checks**, each its own command, green only on exit 0 with a non-zero test count: <the tests at the tiers the plan
names, the type check, lint, format, the performance representative, the size check, the scaffold's checks>. Already red
at the base: name it in the report, don't fix it unless it is in scope.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under <n> words): commits; each item with its proof; each check with its exit code and test count;
assertions deleted, skipped or weakened, and every changed expectation with why; what you left out, refused or fell
back on; drift, gaps and followups (you don't write the work folder); the changelog text; what to look at; for a
rewrite across many files, the identifiers before and after and each one gone.

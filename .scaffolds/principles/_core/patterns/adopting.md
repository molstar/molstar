# Adopting the package {#adopting}

The procedure for setting the package up in a project ([bootstrap.md](../../bootstrap.md) is the entry point that
sends an agent here). Run it in order; each step ends with what to verify. Ask the
author where it says ask, and don't guess: a setup built on a guess spreads into every lane after it.

## 0. Ask first

Before copying anything, settle with the author; note the answers, and they go into the manifest when step 1 creates
it:

- **Where the package goes.** Usually `.principles/` at the root; for part of a larger repository, in that part's
  folder. Keep the package's names inside it, so its rules read the same.
- **Where the instance goes.** The project's own modules, manifest, instance skill, designs, sources, records and
  work folder live in one folder, `.<instance>/`, named after the repository (`.design/` becomes
  `.<instance>/design/`, and so on). Inside a dot-folder, no subfolder needs a dot to hide it; a system folder keeps
  its underscore by the `_` rule ([system.md](../system.md)), so the instance skill lives in
  `.<instance>/_skills/`. The work folder is `.<instance>/execution/`; the merge check assumes that path. The
  package and the instance never share a path: the instance's manifest is `.<instance>/principles/README.md`, so the
  package is never put at `.<instance>/principles/`.
- **What it governs.** Which paths and packages follow the principles, and which follow the host's own rules.
- **Who decides.** The authors whose call a design or rule change is (§shared-ownership); with more than one, which
  areas are whose.
- **What the host already has.** Its own `AGENTS.md`, agent-context folders, checks, test tiers, layout, tracker. Where
  the host and the package speak about the same thing, ask which wins; the usual answer is the host.
- **Which modules.** `_core` always; the others by what the project builds (the manifest's "Read before" column).
- **The section sign.** The id check claims `§` followed by a name across the whole repository, as a citation that
  must resolve. Does the host use it otherwise (legal text, generated docs)? Those files are exempted in step 4.

## 1. Create the instance, copy and map

- **Create the instance**: the folder `.<instance>/` and its manifest, `.<instance>/principles/README.md`, from
  [_template/instance-manifest/README.md](../../_template/instance-manifest/README.md). Write step 0's answers into
  it now; the rest of this step fills in its sections.
- **The adoption is itself work**, with the work id `adopt-principles`: create its work folder,
  `.<instance>/execution/adopt-principles/`, with a plan of these steps, a checklist and a `state.md` (§plan-and-checklist, §state-file), and add
  `.tmp/` to `.gitignore` (the work folder is committed, §permanent-record). Its session record comes at the end (bootstrap, step 4).
- Copy `.principles/` only. Not the home repository's designs, records, sources or prototypes: those are its own
  instance, not the package.
- **Record the provenance** in the manifest: where the package came from and the commit or date it was copied at,
  so a later update from it is a merge from a known base.
- **Map the folders** in one manifest table when they live elsewhere:

  | The package says | Here |
  |---|---|
  | `.principles/`, `.design/`, `.sources/`, `.records/` | ... |
  | `.execution/`, `.tmp/` | `.<instance>/execution/` (committed, never cited; the merge check assumes this work folder), `.tmp/` (gitignored) |
  | `_spec/` | ... |
  | the checks: id, path, merge | their commands |

- **Write the adaptations**, one line each with why: modules dropped; rules where the host wins (its layout, its test
  tiers, its package rules); where open items go if not to issues (§followups).
- **Bootstrap patterns** (`<module>/bootstrap-patterns/`) are starting implementations, to take in where they fit
  ([system.md](../system.md#code-patterns)). One taken in as written stays where it is: the package's copy serves as
  the project's stable pattern, and an adaptation line says so. One the project changes is copied into the instance as
  `.<instance>/principles/<module>/stable-patterns/<new-id>/` with a new id, and its adaptation line says it
  supersedes §<bootstrap-id> here. The project's own new patterns start as sketches, in
  `.<instance>/principles/<module>/sketch-patterns/`.
- **Optionally, the full tree**: the project's own codebase tree, with the rule ids behind each folder, at
  `.<instance>/principles/layout/tree.md` (the package's `layout/tree.md` (not adopted here) is the minimal one).
- **Verify:** the manifest exists at `.<instance>/principles/README.md`, every row of its module tables points at
  a folder that exists, and `git check-ignore .tmp/x` names it while `git check-ignore .<instance>/execution/x` names nothing.

## 2. Coexist with the host

- Add a pointer to the package from the host's root `AGENTS.md` (or its equivalent), and an `AGENTS.md` for the
  governed area that says where things are, how work happens, and that the host's rules apply too.
- List the host's other agent-context folders there, so no agent mistakes them for the package's.
- **Verify:** an agent starting in the governed area finds the package from the root `AGENTS.md` in one hop.

## 3. Link the skills

- Link each skill folder: `.agents/skills/<name>` → `<package>/_skills/<name>`, and `.claude/skills` → `.agents/skills`
  where the harness reads that path. Three skills: `principles`, `principles-router`, `principles-self-check`.
- **Add the instance skill**, named after the repository: copy `_template/instance-skill/SKILL.md` to
  `.<instance>/_skills/<instance>/SKILL.md`, fill it in (what the repository is, the package's provenance, the
  manifest's path), and link it as `.agents/skills/<instance>`. It is what a session here loads first; it points at
  the instance's manifest, which the router reads.
- **Verify:** the agent harness lists all four in a fresh session.

## 4. Wire the checks

A rule a tool can check is a check from day one. Each runs as its own command and gates on its exit code
(§green-is-exit-code). The package ships the first five in `_checks/` ([the checks](../../_checks/README.md): how to
run them, their arguments, the planted error each must fail on); pass them the instance (`--instance .<instance>`).

- **The id check** (`_checks/principles-ids.ts`): ids unique, citations resolve, `_core` stands alone, links stay
  inside the package, the package cites no instance id. The checks are TypeScript run with Node; a project without
  Node ports them or runs them with a pinned runtime, and says which in the manifest.
- **The path check** (`_checks/paths-check.ts`): nothing committed names a path on one machine (§machine-invariant).
- **The module-size check** (`_checks/file-size.ts`): no source file past the hard bound unless the baseline lists
  it, and the baseline only shrinks (the architecture module's rule on module size, if used); a project with files already past it starts its baseline with
  `--update` once, at adoption.
- **The fold-back check** (`_checks/foldback-check.ts`): every session record answers the fold-back question, and a
  closed folder's `foldback.md` carries the author's decision (§reconcile-before-merge).
- **The merge check** (`_checks/merge-gate.ts`): nothing tracked outside the work folder points into it
  (§execution-unreferenced).
- **Host exemptions**: each check takes `--allow <regex>` for file paths it skips (the host's own text that uses `§`
  otherwise, step 0; files that must talk about the work folder). List each exemption in the manifest's adaptations,
  with why.
- **The formatter's ignores**: every agent-context folder, the instance folder (`.<instance>/`: its skill, manifest,
  designs and records) included, and any of the host's where plans, progress notes or session records live
  (§no-auto-format). The one exception is the checks themselves: `.principles/_checks/*.ts` is code, linted and
  type-checked as code (they import `./lib.ts`, so a type-check needs `allowImportingTsExtensions`); everything else
  in the package stays unformatted.
- **The project's own checks**: the finalize step and the drift check (§mechanical-drift-checks) are the project's to
  write; the sketches are in [workflow-files](workflow-files.md#finalize-check) and
  [spec-mirror](spec-mirror.md). Name them, and the host's other repository-wide checks (dependency boundaries,
  exports), in the manifest, so every lane brief lists them (§briefs-name-checks).
- **The gate** (see the testing module, if used) is the build, the test tiers and the repository-wide checks:
  the package's three and the project's own.
- **Verify:** the checks read tracked files (the git index), so stage the package and the instance first
  (`git add`). Each check exits 0 on the fresh setup, and fails on a planted error in a staged file (an unresolved `§`
  citation, a home path, a file forced into the work folder with `git add -f`). A check that can't fail isn't wired.

## 5. Before the first lane

- **The work folder** `.<instance>/execution/<work-id>/`: plan, loop, state, `gaps.md`, `drift/`, `briefs/`
  ([workflow-files](workflow-files.md)).
- **A design names the root nouns** before any lane starts (§name-the-roots): the adoption's own design note, or the
  first task's design.
- **A brief per lane** ([template](workflow-files.md#lane-brief)): model and reasoning effort
  (§model-per-task), files and hot files (the lockfile always, §orchestrator-owns-lockfile), every check including
  the host's (§briefs-name-checks), and that a permission refusal stops the step and comes back (§refusals-to-author).

## 6. While the work runs

- **Gaps and drift in the work folder**, shown at every gate (§gaps-queue, §drift-notes), not scattered across
  reports.
- **A record at every sitting**, compaction or handoff, not only at the close (§permanent-record).
- **No tracked file cites the plan**: a spec, test or comment cites the design or record its content went into
  (§execution-unreferenced).

## 7. Fold back

- What the work taught becomes changes to designs and rules, each through the self-check (§self-check; the
  `principles-self-check` skill), with its delta in the record.
- A change that would help other projects goes back to the package's home: the project is a source the principles
  are tested against ([sources pattern](sources.md)), and its adaptations are the evidence.
- Reconcile and remove the work folder before merge, on the author's approval (§reconcile-before-merge).

## Done when

The manifest names provenance, scope, the folder map and the adaptations; the root `AGENTS.md` points at the package;
the three skills and the instance skill are listed; every check in step 4 exits 0 and fails on its planted error; the
work folder of `adopt-principles` has a plan, a state file and a gaps file, and its session record is written. Lane
briefs are needed only when there is a real first task beyond the adoption.

# Checks

The five checks every project that adopts the package wires into its gate (adopting, step 4,
[../_core/patterns/adopting.md](../_core/patterns/adopting.md)). They are code: TypeScript run with Node, no
dependencies, from the repository root. Each prints every problem with its file and line and exits 1 when there are
any, 0 when there are none (§green-is-exit-code); a bad argument exits 2 with a usage message, never a silent pass.

They read **tracked files** (the git index, so staged files count): what CI will see, not local drafts. Stage the
package and the instance before verifying, and plant each error in a staged file.

| Check | What it rejects | A planted error it must fail on |
|---|---|---|
| [principles-ids.ts](principles-ids.ts) | an unresolved or malformed `§` citation; a duplicate id; a rule, section or pattern without its `{#id}`; a pattern outside its status folder (`bootstrap-patterns/` in the package; `stable-patterns/`, `sketch-patterns/` in the instance); a broken relative link or fragment; a link leaving the package, or a package citation of an instance id (§self-contained); `_core` citing another module | a citation of an id nothing defines, in any tracked file outside the work folder (an archive, §permanent-record) |
| [paths-check.ts](paths-check.ts) | a tracked file naming a path on one machine: a home folder (`Users/<name>`, `home/<name>/...`), a temp folder (`tmp`, `private/tmp`, `var/folders`), a drive path (`D:` and a slash or backslash), also in a file URL or JSON-escaped (§machine-invariant) | a home folder path in a tracked file |
| [file-size.ts](file-size.ts) | a source file over 1500 lines (tests 2500) not in the baseline, or a listed file grown past its baseline count; files over 800 are notes, not failures; `--update` only lowers the baseline (§module-size) | a new 1600-line source file |
| [foldback-check.ts](foldback-check.ts) | a session record whose Learned section doesn't answer the fold-back question (`Fold-back: none` or `Fold-back: proposed (foldback.md)`); a "proposed" with no `foldback.md` in its work folder; a closed folder's `foldback.md` without the author's decision line (§reconcile-before-merge); notes, never fails on, an open folder's `state.md` past 600 words (§state-file) | a record with no `Fold-back:` line |
| [merge-gate.ts](merge-gate.ts) | a tracked file outside the work folder naming a path inside it, or a Markdown file linking into it; files inside the folder may refer to each other (§execution-unreferenced) | a design linking to a plan in the work folder |

## Running them

```
node .principles/_checks/principles-ids.ts --instance .<instance> [--allow '<regex>' ...]
node .principles/_checks/paths-check.ts [--allow '<regex>' ...]
node .principles/_checks/merge-gate.ts --instance .<instance> [--allow '<regex>' ...]
node .principles/_checks/foldback-check.ts --instance .<instance> [--allow '<regex>' ...]
node .principles/_checks/file-size.ts --instance .<instance> [--allow '<regex>' ...] [--update]
```

Node runs TypeScript directly from version 22.18 (or 23.6); with an older Node, run them with `tsx`. Every check
takes the same arguments ([lib.ts](lib.ts)):

- `[root]`: the repository root (one at most), default the working directory.
- `--package <dir>`: where the package is, inside the root; default `.principles`. The id check refuses a folder
  that does not exist (a typo would otherwise check nothing).
- `--instance <dir>`: the project's instance (`.<instance>`), which must exist: its `principles/` defines ids too,
  its `design/`, `sources/` and `_skills/` have their links checked, and its `execution/` is the work folder. Without
  it, the package is checked alone and the work folder is `.execution/` at the root.
- `--allow <regex>` (repeatable): files the check skips, matched against their path. For the host's own files: docs
  that talk about paths in the work folder, a Dockerfile's container paths, a standard cited by section number. List each
  exemption, with why, in the manifest's adaptations.

A project without Node ports them or runs them with a pinned runtime, and says which in its manifest. A project that
changes one keeps the change in the package copy and says so in its manifest's adaptations, like any other edit.

These files are code, the one exception to leaving the package unformatted (§no-auto-format): lint and type-check
them with the project's code. They import `./lib.ts`, so a type-check needs `allowImportingTsExtensions`.

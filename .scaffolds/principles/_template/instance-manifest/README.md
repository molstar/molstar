<!-- The instance manifest: one per project, at `.<instance>/principles/README.md`, created by the bootstrap
(adopting pattern, step 1). The router reads it; the instance skill names it. Replace the <...> placeholders and remove
this comment. -->

# <instance>: the principles here

<One or two sentences: what this repository is and what it builds.>

## Provenance

The package is `<path to the package>`, copied from <where it came from> at <commit or date>. An update from it is a
merge from that base.

## Folder map

| The package says | Here |
|---|---|
| `.principles/` | `<path to the package>` |
| `.design/`, `.sources/`, `.records/` | `.<instance>/design/`, `.<instance>/sources/`, `.<instance>/records/` |
| `.execution/`, `.tmp/` | `.<instance>/execution/`, `.tmp/` (gitignored; the merge check assumes this work folder) |
| `_spec/` | <...> |
| the checks: id, path, merge | <their commands> |

## Modules

The router picks from the "Read before" column. A module of the project's own named like a package module extends
it: read the package's first.

### From the package

| Module | Kind | Read before | Covers |
|---|---|---|---|
| `_core` (in the package) | rule set | any non-trivial work (the `principles` skill loads it) | how the principles work and how work is run |
| `<module>` (in the package) | <kind> | <when to read it> | <what it covers> |

### The project's own

| Module | Kind | Read before | Covers |
|---|---|---|---|
| `<module>` (`.<instance>/principles/<module>/`) | <kind> | <when to read it> | <what it covers> |

## Adaptations

One line each, with why: modules dropped; rules where the host wins; bootstrap patterns taken in as written (the
package's copy serves as the stable pattern) or changed (copied here with a new id, superseding §<bootstrap-id>);
edits to the package's modules; the checks' exemptions (`--allow`), each with why.

- <what was changed>: <why>

## Open items

<Where open items go: the tracker, or a file if the project has none (§followups).>

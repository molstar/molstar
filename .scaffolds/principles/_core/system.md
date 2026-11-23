# How the principles are built {#principles-system}

The principles are a package a project copies and adapts: a stable core, skeletons for new areas, and modules that
each cover one area. This file says how the package is put together; `workflow.md` beside it says how work is run.

## Modules

- **A module is the area, the principle it covers**: architecture, representations, testing, UI. It is a folder
  below `.principles/` (or, for a project's own module, below `.<instance>/principles/`) with a header at the top of
  its `README.md`.
- **Source folders are subfolders by source.** When a module collects material from several sources (one app's UX
  patterns, another's), each source is a subfolder of the module, its source folder. Adding a source adds a folder,
  not a module. A source folder has its own small `README.md` with two fields only: source and status; everything
  else is the module's.
- **Header fields** of a module's `README.md`:
  - **kind**: rule set, pattern collection or layout (below);
  - **status**: stable, directional or draft;
  - **scope**: what it covers, and what it does not;
  - **sources**: where it came from (a prototype, an app, a person), named in words; the source's own folder points
    at the module, not the other way (§self-contained, the [sources pattern](patterns/sources.md));
  - **depends on**: the other modules it cites.

## Kinds

| Kind | Holds | Skeleton |
|---|---|---|
| rule set | `rules.md`: an "in brief", then sections of rules, each rule with an id; its code patterns in `bootstrap-patterns/` (in the package) or `stable-patterns/` and `sketch-patterns/` (in a project's instance; below) | `_template/rule-set/` |
| pattern collection | patterns that are the module itself, by source folder (look and feel, say): one entry per pattern | `_template/pattern-collection/` |
| layout | `tree.md`: how a codebase built on the principles is laid out | `_template/layout/` |

In a text-only pattern collection (look and feel, say), an entry may be a single Markdown file in a source folder
(`<collection>/<source>/<pattern-id>.md`) instead of a folder with implementation files. The file is named
by the pattern's id, which is its anchor, as everywhere else.

## Code patterns

A rule set's code patterns live in the module whose rules they implement. In the package they are
`<module>/bootstrap-patterns/<id>/`; in a project's instance, `<module>/stable-patterns/<id>/` or
`<module>/sketch-patterns/<id>/`. Each folder has a short `README.md` index (id, name, rules). A pattern is a folder named by
its id: a language-neutral `README.md` (the concept and the language-neutral sections) and one implementation file per
language. Today every pattern has `typescript.md`; a codebase in another language adds `python.md`, `rust.md` and so
on beside it, and reads the files for the languages it uses. A pattern spanning modules lives in the module whose rules
it mainly implements and is cited by its id from the others; dropping a module drops its patterns.

**Language-neutral README, one file per language.** A pattern's `README.md` holds the concept: Problem, Shape,
Invariants, Conformance cases, and an "Open" list of the gaps found while writing it, which anyone working on the
pattern can take up. An implementation file (`typescript.md`, and `python.md` and so on when a codebase needs one)
holds only the code and its language-specific notes, and points back at the README.

**Status is the folder**, not a field.

**Bootstrap** (the package's only kind): a starting implementation of the package's rules, to take in where it fits.
It is held to the stable bar below (a wrong bootstrap pattern is a bug in the package), but it is not a rule: a project
may take it, adapt it, or build its own. Taken into a project, it becomes that project's stable pattern:

- **Taken in as written**, it stays where it is: the package's copy serves as the project's stable pattern, and the
  manifest's adaptations say so. Nothing is copied, so no id clashes.
- **Changed**, it is copied into the instance as `.<instance>/principles/<module>/stable-patterns/<new-id>/` with a
  new id; the manifest's adaptations say it supersedes §<bootstrap-id> here, and the router reads the instance copy.

Either way it hardens in the project; changes worth sharing fold back into the package as a merge.

A project's own patterns are stable or sketch:

**Stable**: Stable patterns are what everything else is built on, so they are meant to be correct as written: types check, edge cases are handled, and the conformance cases describe real behavior. They are still dependency-free reference shapes to copy and adapt, not a library. Elided helpers are declared or listed at the end of each file. A bug here is a bug, not a simplification.

**Sketch**: A sketch fixes the *shape*: the concepts, the names, the ports, how the piece connects to the stable patterns, and which conformance cases it will need. It does not claim to be complete or correct in detail: helpers are assumed freely, performance is described rather than proven, and some choices are placeholders. Sketches build on the project's stable patterns and the package's bootstrap patterns and must stay consistent with them (same names, same types, §-references that resolve); they may not quietly redefine a stable pattern's type.

**Several sketches may combine into one stable pattern.** The stable pattern is created in `stable-patterns/` with its
own id, and each sketch it replaces leaves a tombstone in its place that points at the new id
(`{#old-id}` *moved to §<new-id>*), so citations keep resolving and the check reports them for updating.

### How to use them

- Each concept has one owner and one name ("Shared names" below); reuse the name, don't invent a synonym.
- Taking in a bootstrap pattern: copy its shape into the project's code; the code is the project's. The pattern stays
  in the package if taken as written, or moves into the instance with a new id if changed ("Bootstrap" above); either
  way, say so in the instance's manifest.
- Implementing a stable pattern: match its types and behavior; a needed change is drift (§recursive-loop), not a local tweak.
- Implementing a sketch: take its shape as the starting point; where it is wrong or missing, write a drift note with
  the evidence so the sketch can be updated. A sketch never quietly redefines a stable pattern's type.
- Base patterns stay generic: no application roles, vocabularies or policies in `lib/` or a pattern.

### Hardening a sketch

1. **Design** (§design-first): a `.design/<area>/<topic>.md` that settles the open choices for the slice being built, citing the sketch.
2. **Prototype** against real data and the stable patterns, with the sketch's conformance cases made real (see the testing module, if used).
3. **Fold back**: update the sketch with what the prototype taught; when its shape no longer changes under use, move it to the module's `stable-patterns/`. A sketch that the prototype contradicted is rewritten, not patched.

Conventions (TypeScript): strict TS, no `any`, ports are `interface`s named for the capability (not the provider), in-memory implementations ship next to the port.

## Shared names

The patterns build on each other; each concept has one owner and one name. A module lists the names its patterns own
in a "Shared names" table in its pattern folder's `README.md` (`bootstrap-patterns/`, `stable-patterns/` or `sketch-patterns/`; concept, names, owning
pattern); reuse the name, don't invent a synonym, and look the owner up there.

## System folders: the `_` rule

A folder whose name starts with `_` is a **system folder**: its role is defined by the principles system, not by the
project. `_core` (this module: how the principles work and how work is run), `_template` (the skeletons), `_skills` (the skills, one folder each), `_checks` (the checks every project wires into its gate), `_spec`
(the semantic mirror of code, at the root and nested per module) and `.design/_prototypes` (designs of prototypes).
A plain name is the project's own area or module. `_` replaces a leading `.`, because hidden folders are skipped by
some tools and globs.

- **The package is self-contained.** {#self-contained} Nothing in `.principles/` links to or names a concrete file
  outside it: not a project's designs, records, sources, prototypes or code. It may name the folders whose role it
  defines (`.design/`, `.records/`, `.sources/`, `_spec/`) and the companion checks that enforce it. Everything else
  points at the principles, never the reverse, so the package means the same when copied into another project. The
  check rejects a relative link from `.principles/` that leaves it.
- **`_core` is stable, not immutable.** Every project starts with it and rarely changes it; a project may still change
  it when it needs to.
- **`_template` holds skeletons only.** The starter content lives in the modules, which are iterated on in place.

## Ids

- **Every rule and every pattern has a stable kebab-case id**, unique across the package, cited as the section sign
  followed by the id, `§<id>`: §design-first, §machine-invariant. An id means nothing about order; reading order is
  the order in the file. The section sign is used for nothing else, so a checker finds every citation without false
  hits.
- **In a file, the id is an anchor** after the rule's bold lead, `**Requester authority** {#requester-authority}`,
  or on its heading, `# Transactions {#transactions}`. Sections of a rule set have ids too, and so does each rule
  set as a whole (`{#architecture-rules}`), for citing a module's rules at once.
- **Removing a rule leaves a tombstone**, one line in place of the rule: `{#old-id}` *removed (2026-11): reason*, or
  `{#old-id}` *moved to §<new-id>*. Old citations keep resolving. A tombstone that only redirects is reported by the
  check, so its citations get updated; it can go once nothing cites it. An id is never reused.
- **No numbers**: the principles once used numbered ids (a letter, a section, an item); every citation now uses the
  id, and the check rejects an old number with a hint to cite the id.
- **The check** (`_checks/principles-ids.ts`, in the gate; [the checks](../_checks/README.md)): ids are unique; every
  citation in a tracked file resolves to an id; redirect-only tombstones are reported; relative links resolve, and
  none leaves the package.

## Skills: three in the package, one in the instance, and they link

- **`principles`** (`_skills/principles/SKILL.md`): how work is run; loaded at the start of non-trivial work.
- **`principles-router`** (`_skills/principles-router/SKILL.md`): which modules to read, from the manifest's "Read
  before" column.
- **`principles-self-check`** (`_skills/principles-self-check/SKILL.md`): recommends and runs the four-step check
  (§self-check) when work changes the principles or designs: a fold-back, finalizing, a tension that won't dissolve.
- **The instance skill** (one per project, named after the repository; `_template/instance-skill/SKILL.md`): what the
  repository is and where its instance keeps things; loaded first, it points at the instance's manifest. It lives in
  the instance, not the package.
- Skills link into the modules and copy nothing: no rule lists, no summaries. **The "in brief" lives in the module**
  (the top of `rules.md`, or the module README), and the "when to read" line lives in the manifest, so adding or
  dropping a module changes the manifest, not a skill.
- **`_skills/` is a system folder** holding one folder per skill with only its `SKILL.md`, whose links point into the
  modules. `.agents/skills/<name>` links to `.principles/_skills/<name>`, the same way `.claude/skills` links to
  `.agents/skills`, so a skill folder never exposes the whole package.

## The manifest

A project keeps its manifest in its instance (`.<instance>/principles/README.md`, from
`_template/instance-manifest/`): the package's provenance, the folder map, which modules it uses (the package's and
its own, each with its "Read before" line), and one line each on any adaptation (what was changed and why). The
router reads it; the instance skill names it. The module list in `.principles/README.md` is the package's own
catalogue, the starting point a manifest is made from.

## Adapting the package to a project

The full procedure, from the questions to ask first to the first lane and the fold-back, is the [adopting pattern](patterns/adopting.md); [bootstrap.md](../bootstrap.md) is the entry point that sends an agent there.

1. Copy `.principles/`; create the instance, `.<instance>/` with its manifest (adopting, step 1); create the skill
   links and the instance skill (adopting, step 3). Run the checks in `_checks/` in your gate (adopting, step 4).
2. Drop the modules the project does not need and remove their rows from the manifest. `_core` cites no other
   module, so it stands alone, and every citation of a dropped module's id, in a kept module or in your own text,
   fails with its file and line, so you see what to reword (e.g. dropping `ui` names the lines in `architecture`,
   `representations`, `testing` and `layout` that cite its ids).
3. Add the project's own modules from `_template`, in the instance (`.<instance>/principles/<area>/`): copy the
   skeleton for the kind, fill in the header, write the rules or patterns, and add its row (with "Read before") to
   the manifest.
4. Edit any module in place, `_core` included. Git keeps the history; the manifest line says what was changed and why.

**Updating from the package later** is a merge, like any other: modules are separate folders, so a project's edits
stay local to the modules it changed.

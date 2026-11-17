# Configuration and extensibility

## Intent

The author, 2026-11-17: "basically the rule is: what we can make configurable, we should", and "i am also missing notes
on "Extensibility" in the storybook".

Seq\* is a library that applications compose: a structure viewer's sequence panel, a genome browser's tracks, a design
tool, a figure in a paper. Each wants different behavior, looks and data. The defaults must be good; nothing about them
should be fixed when it can be a setting or an extension.

## Precedents

- The previous implementation: a viewer spec shaped like Mol\*'s plugin spec (a registry, behaviors with params, a
  config, render options); extensions with namespaces, registry kinds, a fixed registration order, conflict checks,
  mandatory kernel fallbacks and version upgrades ([sources/seqstar](../sources/seqstar/README.md)). Rich, and heavy to
  learn.
- meso\*: descriptors as data, gathered by the caller; nothing registers itself on import ([sources/meso](../sources/meso/README.md)).
- `@molstar/ui`: commands, keymaps and settings as DOM-free models.

## Decision

### Configuration: what can be configured, is

- **The rule is the direction; this version builds its first step** (the author's rule, "what we can make configurable,
  we should"; the review's cuts accepted, "rest as recommended"):
  - **in this version**: the spec holds what is shown and how it looks as data; the viewer's configuration is one typed
    set of config items given at composition (`ParamDefinition`, as Mol\*'s `PluginConfigItem`s), holding the foundations'
    preset and numbers, the grammar's preset and bindings, behaviors' params, limits and defaults;
  - **later**: user settings layered over it, the effective value of a setting with where it came from, and a generated
    settings form (`@molstar/ui`'s params forms read `ParamDefinition`, so the form needs no second description).
- **What is configurable** (the first list; anything else found fixed is a gap to fix, not a rule): the foundations'
  numbers in presets (compact, comfortable, figure, print); the grammar's bindings and presets; each behavior and its
  params; density thresholds and letter sizes; the data tokens (`--seq-*`) and every palette and residue scheme;
  animation durations; chrome parts shown; limits (path length, issue counts, request concurrency); the codon tables and
  numbering schemes available.
- **Configurable does not mean inconsistent**: one configuration applies to every view of a viewer, and the fixed rules
  of a configured area (the grammar's, the foundations') hold whatever the values ([ux-foundations.md](ux-foundations.md)).

### Extensibility: every registry is open

- **Extension points**, each a typed contract with the core's or the viewer's built-ins as its first implementations:

  | Point | What an extension adds |
  |---|---|
  | node kinds | new kinds in the spec, in any part of the tree (a data kind, a computed node, a view element), each with its params described like the built-ins' and a fallback or an "unavailable" rule for readers without it |
  | custom fields | namespaced fields on any node (`custom: { "<ns>": {…} }`) and extra columns on any table: kept through parse, patch, diff and print untouched by readers that do not own them, validated by the extension that does, readable by themes, filters and queries as fields |
  | formats | readers (a data format onto the one protocol: FASTA, GFF, bigWig, PAF, a SIFTS table, …) and writers (a selection as FASTA, a figure as SVG or PDF, a view as a session for another tool) |
  | layer types | a way to draw (its params, its layout and draw-list output, its hit-testing, its fallback) |
  | ports | the three ports' adapters ([ports.md](ports.md)): assets, text metrics, structure |
  | themes, palettes, residue schemes | how a channel is computed |
  | map runtimes | a way to compute a map of a given kind (align, SIFTS lookup, gene model to codons) ([maps.md](maps.md)) |
  | numbering schemes, codon tables | labels and translation |
  | query functions | generators, modifiers, filters in the query language |
  | commands, behaviors, bindings | what the viewer can do and what gestures do |
  | backends | where draw lists go ([rendering.md](rendering.md)) |

- **Internal until an outside user**: layouts (tracks, text, rings, dot plots) and chrome parts (track header actions,
  context menu items, hover card sections, panels) are written behind internal interfaces, opened as extension points
  when an application needs one; definitions and uses wait likewise ([composition.md](composition.md)).
- **An extension is registry entries**, as a Mol\* plugin's ([viewer-composition.md](viewer-composition.md)): a
  namespaced record of providers for each point, listed in the viewer's spec (one composition root) or registered by a
  behavior at run time with its undo; nothing registers itself on import, and there is no global registry.
- **Built-ins use the same points**: every layer type, theme and runtime the core ships is written as an extension
  would be. Why: a point is real only if its first user is the library itself.
- **In the spec, extension kinds are namespaced** (`ext:kind`) and declare a fallback to a built-in; a reader without the
  extension draws the fallback or makes the node unavailable, with an issue, and the rest of the spec opens (EC-235).
- **Custom fields never break a reader**: a field under a namespace the reader does not know is carried through
  unchanged (a patch, a diff and a print keep it), and a host's own fields go under its namespace, never as new
  top-level params.
- **Registration is atomic, counted and undoable**, as the plugin's: every conflict (two providers under one name, a
  binding claimed twice) found before anything changes and reported as one error; never resolved by order.

## Rejected alternatives

- A fixed registration order with priorities (the previous design): results depended on order; here conflicts are
  errors and composition is explicit.
- Registration on import: an import's side effect decides what a viewer can do, invisible at its composition root.

## Consequences

- The Storybook's Extensibility page shows a small extension end to end (a layer type, a theme, a command and a
  binding) and a composed viewer's configuration.
- The docs and the agent skill list every setting and extension point from their descriptions.

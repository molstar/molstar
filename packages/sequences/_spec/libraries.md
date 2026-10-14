# Extensions and application-level libraries

The packages built on the core, each through what any extension may use: registry entries, behaviors, commands, ports.
The core imports none of them.

## Interactions

`@molstar/seq-viewer/interactions` holds one set of opinions about interactions under its namespace; another vocabulary
would be another extension, and the kernel would not change.

- **The vocabulary.** An interaction is a row of a kernel `pairs` annotation (by convention, semantic type
  `interactions:interaction`); the extension reads only its category, one of `contact`, `hbond`, `salt_bridge`,
  `pi_stacking`, `cation_pi`, `hydrophobic`, `halogen`, `metal`, `water_bridge`, `covalent`, `other` (any category it
  does not know, and a row with none). By convention, `direction` is an `hbond`'s donor side and `distance` (Å) and
  `weight` are number fields; the extension does not read them.
- **Queries** over the kernel's pair index (`kindOf`, `contactsOf`, `atomsContacting`): a row's kind; the positions
  contacting a ligand unit, with their kinds; a ligand's atoms contacting a residue.
- **Fingerprints**: ligands or poses (rows) × positions in a shared coordinate system, a set of kinds per cell. The
  provider `interactions:fingerprint` returns a kernel `matrix` (the ligands as states, the number of kinds as value)
  with the fingerprint as a runtime property, so the kernel places it as any matrix. Contacts of residues with no path
  into the system are a warning; non-polymer partners are left out.
- **The layer type** `interactions:fingerprint` draws a row per ligand, a stripe per kind present; its fallback is
  `matrix`. The theme `interactions:kind` colors and marks by kind. A host with only the format half draws the kernel
  fallbacks, with the builder's kernel fallbacks drawing `other` too.

## Links

`@molstar/seq-core/links` turns links in an application's markdown into invocations. The core has no markdown.

- A link whose target starts with `!` is Mol\*'s syntax, `!key=value&key=value`, URI-decoded; one link may carry Mol\*
  keys and Seq\* keys, each side running the keys it knows.
- **Default keys**: `seq-highlight` (on enter, `seqv:highlight` on the link's own channel, closing by containment;
  cleared on leave), `seq-select` (also selects on click), `seq-focus` (also focuses on click), `seq-snapshot` and
  `seq-view` (on click). Targets are comma-separated refs with range args (`start`, `stop`, or `beg`, `end`,
  `numbering`, and `items`). A SeqQL `query` target is not built.
- The runner is given a function that runs an invocation (`ctx.commands.run(inv, "host:links")`); it remembers the links
  entered, so a link removed while hovered leaves nothing lit.
- An application adds its own keys, each mapped to its registered command.
- Adapters: `/remark` is a rehype plugin (`rehypeSeqLinks`, with `seqLinkSchema` for the sanitizer) that keeps `!` links
  and marks each with an id from its href and its occurrence (`seq-link-<k>:<href>`); a `:` in a target is written
  `%3A`.

## Search

`@molstar/seq-core/algorithms/search` is a search library, composed by applications. What is searched, what a hit means
and how hits rank are the application's; the viewer shows whatever it is given (`ctx.results`), and the format has no
search query.

- **Matchers** find a pattern in code arrays: exact, regular expression, letter classes (a motif over a class table),
  and approximate (Myers' bit-vector infix search, at most `maxEdits` edits).
- Hits arrive through a callback as `(start, stop, score)`, `score` the edit count; nothing is allocated per hit (the
  regex matcher aside). Hits come in order of their ends and may overlap.
- **The window**: hits start inside `[from, to)`, and end by `to` unless the sequence is cyclic: there a hit starting in
  the window may run on past `to`, across the sequence's end when it gets there (`start > stop`), at most one lap long.
- **Approximate ends**: `local_best_last` (default) and `local_best_first` give one hit per site, the last or first end
  of a run of equal minimal cost; `all` gives every end within `maxEdits`. A hit's start is that of the shortest stretch
  of its cost ending there.
- SeqQL's `search` is a thin, explicit use of these matchers.

## Translate

`@molstar/seq-core/translate` maps between coordinate systems, a host's own and Seq\*'s, without the viewer. It depends
on the model only and keeps the registry the stateless model does not; the viewer's `ctx.systems` is one `Translator`.

- **Registration.** `register(descriptor, maps)` gives a key its system (an external one, or a Seq\* one it names);
  registering again replaces its maps and keeps its id; `unregister` drops the maps, the key keeps its id. `addMap`
  joins any two systems.
- **Composite systems**: a descriptor with members is a complex, with no index of its own.
- **Routes**: the fewest maps, at most 4, at most one with span 3; every shortest route is taken and the results unioned
  (the viewer's projection rule, not the model's ambiguity rule). Routes are cached until any registration changes.
- **Extra indices**: `index(system, kind)` builds a host's `IndexKind` (a name and a `build`) over one system on first
  use and holds it until the system is registered again or unregistered. `IndexKinds` has the common ones: a structure's
  residue to atom ranges, a complex's chains to units, a numbering's label to index, a SIFTS table sorted.
- `translate` (runs to runs) and `forEachAt` (one coordinate, every image: allocation-free, hover rate).
- **Map builders** from common identities: `byLabels` (equal labels pair), `byNumbering` (auth numbers with insertion
  codes), `bySegments` (SIFTS-style triples, 0-based, half-open).

## The Mol\* adapter

`@molstar/seq-molstar-poc`, a proof of concept, is the only Seq\* package that imports Mol\*. It reads structures
headless with Mol\*'s parsers (no canvas, no plugin UI), so its format half runs in a page, in Node and on a server.

- **Formats and services**: `mmcif`, `bcif`, `pdb`, `sdf`, `mol2`, `afdb_pae`; `pdbe` (`structure`, `sifts`) and `afdb`
  (`structure`, `pae`). A file is parsed once; a host that parsed it for its 3D view shares the parse. Sequence and
  alignment files (FASTA, A3M, Stockholm, Clustal, Newick) are the kernel's (`seqspec:sources#formats-and-services`),
  not the adapter's.
- **`molstar:structure`**, the complex provider: units are Mol\*'s chain instances (one per `label_asym_id` and
  operator); a polymer entity is a sequence over its full `entity_poly_seq`, with numberings `label`, `auth` (author
  numbers with insertion codes) and, from SIFTS, `uniprot:<acc>`; a non-polymer entity is a ligand over its heavy atoms
  with the `atom` numbering. Options choose the model, assembly, what is included, the chains, unit names and
  numberings. Every unit carries its `mvs` fields. Unit names must be unique (`molstar.unit_name.collision`). An
  insertion code that is not letters is read as none, once reported (`molstar.ins_code.invalid`, a warning).
- **`molstar:plddt`**: per-residue pLDDT from B-factors, colored by AlphaFold DB's bands by default.
  **`molstar:sifts_coverage`**: the listed structures' chains on one UniProt sequence, one `ranges` row per resolved
  stretch with a `structure` field (`instance: true`: rows on each chain's own unit); `molstar:sifts_depth` beside it.
  **`molstar:pae`**: a model's predicted aligned error (`seqspec:sources#formats-and-services`) as a `grid`, a row per
  pair of the complex's polymer units (a residue a token, in the complex's order), side a the scored residue and side b
  the one aligned on (a dot plot draws it as AlphaFold DB does), values in Å; its size must be the units' residues
  together (`svs.source.format` otherwise). Its default color is AlphaFold DB's green, 0 to 30 Å (`PaeScale`).
- **The link** (`StructureLink`): each chain instance but waters is an external system mapped onto the unit with the
  same `label_asym_id` and `instance_id` (else the first with that `label_asym_id`; polymers by `label_seq_id`, others
  by heavy-atom order), so a structure in 3D and its sequences light each other through `ctx.systems`, every copy
  included. Mol\* → Seq\* runs at hover rate without allocating; Seq\* → Mol\* loci are made at edge rate.
- **`linkViews`** links a viewer's hover and selection with a 3D view both ways. What it sends is composed in its
  options (`LinkViewsOptions`, defaults `DefaultLinkViewsOptions`), as a viewer's reactions are in its spec: a Seq\* set
  within `maxItems` items and `maxLength` positions (read from its size before it resolves) goes residue by residue, a
  larger one `over`: the whole chains it touches (the default) or nothing. Defaults: hover 2,000 items and 20,000
  positions, selection 20,000 and 200,000.
- **SeqQL**: `molstar.select` (`molql`) runs a PyMOL (or VMD, Jmol, MolScript) selection on the linked structure and
  returns its residues on their units; `(input structure)` reads Mol\*'s selection passed by the host.
- **MVS**: the structure, its colors (as MVS annotation tables), sticks and a selection's label go to the 3D view as an
  MVS state; highlight, select and focus go as Mol\* loci through the plugin glue.
- **The viewer half** (`/viewer`): a track's colors as drawn, host actions ("Show in 3D"), and linking a viewer to a
  Mol\* plugin; `/plugin` mounts Mol\* as a canvas only.

## Nightingale

`@molstar/seq-nightingale` turns a resolved view into EBI's Nightingale web components and their data. It is the
format's check, not its limit: it shows that a document says what to show and not how our engine draws it. A type
Nightingale cannot draw stays in the format; Nightingale draws its fallback and reports the gap, with the format's
codes. It reads the same resolution as the viewer (resolved with its own supported types), registers no viewer entries
and special-cases no extension; `materialize` is platform-free and only `/mount` imports Nightingale. Of streamed data
it fetches only `values` placed directly, once over the whole axis at the level whose bins are about a pixel; other
streamed layers, and values through a path, are not drawn (`svs.stream.fallback`). It has no linked viewports, and no
two-axis tracks: a `dot_plot` is left out and reported (`svs.render.omitted`, its fallback). It targets small systems:
its tracks take an element or a datum per position, so a section past `NIGHTINGALE_MAX_POSITIONS` (a million) is left
out and reported (`svs.render.omitted`), and the view's other sections render. Highlight and hover go through the shared
projection.

What the check showed about the format: shapes, encodings and heights are hints a viewer may drop without losing
meaning; per-position color is its own type; pairs carry categories and endpoints in invariant coordinates; highlights
are loci, converted at the edge; themes and glyphs carry kernel fallbacks.

### Known losses

What Nightingale drops of every layer of a type, unreported, since a hint is a viewer's to drop: a `values` layer's
`bars` drawn as its line; a `logo`'s `scale`; a `ranges` layer's `display` `union` and `underline` (drawn as `lanes`,
Nightingale's own rows of rectangles) and its `collapse`; a `pairs` layer's `display`, `arcs` and `labels`; the size
channels (`thickness`, `radius`, `width`, `height`), over a field too; `density`; `residue_form`; a section highlight's
`style` (every entry lit alike); `follow` (no linked viewports); a `matrix`'s default palette (the heatmap's own
colors). Everything else it does not draw is reported: a chosen theme it has no slot for (`svs.render.theme`), a
direction (`svs.render.direction`), an `opacity` its element lacks, a ruler of another system (`svs.render.fallback`).

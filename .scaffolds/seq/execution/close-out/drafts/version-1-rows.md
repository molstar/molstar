# Version 1's rows no milestone took: draft decisions (for the author)

Listed on 2026-11-23 from [design/features.md](../../../design/features.md):
- every `planned` row;
- every `deferred` row whose trigger is the Genomes milestone or this close-out;
- every `building` row's part marked "not".

Rows deferred to an outside trigger are not listed: a collaborative host, a React host, the first writers extension,
the typed configuration's users. They stay deferred and need no decision now. The edge cases' "outcome not recorded"
items are not listed either: they are Finalize's reconcile (`ec-ids.ts --worklist`).

Each row gets **build** (a close-out lane), **defer** (past version 1, with its trigger) or **cut**. The column is
the recommendation. Why: version 1 finishes what its views already show and leaves new kinds of work to their trigger.

## 1. Out of `deferred`, already built (the checklist's own item)

Two assemblies on one axis, layers across sections, variants by sample, gene hierarchies, the name index. These are
built by the fifth slice. The rows move once the evidence is linked; there is nothing to decide.

## 2. `planned`, never built

| Row | What | Rec. | Why |
|---|---|---|---|
| Search and results (182) | matchers, stepping, saved as sets | **build**: exact, degenerate, regex; stepping; saved as a set. Defer approximate | three other rows wait on it (Focus and reveal, Highlight channels, Security bounds' search bounds); p53's search story |
| Patches on streamed things (95) | a bookmark or ref naming a dynamic node; waits while not held | **build** | p53's transcripts already stream: without it, a bookmark to one is unsafe |
| Keep this region (96) | streamed data turned into spec content on purpose | **defer**: a host asks to save streamed data in a spec | a new command with no view needing it |

## 3. `deferred` to the Genomes milestone or this close-out, not built

| Row | What | Rec. | Why |
|---|---|---|---|
| Values for genomics (130) | points, split by sample, log scale, reference line (the Manhattan family) | **build**: points, log scale, reference lines. Defer split, with the Ranges split | channels and marks on the existing values layer; split is a layout question |
| Keyed tables and joins (72) | rows keyed by id joined to items (EC-106, EC-224) | **defer**: a dataset with a table keyed by id | a language change this late is the riskiest kind; nothing shown needs it |

## 4. `building` rows, the parts not built

**Build** (small, or a shown view relies on it):

| Row | The part | Note |
|---|---|---|
| Selection (178) | snapping to codons over nucleotides (EC-149) | p53's protein on its gene |
| Keyboard cursor (179), Lazy rows (142) | the cursor across block gaps (EC-154) | one item, two rows |
| Settle and open domains (113) | legends state their domain | fits the legend lane (questions A) |
| Validation (37) | the Antibody view's invalid labels shown | a story's missing half |
| Fallbacks (143) | built-in types' declared fallbacks | with the version-aware validate |
| Encodings (39) | the views' fixtures round-tripped through MessagePack | a test over what exists |
| Security bounds (207) | search bounds | with Search |
| Tracks layout (159), Trees (70), Collapsed clades (141), Track and group editing (185) | groups and summaries, collapse (EC-103) | one design closes four rows; a strong lane |
| Ports with default adapters (206) | "a structure viewer is a library, not a port" | not a gap: reword the row |

**Questions** (could go either way):

| Row | The part | Rec. | Why |
|---|---|---|---|
| Map kinds (50) | `pairwise` and `alignment` | **defer**, with Hover across linked sections' orthologs view (177): the orthologs dataset | both wait on one view no slice built; together they are a milestone, not a lane |
| Ranges (128) | split and collapse | **build** collapse with the track groups; **defer** split with values' split | collapse is the groups design; split is a layout of its own |
| Formats in the core (74) | A3M and Clustal | **defer**: a dataset in either | small readers, but no example uses them |
| Levels of detail (109) | bigWig's own zoom levels | **defer**: a measurement showing bin-width levels too slow | the gate's measurements decide |
| Canvas2D backend (195) | glyph atlas, band copy | **defer**: the same | performance work without a failing bound |

**Defer** (each to its trigger):

| Row | The part | Trigger |
|---|---|---|
| Language as data (34) | docs and the agent skill | agent tools (deferred, the author, 2026-11-17) |
| Patch streams (97) | the agent's editing story | agent tools |
| Archives (40) | every dataset as an archive | an archive writer beyond zip |
| Map sourcing (51) | computed maps | the compute port |
| Direction on every run (53) | reversed pairs rows | a dataset with them |
| Ligands (68) | AlphaFill's transplants, chips in text (EC-99) | an AlphaFill example |
| Services (76) | live | a host asks; recorded fixtures stay |
| Structure files (77) | BinaryCIF and Mol\*'s readers | the structure link on the new mol viewer |
| Session state (99) | open panels | the host owns panels ("engine, not the game") |
| Change kinds (100) | the command palette | a host's palette |
| Self-describing sources (106) | assembly | two assemblies from one source |
| Grid and dot plots (135) | two grids on one block pair (EC-59) | a dataset |
| Query language (183) | extension forms and properties | an extension needing them |
| User-defined refs (184) | refs in data params, a section's `of` and `follow`; the Antibody view's refs | an example needing them |
| Configurable layout (208), Typed configuration (209) | presets, `ParamDefinition`s | settings in layers |
| Host data (213) | sections over a source's complex | in the TODO already |
| Node kinds (214), Registry entries (215) | beyond layer types; the other registries | an outside extension |
| Mol\* structure link (222) | the adapter for the new mol viewer | that viewer |
| Addressing (223) | rects | a link to a rectangle |
| Markdown links (224) | one href driving Mol\* too | the adapter above |
| More genome formats (75) | 2bit, BigBed, PAF (designed only) | each with its dataset (its row's own trigger) |

## What the builds add up to

Lanes, beside the four from the design questions:

| Lane | Rows | Model |
|---|---|---|
| Search | 182, 180, 181, 207 | strong |
| Track groups | 159, 70, 141, 185, 128 collapse | strong |
| Streamed patches | 95 | strong |
| Values marks | 130: points, log scale, reference lines | Sonnet |
| Small items | 178, 179/142, 37, 143 with the version-aware validate, 39 | Sonnet |

113 joins the legend lane. 206 is a row's wording, the orchestrator's.

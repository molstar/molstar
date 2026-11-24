# Search: matchers and commands in the library, the box the application's

<!-- Agreed by the author, 2026-11-24 ("Agree all"): close-out CL5; features 182, and the search parts of 180, 181,
     207. Sections marked **Provisional** say what would settle them. -->

## Intent

A reader looks for a motif: a restriction site on a plasmid, a primer on a 1 Mb genome, a glycosylation motif in a
protein, a pattern across 10,000 rows of an alignment. They see every hit lit in the view, a count ("2 of 14"), step
through the hits with keys, and keep the hits as a set they can come back to. The previous implementation proved the
parts (EC-161); this design rebuilds them on version 1's types (Location, refs, commands, highlight channels).

**This is application-level work on a base built for it** (§application-vs-base; the author, 2026-11-23: "this should
likely be "app port" as well, where the app provides the UI etc"; then "yes, update it and take the rest as
recommended"). The library gains:

- in the core, **matchers** (exact, degenerate, regex) over a system's stored letters, bounded, giving hits and their
  Location; a **matcher registry entry** where an application adds others;
- in the viewer, a **search model** (the query, the hits, the current hit, the status) and **commands** (run, next,
  previous, go, clear, save), a **highlight channel of its own**, keys and the palette (through the commands);
- no chrome: the search box and the results list are the application's. The Storybook composition gets a small box
  built only on the public surface, so a story works and the surface is proved to be enough.

What the library does not ship: approximate search, a server's BLAST, a worker. Each is a matcher an application
registers; the network stays out of the library.

## Precedents

- **The previous implementation** (`packages/sequences`, read only): taken, the kernels' shape (sync, over a byte
  array, hits to a visitor with nothing allocated a hit, in the order of their ends; overlapping hits; a cyclic search
  reading at most one lap past the end: `core/src/algorithms/search/matchers.ts`), Horspool for exact patterns and
  shift-and for letter classes, the results' stepping ("2 / 14", the current again clears it, a set kept: EC-161), the
  capped answer said as `n+`, a box that waits 200 ms and aborts the previous search (`viewer/src/search-box.ts`). Not
  taken: its regex matcher (a JavaScript `RegExp`, a backtracking engine over untrusted text), its approximate matcher
  in the library (Myers: now an application's), its Search part in the viewer's chrome, results that only an
  application could fill (`ctx.results`: here the library runs the matchers it ships).
- **SnapGene**: taken, IUPAC degenerate matching by default, both strands, hits across a circle's origin, "N of M"
  with next and previous. Not taken: a protein pattern searched in the translations of DNA (an application's matcher
  can), mismatches (approximate: the application's).
- **IGV** (Find Motif): taken, a regex or IUPAC pattern, both strands told apart, the whole sequence searched (not
  only the view). Not taken: results as two new tracks (here a channel; saving makes a ref, and a ref can be a layer's
  data).
- **JBrowse 2** (the sequence search adapter): taken, regex, a strand choice, case-insensitive. Not taken: a search
  per region in view, whose answer changes as one pans and has no count to step through.
- **UGENE** (Find Pattern): taken, a cap on results and on a regex's hit length, a search within a region (here
  `within`), results kept (as a ref here, annotations there). Not taken: substitutions and indels (approximate), search
  in translation.
- **Jalview** (Find): taken, an alignment's rows searched ungapped, a hit spanning the gap columns between its
  residues, find next and find all. Not taken: the find selecting its hit (here the current hit has a channel, the
  selection stays the reader's), the find over sequence ids and descriptions (names are `findNames` and the locus box,
  EC-228: an application's box may try a locus or a name first, then a pattern).
- **`include` patterns** (`core/src/lang/include-pattern.ts`): the precedent for untrusted patterns: a safe subset
  matched by a table, never a backtracking search.

What a reader can do there and not here: mismatches and indels (UGENE, SnapGene), a protein pattern in DNA's frames
(SnapGene, UGENE): both an application's matcher through the registry, not built in the library; several named result
sets at once (the previous implementation): not built (one search a viewer, Provisional below).

## Decision

### 1. Two layers: kernels and matchers

- **A kernel** is the hot loop: synchronous, over a byte array window, hits to a visitor `(start, end) => void`, nothing
  allocated a hit, no alphabet knowledge (a code is a byte), as the previous implementation's (`Matcher.find`). It is
  the one callback in the surface, on purpose: a hit a call with no object made (§zero-alloc-hot-loops; the exception
  to §return-observables said in its comment).
- **A matcher** is the registry's port, the surface an application implements:

  ```ts
  interface SearchMatcher {
    readonly name: string;                 // 'exact' | 'degenerate' | 'regex' | '<ns>:<name>'
    readonly title: string;                // 'Exact', 'IUPAC', 'Regular expression'
    readonly alphabets?: readonly SequenceAlphabet[]; // the systems it reads; absent: any
    /** Sync, bounded, no letters read: the pattern checked and compiled, or its issues (Decision 2). */
    compile(pattern: string, options: SearchOptions, bounds: SearchBounds): Compiled | readonly Issue[];
    /** Async and streamed: batches of hits in order; unsubscribing cancels. */
    run(compiled: Compiled, targets: readonly SearchTarget[], ctx: SearchRunContext): Observable<SearchBatch>;
  }
  interface SearchTarget {
    readonly system: string;               // its address
    readonly alphabet: SequenceAlphabet;
    readonly length: number;
    readonly cyclic: boolean;
    readonly within?: Int32Array;          // ranges to search in, from `within` (Decision 5)
    readonly letters: LetterReader;        // read(start, end, signal): Promise<Uint8Array>, by region
  }
  interface SearchBatch {
    readonly hits: SearchHits;             // columnar: target index, start, end, strand
    readonly scanned: number;              // letters read so far, for a progress said once
    readonly issues: readonly Issue[];     // a target skipped (alphabet, too long, failed read)
  }
  ```

  `Observable` is `@molstar/core/reactive`'s (§api-surfaces: many values over time is a stream). Cancellation is the
  subscription's end and the `AbortSignal` in `ctx` (passed to every letter read); a cancelled run ends with no
  error (§cancellation-codes). A server's matcher ignores `letters` and sends what it knows (an accession, the
  sequence's name): the network is the matcher's, never the library's.
- **`kernelMatcher(kernelOf)`** (core) makes a matcher of kernels: it reads each target's letters in chunks of
  `search.chunkLetters` (2^20) through the reader, overlapping by the longest hit less one, runs the kernel on each
  strand, and emits one batch a chunk, sorted by start (a chunk's hits are few; a hit belongs to the chunk it starts
  in). The three built-ins are `kernelMatcher`s. Why: an application's in-process matcher (a PROSITE reader, a
  mismatch kernel) writes a kernel and nothing else; one running in a worker or a server writes `run`.
- **Batches are a prefix of the final order**: targets in order (Decision 4), chunks in order, each sorted, so the hits so far
  never reorder. Why: the count and stepping work while a long search runs.

### 2. What the built-in matchers take

All three are **case-insensitive** (a genome's lowercase is soft-masking, never meaning), read a system's **stored
letters** (the parent's letter of a modified residue, core › Molecules; never the letters a theme shows), and treat
U as T on nucleotides.

| Matcher | A pattern letter means | A sequence letter matches when |
|---|---|---|
| `exact` | itself | it is the same letter |
| `degenerate` (default) | its IUPAC set: nucleotides `R Y K M S W B D H V N`; proteins `X` any residue, `B` D or N, `Z` E or Q, `J` I or L | its own set is inside the pattern letter's set (`N` in the sequence matches only `N`; `R` matches `R`, `D`, `N`, never `A`) |
| `regex` | a letter is itself; classes say sets (Decision 3) | as the automaton reads it |

- **Degenerate is subset matching** (decided here): a sequence's ambiguity letter is a hit only where every base it may
  be is allowed. Why: a site reported on an `N` run is not a site; overlap matching floods a draft genome's gaps with
  hits. Rejected: overlap matching (below).
- **Gaps are not letters**: an alignment row is searched over its residues (`-`, `.` skipped, Jalview's), a hit's
  extent the columns from its first residue to its last, gap columns between them included.
- **The pattern's alphabet**: a target whose alphabet cannot hold the pattern (an `E` in a DNA pattern, a `U` in a
  protein pattern under `degenerate`) is skipped with `seq.search.alphabet` (a warning naming the system), the others
  searched. The system's declared or detected alphabet decides (EC-37), never a guess from the pattern.
- **Kernels**: `exact` Horspool (sublinear, any length to the bound); `degenerate` shift-and over a table of 256
  masks, multiword past 32 letters; `regex` Decision 3 (below). Each pattern position costs a bit: work is letters × ⌈positions / 32⌉
  word steps, the unit of the step budget (Decision 6).

### 3. A regular expression's safety: an automaton, never backtracking

- **The language** (a subset of POSIX ERE, written in its own module like `include-pattern.ts`): letters, `.`, classes
  `[…]` and `[^…]`, groups `(…)` and `(?:…)` (the same: nothing is captured), alternation `|`, `?`, `*`, `+`, counted
  `{m}`, `{m,}`, `{m,n}`, anchors `^` and `$` at a system's ends (a cyclic system has none: `seq.search.pattern.anchor`),
  `\` before a character to take it literally.
- **Refused with a reason** (`seq.search.pattern.unsupported`, error, its position in the pattern and the hint "a search
  runs in linear time: no backreferences or lookaround"): backreferences, lookaround, lazy and possessive quantifiers
  (meaningless under Decision 4's longest-per-start rule), `\b`, `\d`, `\w`, Unicode classes, flags.
- **Compiled to a Glushkov automaton** (one state a letter position, no ε-moves), simulated bit-parallel (Navarro and
  Raffinot's tables), counted repeats expanded first. Bounds checked at compile, before any table is allocated: the
  pattern's characters (`maxRegexChars`, 256), the expanded positions (`maxRegexPositions`, 128: four words), nesting
  (64, the queries' depth, EC-76). Why this engine: its cost is the letters times the words, whatever the pattern
  (`(A|AA)*B` over a million `A`s costs what `AB` does); a `RegExp` with the `v` flag is still backtracking, and V8's
  linear engine is a flag no browser sets.
- **Hits**: one a start, the **longest** match from it, at most `maxHitLength` (1,000) long: a forward pass over the
  reversed automaton finds the starts, a bounded pass from each start finds its end. An empty match is no hit; a
  pattern that matches only the empty string is `seq.search.pattern.empty`. The second pass counts against the step
  budget, so a pattern that is long at every start (`A*` over poly-A) ends at the budget, said, not at a hang.

### 4. Hits: strands, pieces, circles, maps, order

- **A hit** is `{ target, start, end, strand }`, held columnar (`Int32Array` columns grown by doubling, `Int8Array`
  strands): no object a hit. `strand` is `1`, `-1` or `0` (a protein, or a site equal on both strands).
- **Both strands on nucleotides by default** (`strands: 'both' | 'forward'`): the minus strand is the pattern run over
  the reverse complement of the letters, read backward through `complementCode` (`resolve/codon.ts`, IUPAC; RNA's A
  pairs with U) as a view, nothing copied; its hits mapped back to forward positions. A site the same on both strands
  (`GAATTC`) is one hit, strand `0`. Proteins have one strand. Why the reverse complement of the letters, not of the
  pattern: one kernel, one compiled pattern, any matcher.
- **Pieces are layout, not data** (EC-31): a reversed section, a rotated plasmid or collapsed introns change nothing a
  search reads or finds. A reversed section shows a forward hit drawn right to left; reversing never complements
  (EC-233): the minus strand is the `strands` option's, never a section's direction. A hit in a collapsed stretch is
  kept; stepping to it says where it is (Decision 5).
- **A circle**: a hit may run across the origin, at most one lap (a pattern longer than the circle has none); held
  unrolled (`end` past the length), its Location two ranges (EC-9), one hit in the count.
- **Blocks and assemblies**: each sequence is its own system and its own target; no hit spans two sequences or two
  blocks of a section.
- **Through a codon map: projected, not translated.** A protein pattern searches proteins' letters; its hits on p53's
  protein light p53's gene and transcripts through the codon map as any highlight does (the `Projector`: on the minus
  strand high to low, split at an exon junction inside a codon, EC-15). A protein pattern is not searched in a DNA
  sequence's translations, and a DNA pattern does not reach a protein through its codons: either is an application's
  matcher (SnapGene's six frames). Why: projection is built and exact; translating a genome in six frames is a new
  kind of work with its own bounds, and no version 1 view asks for it.
- **Overlapping hits are all kept** (`AA` in `AAA` is two hits, as the previous implementation's), and the channel
  draws their union. Why: a count that hides overlaps under-reports repeats, and the union draws them right anyway.
- **The order**: targets in the view's order (sections top to bottom; in a section, its blocks', or its rows' display
  order), then start, then end, then strand (`1`, `0`, `-1`). Hits on systems the view does not show come last, in
  address order.

### 5. The viewer: the search model, commands, channel

- **`SearchModel`** (`viewer.search`, DOM-free, `@molstar/core/reactive`; §frontend-state, the UI works headless):
  read-only states `query` (`{ pattern, matcher, strands, within } | undefined`), `status` (`idle`, `running`, `done`,
  `cancelled`, `failed`), `hits` (the columnar hits and their targets' addresses), `count` (hits found, past the held
  ones: Decision 6), `complete`, `current` (an index, `-1` for none), `issues`, `stale`. A hit's address and description on
  demand (`addressOf(i)`, `describe(i)`), never stored a hit.
- **What is searched**: every system a letters, rows or cells layer of the resolved view draws, whole (not only its
  window: IGV's, not JBrowse's), or what `within` names: an address, or `selection` (the selection's Location, its
  parts the targets and their ranges the windows; Jalview's and UGENE's region). A system no section shows is searched
  only when `within` names it.
- **Commands** (`@molstar/ui/commands`, the `search` behavior's; each in `CommandChangeKinds`; the palette lists them
  from the registry, §command-palette):

  | Command | Args | Kind | Does |
  |---|---|---|---|
  | `search.run` | `pattern`, `matcher` (`degenerate`), `strands` (`both`), `within` | `session` | compiles (a compile error is the result's message and issue: nothing changes), cancels a run in progress, runs; resolves when done, cancelled or failed |
  | `search.next`, `search.previous` | `focus` (`false`) | `navigate` | the next or previous hit becomes current and is revealed (pans only, EC-159; `focus: true` zooms to it); round, the wrap said |
  | `search.go` | `index` | `navigate` | that hit becomes current and is revealed; the current one again clears it (EC-161); an index outside the hits refused (EC-160) |
  | `search.clear` | | `session` | cancels a run, clears the query, the hits and both channels |
  | `search.save` | `id`, `label` | `patch` | the hits as a ref (below): one patch, one undo step (EC-162) |

  From no current hit, `next` takes the first hit at or after the start of the window of the first leading section
  showing hits, `previous` the last before its end; then they step in order. Why from the view: on a genome, "next"
  means next from where the reader is (IGV, browsers' find); the count still says "k of n".
- **Reveal through the address model**: stepping runs `AddressModel.reveal` on the hit's address, so every leading
  section showing it pans (a circle by the shortest arc, EC-12), a folded clade opens to its row, a hit in a collapsed
  stretch lands at its break with `HiddenMessage`, one on a hidden track or row says `HiddenTrackMessage`; each a step
  of the view's history (Alt+Left returns).
- **Two highlight channels of their own**: `search` (every hit, a tint, `highlight.colors.search`) and
  `search-current` (the current hit, style `outline`). Set through `HighlightModel`, so they draw, project onto linked
  sections (a protein's hits on its gene) and export as any channel does; listed in the legend's channels row, where
  removing `search` runs `search.clear` (the model listens: a channel and the search never disagree). The selection is
  left alone: the reader's selection is not the search's (Jalview's find selects; not taken). `address.select` on a
  hit's address selects it when a host wants that.
- **Keys** (`view.keys`, a canvas scope, so typing in a host's box is never taken): `search.next` `Mod+g`, `F3`;
  `search.previous` `Mod+Shift+g`, `Shift+F3`. No key for `search.run` (it needs text: the host's box, or the palette's
  argument) and none for `Mod+f`, the browser's find, which a host binds to its own box. `Esc` keeps its meaning (a
  gesture, then the selection); the search is cleared by `search.clear`, the legend, or the host's box.
- **What is said** (`viewer.notify`, the announcer's `notice`; the commands' messages are the same words): on done,
  "14 hits for GAATTC on 2 sequences", "No hits for GAATTC", "10,000 of 182,331 hits held"; a cancelled run says
  nothing; on a step, "Hit 2 of 14: pUC19 396 to 401, minus strand" (the system's label and its ordinal numbering, the
  hover card's words for a place), "wrapped to the first hit" when it wraps. Progress is said once, after 1 s, never per
  batch. Why: a live region queues, and a count read out per chunk drowns the reader.
- **A saved set is a ref** (addressing.md, Bookmarks and location sets are refs): `refs.add` with `to` the list of the
  hits' addresses, `(range :from a :to b)` in the system's ordinal labels (copy link's numbering), a wrapped hit as one
  range with `from > to` on its cyclic system (EC-10), and `label` `Search: GAATTC (14)`. The strand is not in a
  Location, so not in the ref (said in its label when every hit is on one strand). Refused, never cut: past
  `search.maxSaved` (1,000) hits (`seq.search.save.limit`), an incomplete result (`seq.search.save.incomplete`: a
  set must mean what it says). A saved ref is a ref like any: a layer's data, a link's target, `refs.remove` its inverse.
- **Results follow the spec, never re-run on their own**: a patch that removes a searched system drops its hits
  (`seq.address.missing`, as the session's other addresses); new letters for a searched system (a source reattached, a
  sequence's letters edited) mark the result `stale`, still drawn and said ("results may be out of date"), until the
  next run. Why: a search that re-runs itself moves the current hit under the reader.
- **Letters through the store**: the viewer's `LetterReader` asks the store for a target's letters by region (a
  `want` through a wanter, as `selection.copy-letters` does) and waits for them; a region pending is never searched as
  empty (EC-169); a failed read skips that target with `seq.search.letters.failed`.

### 6. Bounds (security.md's "regular expressions in search", now built)

Every bound is checked before the work it sizes, refused with a code (§bounded-work), and a setting (Decision 7):

| Bound | Default | Past it |
|---|---|---|
| a pattern, `exact` and `degenerate` | 1,000 letters | `seq.search.pattern.size` (error) |
| a regex's characters, positions after expansion, nesting | 256, 128, 64 | `seq.search.pattern.size` (error) |
| a regex hit's length | 1,000 | the longest match from a start stops there |
| letters a run reads, all targets | 2^25 (33.5 M) | the target that would pass it is skipped, `seq.search.letters.limit` (warning, hint: `within`) |
| word steps a run takes | 2^28 | the run ends, `complete: false`, `seq.search.steps.limit` (§cancellation-codes' `limit_exceeded`) |
| hits held | 10,000 (`address.limit`) | counting goes on, holding stops: `count` the true number, `complete: false`, `seq.search.hits.limit` (warning) |
| a saved set | 1,000 addresses | `search.save` refused, `seq.search.save.limit` |

- **Counted past the cap, never silently cut**: the count is exact within the step budget, the held hits are the first
  in order, and both are said. Why: "182,331 hits" tells the reader to narrow the pattern; a quiet 10,000 lies.
- **No clock in the core**: the budget is steps, not time (tests have no clock, root `AGENTS.md`); the viewer's run
  yields between batches so a frame is never held, and its signal is read before each chunk.

### 7. The registry and composition

- **A `search-matcher` registry entry** (core's `SeqRegistryEntry`): `{ kind: 'search-matcher', matcher }`, the
  built-ins unnamespaced and in the default registry, an application's namespaced (`<id>:<name>`), checked at
  composition with the other entries (no built-in shadowed, conflicts one error; extensibility.md); a run naming an
  unregistered one is `seq.search.matcher.missing`, naming it, no fallback. An extension adds one as registry entries
  (`defineExtension`), never on import.
- **What an application builds on**: `viewer.search` (states), the commands, `describe(i)` and `addressOf(i)` for a
  results list, the registry for its matchers. Nothing else is needed for the Storybook box (Decision 8), and a test holds it
  (Consequences).

### 8. The Storybook box

`storybook/src/_stories/support/search-box.ts`, beside `support/toolbar.ts` (the stories' host composition, as
`mountHiddenTracks` is placed there): an input, a matcher choice (IUPAC, exact, regex), a strand switch on
nucleotides, the count ("2 of 14", "14+" while running, "No hits"), previous and next buttons, Save and Clear;
Enter and Shift+Enter step, Esc in the box runs `search.clear`; the text runs `search.run` after it rests 200 ms, the
previous run cancelled (the previous implementation's). Built from `@molstar/ui` components and tokens (ux-foundations
section 9) on the public surface only. Used by the p53 example's "search" story (features.md's "p53: search") and a plasmid
story (a site across the origin, both strands). It is the stories' box, not the library's: no package imports it.

## Defaults

What a reader sees with no setting: the hits a tint of `highlight.colors.search` (`token:cat-3`, distinct from the
links' `cat-5` and from the accent: data and chrome never share a role) under the data, at 0.3 as every channel; the
current hit an outline of the same color at full strength, concentric with what it surrounds (ux-foundations section 5);
where the search and another channel light one position, the stripes rule (section 6 there), the current hit's outline drawn
over them, never a stripe of its own. Empty: nothing drawn, the legend's channels row as its mode says. Crowded (a 1 Mb
genome zoomed out, thousands of hits): the union's tint at the density step's resolution, as any channel; no tick
marks on the navigator (not built). Dark: the tokens' dark variant.

**Mock**: `Seq* (next)/Design/Mocks/Search`: a fake plasmid of 2,686 letters and a crude box (input, count, previous,
next) driving a stub of the model; hits tinted, the current outlined, one across the origin; played at the gate,
deleted by the merge that lands the built box. **Precedent**: SnapGene's find bar.

**Provisional**: the color token and the current hit's look, settled by the mock at the gate.

## Configuration

| Item | Default | Why | Alternatives |
|---|---|---|---|
| `search.matcher` | `degenerate` | IUPAC is what a nucleotide pattern means (SnapGene); an ACGT pattern is exact anyway | `exact`; a registered name (not built) |
| `search.strands` | `both` | a site on the minus strand is a site (SnapGene, IGV, UGENE) | `forward` (not built) |
| `search.maxPattern` | 1,000 | primers, guides and motifs are far shorter; a pasted gene fits | larger, raising the step cost (not built) |
| `search.maxRegexChars`, `maxRegexPositions` | 256, 128 | four words a step; `include` patterns' 256 | more words (not built) |
| `search.maxHitLength` | 1,000 | UGENE's regex result length | (not built) |
| `search.maxLetters` | 2^25 | a bacterial genome and a 10,000-row alignment fit; a human chromosome does not, and says to use `within` | larger, or a source's own search (not built) |
| `search.maxSteps` | 2^28 | about a second of work at worst on a laptop, to measure (Provisional) | (not built) |
| `search.maxHits` | `address.limit` (10,000) | what a channel and stepping hold; the count goes on | (not built) |
| `search.maxSaved` | 1,000 | a ref's addresses stay readable and its patch small | (not built) |
| `search.chunkLetters` | 2^20 | a batch a frame's work | (not built) |
| `search.channel`, `search.currentChannel` | `search`, `search-current` | channels of their own (features 181) | (not built) |
| `highlight.colors.search` | `token:cat-3` | apart from the links' channel and the accent | (not built) |
| `view.keys['search.next']`, `['search.previous']` | `Mod+g`, `F3`; `Mod+Shift+g`, `Shift+F3` | the platforms' find next and previous, in the canvas only | (not built) |

## Rejected alternatives

- **A JavaScript `RegExp`** (the previous implementation's): backtracking, so a spec's or a link's pattern can hang the
  page (`(A+)+B`); a timeout needs a worker and a clock. The automaton bounds the cost by the letters.
- **Overlap matching for degenerate letters** (a sequence `N` matching any pattern letter): floods gaps and drafts with
  hits.
- **Searching only the window** (JBrowse): the count and stepping change as one pans; there is no "next" past the
  view.
- **The current hit as the selection** (Jalview): it costs the reader's selection on every step; `address.select` on a
  hit is one command for a host that wants it.
- **Results as a layer or new tracks** (IGV): a session's search would edit the spec. Saving makes a ref, and a ref is a
  layer's data when the reader wants a track.
- **Translating DNA in six frames for a protein pattern** in the library: a new kind of work for no version 1 view; an
  application's matcher.
- **A callback-taking matcher port**: §return-observables; the kernel's callback is the hot loop's exception only.
- **A deadline in the core**: tests have no clock; steps bound the work, the viewer yields.

## Consequences

- **Core** (`@molstar/seqx-core/search`, a new entry): the kernels (exact, degenerate, the regex compiler and
  automaton), `kernelMatcher`, the IUPAC sets (`resolve/codon.ts`'s `Iupac` and `complementCode` shared, not copied),
  the `search-matcher` entry, `SearchBounds`, the issue codes (`seq.search.*`, cataloged with hints). The `_spec` gains
  a "Search" section; its Language › Bounds row for regular expressions now points there.
- **Viewer**: `SearchModel` on the viewer, the `search` behavior (commands, keys), the two channels, the store's
  `LetterReader`, the `search.*` items in the configuration and its table (`config-table.test.ts`); the `_spec` gains a
  "Search" section and the commands' rows.
- **Designs amended in the same change**: security.md (the regex row built; its In the specs row), locations-and-
  queries.md (`search` there stays a later form, now pointing here), ux-foundations.md (the matrix row below; the
  channels), extensibility.md (the matcher entry in its table), features.md (182 to building, then done; 180, 181 and
  207's "search not" removed; §features-move-with-merges), the Storybook's Design pages ("Search", the Features line
  with its story).
- **TODO.md**: the "Bounds" line's regular expressions removed; the `(search …)` form, several result sets and
  translated search recorded as followups.
- **The interaction matrix**, a new row (a test names each cell's gesture):

  | | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
  |---|---|---|---|---|---|---|---|---|---|
  | **search hits** | no: the selection is the reader's (`address.select` on a hit) | — | no: a new run replaces | `search.next`, `search.previous` (`Mod+G`, `F3`; reveal), `search.go`; `focus: true` zooms | Alt+Left | — | — | `search.clear`, the legend's channel | `search.save` (one patch); `refs.remove` its inverse |

- **Edge cases**, each a test by name: EC-161 (stepping "2 of 14", the current again clears, saved as a set), EC-9
  (a hit across the origin, two ranges, one hit), EC-10 (its saved address wraps on the cyclic system only), EC-12
  (revealed by the shortest arc), EC-15 (a protein hit projected onto a minus-strand gene, split at a junction), EC-31
  and EC-233 (a reversed or collapsed section finds what a plain one does; reversing never complements), EC-37 (the
  system's alphabet decides), EC-155 (the keys never a browser's zoom, never in a host's input), EC-157 (the search's
  channels cleared apart from the others), EC-159 (stepping pans only), EC-160 (`search.go` outside the hits refused),
  EC-162 (a save one patch, one undo step), EC-169 (a pending region never searched as empty), EC-76 (a regex's
  nesting bounded).
- **Tests**: each matcher against a naive oracle on seeded random letters and patterns (both strands, a circle,
  overlaps), the regex compiler's refusals (backreference, lookaround, lazy, size), a backtracking pattern (`(A|AA)*B`
  over 100,000 `A`s) costing what `AB` costs in counted steps (quick lane), the bounds each refused before allocating,
  cancellation mid-run, stepping and saving headless, the Storybook box importing only public exports, the box's DOM
  test.
- **Performance representatives** (§performance-tested; quick representatives counted, heavy bounds in `itWithin`):
  - **a 1 Mb genome** (seeded, `core/src/search/_test/`): `degenerate` `GANTC` and `exact` `GAATTC` on both strands, and
    `regex` `GG.{2,8}CC`: quick, 100 kb with counted steps and no allocation a hit (an allocation guard on the kernels,
    `*.quick.ts`); heavy, 1 Mb under 30 ms (exact, degenerate) and 100 ms (regex), Provisional until measured;
  - **a 10,000-row alignment** of 300 columns (seeded): a motif across ungapped rows, about 3 M letters: quick, 1,000
    rows counted; heavy, 10,000 rows under 150 ms, and a step after the run under 1 frame (the reveal, not a rescan).
- **Changelog**: an entry file in `.v6/changelogs/` and the `[Unreleased]` line (root `AGENTS.md`).

## Tensions

- **A cap against "never a truncation"** (addressing.md, Bounded): the held hits stop at 10,000 while the count goes
  on. Settled here as not a truncation: the count is exact and said, the result says it is incomplete, and an
  incomplete result cannot be saved. Standing: a reader cannot step past the 10,000th hit without narrowing.
- **The current hit's outline against the stripes rule**: an outline over stripes is a second look for one rule;
  carried until the mock is played.
- **Whole systems against streamed genomes**: searching a chromosome read by region reads it whole, bounded by
  `maxLetters`; past it the reader must narrow (`within`). A source's own search (as `findNames` is a source's
  capability) would lift it: not built, carried.
- **A strand the Location cannot hold**: a saved set loses which strand each hit was on. Carried until a ref can carry
  fields (a strand column) or the set is saved as a ranges annotation.

## Provisional

- **Kernels in `@molstar/core`?** The kernels are generic (bytes, no alphabet): root `AGENTS.md` says to suggest such
  code for `@molstar/core` (`algorithms/sequence/`, beside Needleman–Wunsch) with a test, and to ask before changing a
  public API. Recommended: start in `@molstar/seqx-core/search`, move when a second package (Mol\*'s sequence panel)
  wants them. Settled by the author's word.
- **Saving an alignment's hits**: they are rects (a row by its key × the columns), and rect addresses are past
  version 1 (the author, 2026-11-23, trigger "a link to a rectangle"). Recommended: version 1 searches, lights and
  steps an alignment's hits, and `search.save` on them is refused (`seq.search.save.rects`) until rect addresses land;
  the alternative is to take this as the trigger and build rect addresses now. Settled by the author.
- **The `(search …)` query form** (locations-and-queries.md, "later, with the examples that need them"): a live set a
  link or a layer could name, run by the same matchers. Recommended: not in version 1 (a form evaluated at hover time
  over a whole genome is a different bound); trigger, a link or layer that must name hits live.
- **One search a viewer**: the previous implementation held named result sets. Recommended: one; several when a host
  asks (trigger).
- **From where `next` starts** (the window's start, Decision 5) and **reveal over focus** by default: settled by the
  novice pass on the p53 and plasmid stories.
- **The numbers**: `maxSteps`, `maxLetters` and the heavy bounds, settled by measuring the representatives on the CI
  runner.

## Self-check

- **The claim**: an application builds its own search UI on the library's model, commands and registry, and the
  library's search is safe on untrusted patterns.
- **Negated**: an application needs more than the public surface (a box reaching into the viewer's internals), or a
  pattern hangs the page. **What the negated world explains**: the previous implementation's search part lived in the
  viewer because its results needed internal hooks; its regex trusted the engine. Here the Storybook box is built
  outside the viewer and a test holds its imports to the public exports; the regex engine is an automaton with a
  counted budget, and a test runs a classic catastrophic pattern.
- **The inverse and the symmetric cases**: run and clear; next and previous; a save and `refs.remove`; the forward and
  minus strands; a hit on the protein lighting the gene, and a search on the gene (its own letters) — not lighting the
  protein through translation (said, an application's matcher); removing the channel in the legend clears the search
  (the model listens), not only the reverse. **The plural**: many hits (the cap, counted), many systems (the order),
  many searches (one, Provisional).
- **Built or planned**: everything in Decisions 1 to 8 is this lane's; the rect save, the query form, translated and
  approximate search, several result sets and a source's own search are carried above, each with its trigger.

## In the specs

| Decision | Where |
|---|---|
| Kernels sync, a visitor, no allocation a hit | core › Search › Kernels |
| The matcher port: compile sync, run an observable of batches, cancellation | core › Search › Matchers |
| `kernelMatcher`: chunks, overlap, batches a prefix of the order | core › Search › Matchers |
| The built-in matchers: case, stored letters, U as T, degenerate subset matching, gaps skipped, the alphabet check | core › Search › Built-in matchers |
| The regex language, its refusals, the Glushkov automaton, longest hit a start | core › Search › Regular expressions |
| Hits: columnar, strands, palindromes once, circles, blocks, overlaps, the order | core › Search › Hits |
| Through a codon map: projected, not translated | design-only: a non-feature (said in the viewer's Search section's projection line) |
| The bounds and their codes; counted past the cap | core › Search › Bounds |
| Regular expressions in search (security.md's row) | core › Search › Regular expressions |
| The `search-matcher` registry entry | core › Context › Registry entries |
| `SearchModel`: states, what is searched, `within` | viewer › Search |
| The commands, their args and kinds | viewer › Commands |
| From where `next` starts; reveal through the address model | viewer › Search |
| The two channels, the legend clearing the search | viewer › SelectionModel and HighlightModel |
| The keys | viewer › Bindings |
| What is said | viewer › Input › Announcements |
| A saved set is a ref; refused past the bound or incomplete | viewer › Search |
| Results follow the spec; stale on new letters | viewer › Search |
| Letters through the store, pending never empty | viewer › Search |
| The configuration items (the bounds' defaults are core's Search › Bounds) | viewer › Configuration |
| The Storybook box on the public surface | design-only: the stories' host composition, not a package's |
| The representatives | core › Search › Performance |
| Rejected alternatives | design-only: rejected alternatives |

# The previous implementation's edge cases

An inventory of the edge cases and rich behaviors the previous implementation specified, tested or demonstrated
(2026-11-17, from its `_spec/`, `todo.md`, test names and stories). It is the checklist that keeps the second
implementation from smoothing them over ([design/edge-cases.md](../../design/edge-cases.md)). Each item: the
behavior, why it is tricky, and where it was evidenced. Paths are in `packages/sequences/`: **S/** `_spec/`, **C/**
`core/_test/`, **V/** `viewer/_test/`, **SB/** `storybook/src/_stories/`, **TD** `todo.md`. Items are numbered so a
design or a test name can cite them (`EC-23`). Section 9 adds the cases the JBrowse reply raised.

## 1. Coordinates and maps

1. Indices 0-based half-open; labels inclusive, converted once at load; never mixed in one table. S/format.md#positions
2. Label order: integer, then no suffix, then letters (shorter first, alphabetic, case-blind), then `.n`; `100AA` after `100Z`. *Labels a numbering lacks must still be placed.* C/model/numbering
3. Four numbering methods: `offset` (may skip zero), `segments` (holes, no overlap), `labels` (sparse suffixes, Kabat/IMGT), `names`. C/model/numbering
4. A label outside int32 is refused, never wrapped; a repeated label names its first position. C/model/numbering
5. Insertion-code columns complete positions; `""`, `?`, `.` mean none; a code on a range is a form error. C/format/insertion-codes
6. A range end the numbering lacks snaps inward; an unknown single label drops the row and lists labels that exist. *Errors must teach.* S/format.md#positions
7. IMGT labels are not in index order; a backwards range is an error saying to write beg ≤ end. C/format/insertion-codes
8. A unit's own numbering beats its molecule's of the same name; section numbering is looked up along paths, never across a codon step. S/format.md#sections
9. Cyclic: an extent never wraps; a wrapped range is two runs; a wrapped annotation row is two rows of one item. C/model/cyclic
10. `start > stop` wraps only on a cyclic target and only on side a of ranges, values, matrix; pair ends and grid sides never wrap. S/model.md#cyclic-systems
11. A map segment crossing the origin is split when built; a CDS may cross the origin on either strand. C/model/cyclic
12. Focus and scroll-to on a cyclic system use the shortest arc. *The naive hull is the whole circle.* V/cyclic
13. A wrapping unit range is unsupported, a wrapping block range is allowed (deliberate asymmetry). C/model/cyclic
14. A cyclic label range must be less than a whole lap. C/format/insertion-codes
15. Codon maps: CDS ranges in transcript order; minus strand highest first, read downward and complemented; an exon junction may split a codon. S/model.md#codon-maps
16. Codon status: T/U alike, ambiguity → X, first codon M, TGA as Sec and TAG as Pyl, NCBI tables, stop codons pair with nothing. C/model/codon
17. A CDS cannot start mid-codon, its length is a multiple of 3; overlapping or out-of-order ranges are reported. C/model/codon
18. Through a codon step: residue → 3 nucleotides, nucleotide → 1 residue, UTR nucleotide → nothing. C/model/codon
19. At most one codon step per path; numbering labels never cross one. S/model.md#paths
20. Maps are one-to-many backward (a GFP–GFP fusion), and that is not ambiguous. C/model/coordinates
21. Partial maps and chimeras: what no map covers maps nowhere, never an issue; a mapping that pairs nothing is. C/model/mapping
22. Paths: shortest only, at most 4 counted steps, unit↔molecule steps free, at most 16 candidates. S/model.md#paths
23. Several shortest candidates that compose to different maps are skipped as ambiguous; equal compositions count once; empty paths dropped first. C/model/mapping
24. Two mappings joining the same two systems make every path through them ambiguous. S/format.md#mappings
25. Copy rule: no path climbs from a unit to its molecule and down to a sibling; a longer legal path beats a shorter illegal one. *Part of the search, not a filter.* C/model/mapping
26. A path ending at a ranged unit is clipped to the unit's range. C/model/coordinates
27. Composed paths are sparse (a protein through a chromosome is the size of its CDS) and composed at resolve, never on hover. C/model/paths
28. A boundary maps before its element's first image; if unmapped, after element k−1's image; else nowhere. S/model.md#coordinate-maps
29. An insert or deletion seen from the other side lights a caret at the next pair; a boundary inside a codon lands after its residue. V/insert-carets
30. Closure is one step: a nucleotide gains its residue, a residue its same-phase isoform's residue; a cell never lifts to its column; closure reads the whole view's maps. S/model.md#closure
31. A collapsed intron, a reversed strand, a rotated plasmid are layout maps, not data; two pieces of one unit are two blocks over one system. S/model.md#pieces-and-layout-maps
32. A reversed piece: shown position k is index `hi−1−k`, runs in flipped order; filter, reversal and origin compose seamlessly. V/axis
33. A filtered block of 150,000 single-position runs lays out in either direction. V/layout-edges
34. A block lays out at most 2^31−1 positions; past it, left out with a message to use a range. V resolve
35. Units, rows, columns are systems; a complex adds no arithmetic; `item` has no domain; `external` is the host's. S/model.md#coordinate-systems
36. Ids append-only within a lineage; a key that returns gets its old id; across lineages loci remap by key. C/model/set
37. Alphabet detection: nucleotide only if all IUPAC and ≥ 90% ACGTUN; HHHHHH is protein; non-ASCII rejected. C/model/molecule

## 2. Data

38. `-` and `.` are gaps; insert columns from `#=GC RF` or all-lowercase columns; all-lowercase input is plain letters; A3M lowercase maps to no column. *Three readings of lowercase.* C/model/alignment
39. A header `NAME/27-137` numbers the row from 27. C/model/alignment
40. Ragged rows, bad letters, descending headers, a bad reference row are reported; Clustal blocks must list the same rows. C/format/load
41. Lazy alignments: rows not held answer none and never throw; only endpoints, annotated and reference rows are held. S/sources.md
42. Path search over 20,000 rows goes through bound rows only; a row bound twice is ambiguous. C/model/mapping
43. Tree leaves match by label, else up to the last `/`; a failing tree keeps its alignment; nodes in any order. C/model/alignment
44. Scale: 10,000×300, 1,000×2,000, 35,000 residues, 3M nt, 10M lazy rows. SB/large-documents
45. Rows sharing an `id` are one item; where rows overlap, the last covering row wins. C/model/annotation
46. Zero-length (boundary) items have their own draw and hit rules. S/model.md#annotations
47. Rows that cannot be placed are left out with instructive issues; an empty column set is valid. C/model/annotation
48. A target that is a sequence reaches every unit over it; an alignment its columns; a row its residues. S/format.md
49. `fields` remaps file columns (swaps, chains); text cells keep their text (`01` stays `01`). C/format/tables
50. `__proto__` and `constructor` are ordinary column names. C/format/proto-keys
51. UniProt features past the sequence end dropped, running past it clipped, counted in one note. C/format/load
52. NaN means missing, never zero; a comparison with missing is false, `not` of it true. S/seqql.md#types
53. Numbers are exact JSON doubles; Float32 columns compare as their shortest decimal. V/values-exact
54. Dense values: one array per target; more values than positions is reported. C/model/annotation
55. Values draw at their value, `domain` bounds only the height (the old docs conflicted: one said clamped). S/format.md
56. Shared domains: union per (group, field); an explicit domain wins in its layer; a signed member makes the group diverging. C/format/domain-group
57. A line changes color halfway between positions and holds over half cells at block ends and gaps. S/format.md
58. Grid values a-major, `null` missing, typed arrays by reference; a wrong length drops the row and says how the count is formed. C/format/grid
59. Two grids on one block pair merge, later over earlier; a size cap per block pair. V/dot-plot-more
60. Matrix density `majority` (tie to the first color, alpha the share) or `maximum` (rare signal shows). V/matrix-dense
61. Blocks with different state counts keep their own stride. V/matrix-stride
62. Pairs: ends on different targets in their own numberings, self-pairs, open categories, found from either end, contacts fold by copy group. C/model/annotation
63. Several pair kinds at one cell split into glyphs, "+" when they don't fit; hover lights that kind only. V/pairs-kinds
64. A ligand is one position, or one per heavy atom; flows as chips in text; hydrogens rejected. V/projection
65. Copy groups derived; `copies: one` merges into one block labeled `×N`; hover on it stands for every member. C/model/set
66. Streamed data on a merged block through a path per copy is not placed, with an issue. C/format/streamed-paths
67. 200 identical copies of secondary structure draw what one does. V/secondary-copies

## 3. Format, validation, patches

68. Ref syntax; only the root and data may omit one. S/format.md
69. Derived refs: the prefix kind checked statically; existence known only after load. S/format.md
70. Reference cycles reported on every node in them; the rest runs. V/graph
71. Nothing cascades from a parent; dependencies are exactly the typed ref params. S/format.md
72. Missing pieces never fail a document: unregistered provider, format or theme → unavailable or fallback, reported once. C/format/load
73. A newer minor version: unknowns are warnings, drawn by fallback; another major is rejected. C/format/validate
74. Fallback chains end at letters, ranges or omit, checked against the section, reported. C/format/resolve
75. Theme nesting capped; extension types state a kernel fallback. C/format/validate
76. Nesting capped (256, queries 64): a document 20,000 deep errors at its path, no stack overflow. C/format/untrusted
77. Malformed children give a type error, never an exception. C/format/malformed
78. Patches: shapes checked, then ops applied, then the result validated; errors pathed at the op. C/format/malformed
79. Anchors are refs, `first`, `last`; `set` merges shallowly, `null` restores the default, never changes a kind. S/format.md#patches
80. Removing a node with refs to it fails the patch (no cascade). S/format.md#patches
81. A patch's `base` must match the revision; `null` applies to any. C/format/patch
82. Patches share untouched subtrees and return exact inverses; `diff` then apply is tested both ways on random edits, including a ref reused under another kind. C/format/patch
83. Diff must not hang on a duplicated ref. C/format/malformed
84. Snapshots apply to the base; a base change that breaks one marks it invalid and the player skips it. V/snapshots
85. Entering a snapshot sets its view; only sections whose focus changed animate. V/snapshots
86. Issues of one cause reported once with a count; first ten paths kept; at most 100 per node, most severe first. S/format.md#validation
87. Severity means one thing everywhere: error invalid, warning a loss, info never changes the picture. S/format.md#validation
88. An unknown key lists the valid keys; a bad value lists allowed values. SB/style InvalidKey
89. A stale query text is reported; the expression wins. S/seqql.md
90. Binary encoding: 32-bit integers stay integers, non-finite → nil, typed columns as extensions; round-trips to the same text. C/format/msgpack
91. Archives: entry, size and path limits, checksums. C/format/svsx
92. TSV/CSV: trimmed cells, quoted kept, duplicate header an error, any line ending; fuzzed. C/format/load
93. A host's echo of its own commit keeps history even reformatted. V/snapshots

## 4. Layout and views

94. Text layout wraps whole groups of 10; a short last group on its own line where needed; labels above or inline. V/text
95. Inline group labels take part in line breaking; continuation lines have none. V/text
96. Text lanes per layer kind, compacted per line, empty bands collapsed (not values or color). V/text-tracks
97. Alignment rows in text (ESPript style): label column capped with ellipsis, row's last residue number per line, group numbers from the reference. SB/text-view AlignmentFigure
98. Centered and max-width text move hits and margins with it. V/wrapped-spacing
99. Ligand chips flow like words. V/wrapped-spacing
100. Frames: open sides where a line or the view cuts the item; title on the starting line only; a collapsed clade frames its wedge. V/frames
101. Default blocks: polymer units, then ligands, then own units, then the first alignment. S/format.md#blocks
102. Side-by-side vs stacked blocks; 24 blocks share gaps, a quarter of the area kept for positions when zoomed out. V/layout
103. Groups collapse; a summary track shows only while collapsed. V/union-summary
104. Split rows: fixed order, empty rows shown or hidden, guides only on rows ≥ 4 device px, one split layer per track. V/split
105. A split row on a unit of its own takes that unit. V/split-units
106. Rows layer: order by data, tree or key table (missing last), filters, subset first, tree only in tree order. V/keyed-rows
107. Data columns (number, fraction, text, category) yield to a narrow width down to a floor. SB/style DataFloor
108. Rings: tracks outside in, polar hits, sectors reduced per angular column, labels outside with leaders, dropped by priority. V/ring
109. Dot-plot aspect bounds, with a rule that stops flicker from a scrollbar. V/dot-plot-aspect
110. Navigators: strip, miniature, ideogram, vertical overview; idle when there is nowhere to go; a floored window drawn wider. V/navigator
111. Viewport across updates: an equal axis keeps it; a changed one keeps center and zoom when still shown, else fits. V/viewport
112. A move before the section has a width is kept as an intent; a section hidden for lack of width keeps its zoom. V/viewport
113. Scroll groups need equal axes and join across contexts; follow goes through a map at scale 1, 3 or 1/3. V/follow
114. Follow chains; followers never lead; a ring leads only through its window. V/follow
115. A stacked section anchors a zoom on the block under the point. V/layout-edges
116. Text caps at the label size; a row shows text for all its columns or none; then dots per residue, then bars per run. V/letter-cap
117. Three-letter codes need room inside the narrowest piece in view (a codon split by an intron). S/viewer.md
118. Unobserved residues gray in cells, text and zoomed out. V/unobserved-gray
119. A track ordered before the letters draws above them in both layouts. V/track-order

## 5. Rendering

120. Glyph atlas matches fillText at each DPR; least recently used page evicted, never one in use; dropped on font or token change. V/canvas/atlas
121. An animated zoom uploads each atlas page at most once; nothing re-measured after zooming in and out. V/canvas/zoom-text
122. Letters centered on their cells per role, size and DPR. V/canvas/letter-centering
123. A device pixel hits the cell drawn on it; a glyph's overhang into a gap is hit too. V/hit-drawn
124. SVG/PNG parity with Canvas2D: same ink box (offsets, no baseline attributes), committed font metrics, valid PNG checksums. V/svg/metrics
125. Export never changes the live viewport; streamed regions awaited at export size; a timeout says unsettled. V/canvas/export
126. Render CLI: files confined to a root, no network, distinct exit codes. V/svg/cli
127. Boundaries over letters: a pinned triangle and hairline, a hairline alone when tight, end boundaries inside the data area. V/boundary-pins
128. Line styles (solid, dashed, dotted, wavy, double, box, highlight), dash phase anchored to content; an invalid style draws solid, reported once. V/line-styles
129. Arcs stay in their band, move continuously on pan, reach off-window partners a little past the edge, nest by span, label only under an apex in reach. V/arcs
130. Secondary structure: arrowhead on the piece holding the C-terminus (reversed, across the origin); square ends where cut; minus-strand codon maps pointed wrong (a known flaw). V/secondary-strands, TD
131. Density: stacked or majority when a device column holds over 64 items. V/density-band
132. Legends: ramps for continuous, steps for discrete; identical groups merge; "+n more"; hashed values can't be listed. V/legend
133. One palette per layer ordered by first occurrence over all rows (a manifest dependency the new design drops). C/format/palettes
134. A value not in the vocabulary takes the FNV-1a hash of its UTF-8 bytes modulo the list; dictionary palettes give it the missing color. C/format/palettes
135. Rule themes: a measure that doesn't apply is NaN and fails its clause; Clustal X and ESPript as rule tables (ESPript ignoring gaps, a known flaw). V/themes, TD
136. Color schemes: dark variants, only fitting schemes offered, a pick survives patches, cleared by a new document. V/color-schemes
137. Size channels change the fill, never the slot or hits; an undrawable channel is reported. V/size-channels
138. Hover emphasis: an outline with ≥ 3:1 contrast; none for external hovers or in exports. V/hover-emphasis
139. Skin: dark below luminance 0.179, AA contrast for ink, conflict color ≥ 4.5:1. V/skin
140. Hide then show, or resize, equals a fresh mount; a band copy on pan matches a full redraw. V/canvas/resize

## 6. Interaction and editing

141. One hover owner at a time; a stale owner's clear does nothing (a late leave never clears a 3D hover). V/channels
142. The hovered element is named in its own system; its locus resolved lazily once per version. V/channels
143. A resting pointer is re-hit once per frame when the view moves under it. V/resting-hover
144. Hovering a nucleotide lights its residue's whole codon; a residue split by an intron lights every piece. V/projection
145. A dot-plot cell is two positions: crosshair hover, rectangle drag, mirrored on a reversed piece. V/dot-plot
146. Header hover lights a track's shown items up to a limit, else the row band. V/header-hover
147. The data-hover row band is drawing only: it changes no state. V/row-hover
148. Clicks: Shift toggles or extends (anchor kept, across cyclic runs), Ctrl adds; a single click waits when a double click is claimed. V/selection
149. A drag over nucleotides snaps to whole codons on release; a drag over a million positions closes by runs. V/input-findings
150. Autoscroll past an edge speeds with distance, time-based, stops at the end. V/text-select
151. A text drag selects line by line within its block and never pans; cross-block selection deferred. V/text-select, TD
152. A press always ends: a lost release is synthesized; one pointer at a time; outside points clamped, overshoot kept. V/input
153. Input claims by priority; an unclaimed drag does nothing. V/input
154. Keyboard cursor crosses block gaps, steps through filtered, reversed and cyclic pieces and a million lazy rows, clamped on relayout. V/key-cursor
155. Key bindings match modifiers exactly; browser zoom keys and Tab never captured. V/keyboard
156. Context menu from the keyboard, deduplicated against the browser's. V/context-menu
157. Highlight channels cleared independently; multi-kind targets get numbered channels. V/channels
158. Reveal-on-highlight waits for a dwell, never under the pointer or during a press. V/reveal
159. Focus zooms; reveal only pans, and does nothing if the target is shown. V/reveal
160. Range arguments outside the domain are refused, never clipped. V/commands
161. Results step through hits ("2/14"); the current again clears; results save as a location set. V/results
162. Bookmarks are one patch, one undo step, with optional color and line style. V/location-sets
163. Document layers base → snapshot → local → preview; local edits tagged to their snapshot; an invalid preview is not drawn. V/document
164. Undo and redo round-trip random sequences; an inverse that no longer applies is refused. V/document
165. Removing a behavior mid-gesture leaves nothing lit or stuck. V/input-findings
166. The hover card bridges the gap into itself, closes on zoom, survives a reassembly. V/hover-card
167. Large selections linked to 3D in size-gated chunks. S/libraries.md
168. Switching views shows the extents, preferring a view that holds them directly. V/show

## 7. Streaming and sources

169. Regions pending, loaded, empty or failed; pending is never drawn or hit as empty. V/streams-port
170. Levels: the coarsest bin that fits a device column; unlisted falls back to the coarsest no wider. V/streams-levels
171. A zoom gate: bins show density; "load anyway" forces. V/streams-gate
172. Requests view → margin → prefetch (pan direction, next zoom, neighbor snapshots), at most 6 in flight. V/streams-prefetch
173. Cancellation with a grace period; the oldest stale request goes first. V/streams-jitter
174. Retries with backoff then a cooldown, reported once; "not served" never retried; a malformed answer is a failure. V/streams-port
175. Streamed data never enters the document or its history. V/streams-port
176. Dynamic nodes identified by (source, key); a conflict keeps the held node; eviction quiet, ids return. V/dynamic-nodes
177. Item coordinate systems made on demand, length from the map. V/item-systems
178. Templates rerun only themselves on arrival. V/templates
179. Nested sources stop at a max depth. V/nested-parallel
180. A streamed codon label shows "loading…", high to low on the minus strand, `|` at an intron. V/streams-paths

## 8. Integration

181. Mol\* structure links with ids per chain instance (repeated ids across structures crashed: a known flaw). S/libraries.md, TD
182. Markdown links: a leave clears only its own channel; a removed link releases once. C/links
183. Nightingale: fallbacks and known losses; huge sections omitted. S/libraries.md
184. JBrowse demo: gene folding, navigation across chromosomes, zoom levels at whole-chromosome views. V/jbrowse

## 9. From the JBrowse reply

Cases the JBrowse reply raised by drawing SeqViewSpec documents in JBrowse and probing the validator
([jbrowse-reply](../jbrowse-reply/README.md)). Paths are in the seqviewspec-reply repository: **R/** `README.md`
sections, **X/** `examples/README.md` (Gaps table), **P/** `spike/jbrowse-consumer/README.md`, **Q/** `scratch/probes/`.
They are cases to handle, not only requests: each gets an outcome like any other item.

**Portability: what a document must say for itself**

185. Every example drew nothing once its host-private block was removed (the stripped test). *A spec must open with no host code.* P/The stripped test
186. Where a file is, its format and its index (`.tbi`, `.fai`, `.gzi`), never assumed as `<url>.tbi`. P/What the document did not tell JBrowse
187. Which assembly a genome is (an accession or name, not only a ref). P/same
188. Sequence names and lengths without a manifest; whether `chr17` exists is knowable from the spec or the file. P/same
189. Sequence-name aliases (`1`, `chr1`, `NC_000001.11`; `MT` vs `chrM`). P/same
190. Which file column is which field (a BigBed's `geneName` as `label`). P/same
191. A categorical vocabulary and its order without a manifest; palette order by first occurrence cannot be computed by region. P/What a reader … has to reimplement or guess
192. A size gate: how much a reader may load before asking (`max_bytes`). P/What the document did not tell JBrowse
193. Pixels per row, label fonts, lane heights are the viewer's, and a second reader picks its own; track height in px is asked for. X/Gaps (height_px)
194. Color names and named palettes need values a reader can use, not only names. P/What the document did not tell JBrowse
195. The default numbering of a streamed sequence (1-based labels) must be stated, not learned from the resolver. P/What a reader … has to reimplement or guess
196. What a streamed item's columns are (an id as `name@start-end` or `pos:ref>alt`, a row per BED12 block, a row per ALT allele). P/same
197. Dictionary-coded columns: a code is a 0-based index into the dictionary (stated, not guessed). P/same
198. An annotation's default label: from the source, else its id. P/same
199. How a fallback draws data of another shape (pairs as points, a matrix as values) is stated per pair of types. P/same
200. Derived-ref and id characters must cover real names: HLA contigs such as `HLA-A*01:01:01:01`. R/Notes on the format itself

**Value ranges and colors without the whole data**

201. An open value domain means the data in view (the reply's question; JBrowse also clips outliers past a quantile). R/Where Seq\* and JBrowse will differ
202. A domain with one open end (`[0, null]`). X/Gaps
203. A value not listed takes one color everywhere (a hash of its text), the same across tracks and readers; the old docs disagreed (missing color vs hash). R/Where Seq\* and JBrowse will differ
204. A legend of what was painted vs of the whole vocabulary: say which, and list hashed values when they appear. R/same
205. A palette's value domain separate from the y domain of a bar or line plot. P/What JBrowse could not honour

**Comparing genomes**

206. Strand on every map and on pairs rows: an inversion between genomes; a reversed side on a pairs row is valid, not an error. R/Comparing genomes, X/Gaps
207. Several maps between one pair of systems (a mapping per region or per chain file block) without every path becoming ambiguous. X/Gaps (mapping per region)
208. Pairs read by region (PAF, PIF, BEDPE), not whole with a manifest. R/Comparing genomes
209. Pairs drawn between two sections (ribbons) and a dot plot whose y axis is another section. X/Gaps
210. Two assemblies as blocks of one axis or one circle (peach and grape on one ring). R/A probe: two genomes on one circle
211. A follow through a map that arrives by region (a dynamic map). X/Gaps
212. Liftover through a chain file, by reference; CIGARs with S, H and N operations, and CIGARs that cover only part of each sequence. R/Smaller points
213. A multiple alignment (MAF) laid over a genome axis. R/Smaller points

**Many samples, signal and statistics**

214. Many signals from one file: `split` on values (a column per sample, or a column naming the sample). R/Signal and statistics
215. Log scale on values and on continuous palettes (TPM over 0 to 10^4). R/Signal and statistics, Q/probes-data (E3)
216. Points at a height (a Manhattan plot): a zoomed-out bin must not average away one significant SNP. R/Signal and statistics, Q/probes-data (M1 to M3)
217. A reference line (a significance threshold) and a baseline bars rise from; a scale's title and units. R/Smaller points, Q/probes-data (M5)
218. Genome-wide values colored by chromosome over a whole-assembly block. Q/probes-data (M4)
219. LD to a chosen lead variant from a per-pair table, and LocusZoom bins. R/Smaller points, Q/probes-data (L1 to L3)
220. A sample-by-position matrix by region (genotypes, methylation per CpG with missing values), wide and long forms. X/Gaps, Q/probes-data (ME1 to ME4)
221. A gene-by-sample heatmap (expression) as split ranges or a matrix of ranges. Q/probes-data (E1, E2)
222. Many samples' variants: one row per sample, a long deletion solid across the SNPs it covers. R/Many samples
223. A binned contact map (Hi-C) by region and level. R/Signal and statistics, X/Gaps
224. Results keyed by id joined to features (differential transcript usage: dIF and padj per transcript), gated by a rule. R/Notes on the format itself, R/A probe: differential transcript usage

**Features and navigation**

225. Numeric filters (MAPQ at least 20), not only value lists. X/Gaps
226. Lane order on ranges (features or reads sorted). X/Gaps
227. A gene hierarchy deeper than one level (gene → transcript → exon). R/Smaller points
228. Opening on a gene by name from a file read by region (a name index). R/Smaller points
229. Queries over items read by region (flanks of streamed genes). X/Gaps
230. An annotation over a nested source must be placed or say why. X/Gaps
231. A color on a highlight entry. X/Gaps
232. Arrowheads only when asked (`direction_field`), not from the file's strand by default. P/What JBrowse could not honour
233. One strand of letters across forward and reversed blocks: reversing never complements by itself. P/same
234. A ruler that is absent or below other tracks. P/same

**Format behavior**

235. An unknown kind in a newer minor version makes the nodes that name it unavailable, not errors. R/Notes on the format itself
236. `diff(a, b)` writes a param back to its default correctly (not `null` for a non-null default) and includes snapshots and metadata. R/same
237. A session as a spec: track order, hidden tracks, highlights, bookmarks, a hub of tracks; and its patch (reorder, unhide, add). R/Examples (11)

**What another renderer may lack** (each a case the resolved view and fallbacks must carry: a reader that lacks it degrades
predictably, P/What JBrowse could not honour)

238. Sections whose axis is not a genome (alignment columns, structure chains, an item's own positions).
239. Placement through a mapping path (protein features along exons, variants on residues).
240. `follow` and scroll groups.
241. Several layers composited in one track row; `blocks` restricting a layer.
242. Size channels (glyph, thickness, radius, width, line style, opacity).
243. Collapse by an arbitrary field; split order and empty rows shown.
244. Density modes; templates and dynamic nodes; computed nodes and services; snapshots; collapsed groups and summaries;
     the text layout.

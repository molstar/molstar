# Themes and legends: what colors and sizes a mark, and the key that reads it back

## Intent

The third slice's themes: "field palettes (binned, continuous, log), rules, the residue schemes, legends, shared
domains", with the small layers (the ruler and numbering, secondary structure, variants, values' line and area, ranges'
line styles and arrowheads) and the p53 themes view; the theme sublanguage (the `layer.color` slot's grammar, legends'
place) designed first. The building blocks:
"**Theme.** A channel (color, size) and how its value is computed: `uniform`, by a field's values (a palette, binned or
continuous, a domain open to the data in view), by rules (thresholds), by a named scheme. Data, never code." A known
risk: "legends, navigators and the animation loop stay host-neutral (parts a host places, a port the host drives)". The
author's rule throughout: "what we can make configurable, we should".

A theme answers "why is this mark this color (or this size)?" as data in the spec, so the canvas, the SVG, a headless
render and another viewer all answer it the same way, and a legend can answer it back to the reader.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| The current Seq\*'s themes | a channel per layer param (`color`, sizes, `line_style`), each a theme; palettes as MolViewSpec's (categorical, discrete, continuous), evaluated by Mol\*'s `value-palette` so a value has one color in both; one `rules` evaluator, named schemes being tables for it; a measure that does not apply is NaN and fails its clause (EC-52, EC-135); a hash for values the palette does not list (EC-134, EC-203); domain groups as a union per group and field, an explicit domain winning (EC-56); sizes change the fill, never the slot or hits (EC-137); seven line styles (EC-128); a dark variant per scheme and a legend of the variant drawn (EC-136) | `ThemeRef` nesting with a fallback per theme, four levels (here one grammar, the layer type's default as the only fallback); palette order by first occurrence over all rows (EC-133, impossible by region); a group named by a free string (here a scale node with an address) |
| The current Seq\*'s legend and scheme picker | a ramp for continuous, hard steps for binned, swatches for categories; identical groups merged and titled by what they share; "+n more"; hashed values cannot be listed, said once (EC-132, EC-204); the picker in the legend group's title and the track menu, offering only schemes that fit | the legend as one viewer-owned panel, collapsed by default (here a part a host places); a scheme pick as UI state outside the spec unless asked to write (here a pick is a patch, one undo step: ux-foundations section 8, "edits to the spec (… a theme) are patches") |
| Jalview's colour schemes and "colour by annotation" | the scheme catalog (Clustal X, Zappo, Taylor, hydrophobicity, helix, strand and turn propensity, buried index, nucleotide, purine/pyrimidine, percent identity, BLOSUM62) as data with provenance; Clustal X and percent identity as threshold rules over a column; colouring by an annotation row with a threshold or a gradient between two colors; "colour text" (letters colored, not cells) | the scheme applied per selection group as a global mutable state; user-defined colours as a modal dialog (here a command per value, one patch) |
| ESPript | identical columns as a filled box, similar residues as colored letters (a cell channel and an ink channel); its black-and-white figure as a mode of the same rules (`espript :mode bw`: identical black cells, similar grey); gaps counted against identity (EC-135's flaw fixed) | its global similarity score (here an approximation by residue groups, stated); frames around similar blocks (the Frames feature, [layouts.md](layouts.md)) |
| Nightingale's tracks | a color per feature from the data (a `color` field drawn as is); a colored sequence as a heat strip of per-residue values | scales and color ranges as strings (`"A:0,B:1"`); no legend component |
| Vega-Lite's scales and legends | a scale apart from the encoding: type (linear, log), domain, range, shared by name across marks (`resolve: shared`); legends derived from scales, never written by hand; conditional encodings (our rules); a legend click selecting the items of a value | legend placement inside the chart (`orient`); expressions as strings for tests |
| Mol\*'s color themes (`graphics/src/theme/color.ts`, `core/util/color`, `core/util/legend`) | a provider registry with `isApplicable` (a scheme offered only where it fits); legends as data, `TableLegend` and `ScaleLegend` (a Seq\* legend converts to them, so a Mol\* host's legend control shows it); the color lists (ColorBrewer, Viridis, Cividis); `hashFnv32a` | a theme's params as a `ParamDefinition` per provider (here the language description, which also generates the settings form later) |
| This codebase ([ux-foundations](ux-foundations.md) section 1, [layout-and-draw-lists](layout-and-draw-lists.md)) | data and chrome never share a color role; palettes are tokens per theme; ink follows its fill; missing is not zero; conflict is one color; colors resolve through the token table, whose version keys the draw lists; `legend?(lay)` on the layer type | — |

### Interactions in the precedents

| What a user does | Where | Here |
|---|---|---|
| picks a color scheme for a layer from a menu | Jalview, Mol\*, the current Seq\* | **built**: the track menu's Color and the legend's title list the schemes that fit; each previewed while hovered; a click is one patch |
| applies one scheme to several tracks at once | Jalview ("apply to all groups") | **built**: `layer.color` takes several layers, one patch |
| colors by an annotation, with a threshold or a gradient | Jalview, IGV, Nightingale | **built**: a field theme (binned, continuous) or rules |
| sets a track's data range, or autoscales a group of tracks together | IGV, genome browsers | **built**: a scale node's `domain` (one patch, every layer naming it follows); open when unset, unioned over its layers (EC-56, EC-201) |
| switches to a log scale | IGV, Vega-Lite | **built**: `scale.type: log` on y and on colors (EC-215) |
| edits the color of one category | Jalview (user colours), plotting tools | **built**: `theme.value-color`, one patch, from a color input in the legend part |
| drags a continuous legend's ends to change the domain | some plotting tools | **later** (decided: not this slice); the command and a numeric field are built |
| moves an identity or conservation threshold with a slider | Jalview | **planned** (the column measures arrive with the rows layer, [grids-and-rows.md](grids-and-rows.md)); the threshold is a rule's number, so the slider is a preview ending in one patch |
| shades by conservation over an alignment | Jalview, ESPript | **planned** (the rows layer's column measures); the schemes are listed now and offered once they fit |
| clicks a legend entry to select its items | Vega-Lite | **built**: the items of that value, by `(item :<field> v)` |
| hovers a legend entry to see its items | plotting tools | **built**: a highlight channel over its items in view |
| clicks a legend entry to hide its items | plotting tools | **built**: `layer.filter` (one patch on the layer's `layer.data` filter, the query sublanguage's), the track header showing an "n filtered" chip whose click clears it; a theme still never hides |
| sees the active highlights in the legend and removes one | Vega-Lite's selections, some genome browsers | **planned** (the legend part's highlight rows) with `highlight.clear` |
| sees the hovered item's value on the legend | some plotting tools | **built**: the hovered item's entry is outlined, its value ticked on a ramp |
| colors the letters instead of the cells | Jalview ("colour text") | **built**: the letters layer's `ink` channel |
| colors a selected region on a track | Jalview (groups) | **planned**: `theme.color-selection` adds a rule `(in <the selection>)` at the top of the layer's theme, one patch |
| switches light and dark, or the platform's forced colors | Mol\*, the current Seq\* | **built** (light and dark: the token table); forced colors: **planned** (with accessibility) |
| asks for a color-blind safe rendering of someone else's figure | — (no precedent offers it per reader) | **built** as a reader's setting (Self-check) |
| moves the legend into the figure | Vega-Lite (`orient`), Mol\* (viewport overlay) | **no** on the canvas: never over data (ux-foundations section 9); **built** in the SVG export when asked (`legend`) |

## Decision

### The author's decisions

The author, 2026-11-18: "ok as recommended; if there are gaps in mol*'s themes, add them in non-breaking
way."; in general: "a lot of this is "configurable behavior" and where appropriate, the code should allow the
configuration via spec ports etc."; on sizes: "don't overdo it on the limits (they were sort of self-imposed and
estimate), but obviously try to minimize the amount of code". On the review of the decisions: "ok" to its recommendations,
here: a layer's default theme is fixed by the language, not the host (one spec draws the same in canvas, SVG, headless
and another viewer); a host recolors through `--seq-*` tokens.

1. Shared domains are `scale` nodes under `views`, not a group name string.
2. Palettes are evaluated by Mol\*'s value palettes, spelled in this language's words, with a converter to and from
   MolViewSpec's; `ValuePaletteLists` becomes a lazy getter in `@molstar/core/util/color`, its API unchanged.
3. Unlisted values are hashed with FNV-1a over UTF-8 (EC-134 restored; five lines of two records change), and the
   builder lists an inline column's values.
4. The `seq-cat-*` tokens are redesigned in both tables to pass the color-blind check (ΔE ≥ 0.06 in oklab, normal and
   three simulations), with item 3; every record that draws a category changes once.
5. The default continuous palette is the token ramp (`seq-ramp-*`).
6. The color-blind substitution is a reader's viewer setting, never the spec.
7. Sizes: no share of headless rendering's room; the code kept minimal (the author: "don't overdo it on the limits
   ... but obviously try to minimize the amount of code"), measured when themes land, the target moved then if needed.
8. No ramp-end drag in this slice: the command and a numeric field now.
9. The p53 themes view's per-residue values are AlphaFold DB's pLDDT for P04637, within the fixtures' licence rule ([examples.md](examples.md)).
10. A layer's default theme is the language's (uniform tokens, the token schemes); a default that would read a column
    is written by the builder as a `field` form; a field's palette is continuous unless the form says `:mid` or names
    the `diverging` scheme. A host recolors through `--seq-*`.
11. What Mol\*'s themes lacked, added without a break: `ValuePaletteLists` lazy, the optional palette fields, oklab, the
    color-vision check, the residue schemes' data, the hydrophobicity scale moved to `@molstar/core` with a re-export
    from `@molstar/model`; `okabe-ito` is core's already.
12. The author, on the remaining questions: "others as recommended": `layer.filter` writes the layer's data
    filter (one patch) and the track shows an "n filtered" chip that clears it; export choices (`legend`, `navigators`,
    `window`) are optional spec params of a view, `renderSvg`'s options overriding them.

### The writer says what is colored by what; the reader's settings say how it looks to them

- **In the spec**: each channel of a layer (which field, which palette, which scale, which rules, which scheme). Why:
  the meaning of a color is the figure's; it is saved, diffed, patched and drawn the same by every reader.
- **Not in the spec**: light or dark, forced colors, a color-blind safe substitution, the legend's place. These are
  the reader's (settings and the host's chrome); they change how a meaning is shown, never which. Why: a published,
  read-only figure (a full-state story step refuses edits) must still read for its reader.

### Themes color and size; they never hide

- **A theme never removes a mark**: hiding items by kind, score or enzyme (tRNA's tertiary pairs, AlphaFill's
  transplants under a score, a plasmid's enzymes) is a filter on the layer's data (the query sublanguage's filters in
  `layer.data`), not a theme rule; the datasets' views use data filters. Why: a mark at opacity 0 keeps its slot and
  its hit box (EC-137), so a "hidden" item would still be hovered and clicked, and a legend would still list it.

### The theme sublanguage: a slot per channel

- **Slots**: a layer's channel params are `type: 'node'` slots naming the theme sublanguage (`seq.theme`), read and
  printed in place, as queries' are ([language-description](language-description.md), Slots). A plain string in a slot
  is the `uniform` form (`"#d62728"`, `"token:range"`), so every spec written so far stays valid.
- **Channels** (each layer type declares which it reads; a theme on a channel it does not read is
  `seq.render.theme.channel`, a warning, drawn without it, EC-137):

  | Channel | Value | Read by |
  |---|---|---|
  | `color` | a color | every data layer (letters: the cell), `grid` (a cell's value), `pairs` (an arc or chord), `rows` (a row's band or label chip) |
  | `ink` | a color, or `auto` (follows its fill, 4.5:1) | `letters` (Jalview's "colour text", ESPript's colored letters) |
  | `size` | a number in the layer's unit (`band`: a share of the slot; `px`) | `ranges` (height in its lane, band), `variants` (radius, px), `values` line (thickness, px), `pairs` (thickness, px) |
  | `opacity` | 0 to 1 (never a way to hide: above) | `ranges`, `variants`, `values` area, `pairs`, `grid` |
  | `line` | `solid`, `dashed`, `dotted`, `wavy`, `double`, `box`, `highlight` (EC-128) | `ranges` |
  | `arrow` | `none`, `end`, `start`, `both` | `ranges` (only when asked, EC-232) |
  | `shape` | `circle`, `square`, `triangle`, `diamond` | `variants` (EC-242's glyph) |

- **A channel changes the mark, never its slot or its hit box** (EC-137): lanes, rows and hit boxes are laid out
  before any channel is read. Why: a size by score must not reorder lanes or move what a click finds.
- **The forms** (kinds of `seq.theme`; printed as S-expressions, JSON as nodes):

  | Form | Says | Example |
  |---|---|---|
  | `uniform` | one value | `(uniform :value "token:range")`, or the string |
  | `field` | a field's value through a scale | `(field :name label)`, `(field :name score :scale /views/conservation)` |
  | `rules` | the first rule whose clauses all hold gives its value or a nested form; else `:else` | below |
  | `scheme` | a named scheme from the registry, with overrides | `(scheme :name clustal-x)`, `(scheme :name zappo :overrides {"C" "#e0c000"})` |
  | `letters` | a value by the letter (a dictionary, with groups that carry labels) | `(letters :colors {"C" "#e0c000"} :groups [(group :letters "AILMV" :value "#ffafaf" :label aliphatic)] :default "token:missing")`: `:colors` first, then the groups, then `:default` (absent: the layer's default) |

- **The JSON shape of a slot body**: a string (`uniform`), or a node of `seq.theme` without an id,
  `{ "kind": "field", "params": { "name": "score", "scale": "/views/conservation" } }`, nested forms and palettes
  being nodes too (`"palette": { "kind": "continuous", "params": { "colors": "viridis" } }`); a scale is named by its
  address as a string.

  A scale node is named from inside a `field` form, never as a slot's whole value; a named theme reused across layers
  waits for an example ([composition](composition.md), "Definitions and uses wait"), since a scale node already
  shares what two tracks need to share.
- **`field`**: `:name` (a column of the layer's data; `value` for a values layer's own values; `letter` for the cell's
  letter), `:data` (another data node on the same system or through a path, read per position: residues colored by
  a values node such as a confidence track), and the scale: inline (`:palette`, `:domain`, `:type`, `:mid`) or a
  `scale` node by address. A field with no palette takes, by its column's declared type (never its values), the default
  categorical palette (text) or the continuous ramp (`seq-ramp-*`, numbers); a palette diverges only when the form says
  `:mid` or names the `diverging` scheme. Why plain: a default read off the data in view would flip a track's palette
  while panning, and differ between two readers of one spec.
- **`rules`**:

  ```
  (rules
    (rule :when [(>= (field plddt) 90)] :value "#0053d6")
    (rule :when [(>= (field plddt) 70)] :value "#65cbf3")
    (rule :when [(in /refs/dbd)] :value (field :name score))
    :else "token:missing")
  ```

  Clauses: `=`, `!=`, `<`, `<=`, `>`, `>=`, `one-of`, `in` (a location: an element whose position or extent overlaps
  it by one position or more) and `within` (an extent contained in it). Measures: `(field <name>)`, `value`, `letter`,
  and the column measures (`share`, `identity`, `matches-reference`), which are the rows layer's ([grids-and-rows.md](grids-and-rows.md)): built with its
  column summaries and registered as measures for themes. A measure that does not apply is NaN and fails its clause
  (EC-52). **Types**: `=`, `!=` and `one-of` compare text with text (case-sensitive) and numbers with numbers; `<` to
  `>=` take numbers only; a comparison across types is false, and `seq.lang.theme.type` when the check can see it (a
  literal against a column's declared type). A rule's value is a channel value or a nested form, nesting bounded
  (`SpecBounds.themeDepth`, 4). Several clauses are an `and`; an `or` is two rules. No `:else` gives the layer type's
  default. Why rules, not expressions: data that validates, diffs and prints, as Vega-Lite's conditions without its
  expression strings.

### Palettes and scales

- **Three palette kinds**, evaluated by one implementation with Mol\*'s value palettes (`@molstar/core/util/color`),
  spelled in this language's words:
  - `categorical`: `:colors` (a list, a named list, or a dictionary value → color), `:values` (the listed values, in
    order), `:repeat` (true), `:unlisted` (`hash`), `:others` (the list unlisted values hash into; the dictionary's
    palette default, `categories`);
  - `direct`: the field already holds a color (Nightingale's `color` attribute, a file's `itemRgb`): read as a spec
    color; one that does not read is missing, said once (`seq.render.theme.color`);
  - `binned` (Mol\*'s and MolViewSpec's `discrete`): `:edges` (n+1 numbers, a bin closed below and open above, the last
    closed) and `:colors` (n), or `:bins n` over the domain; `:under`, `:over`;
  - `continuous`: `:colors` (a named list or stops `[color, at]`, `at` in 0 to 1 of the domain; `:stops absolute`
    reads `at` as values, MolViewSpec's `mode`), interpolated in oklab (fixed), `:mid` (a diverging center, a number
    the writer gives), `:under`, `:over` (the ends' colors by default).
- **A value the palette does not list takes one color everywhere** (EC-203): the color at the 32-bit FNV-1a hash
  (`hashFnv32a`) of its UTF-8 bytes, modulo the palette's list (a dictionary palette hashes into its `:others`). This
  restores EC-134's rule (FNV-1a over UTF-8); EC-134's "a dictionary gives the missing color" is changed by EC-203.
  Why fixed and not a setting: another reader must reproduce it from the format's documentation alone. `hashFnv32a`
  takes bytes: a value is encoded once per dictionary code of its column (or once per distinct string in a tile),
  never per mark. Today's `categoryToken` (a UTF-16 string hash, `render/tokens.ts`, used by `render/ranges.ts` and
  exported from `render/index.ts`) is replaced, so today's category colors change once: the committed records
  `render/_test/records/p53-gene.txt` (lines 23, 118, 119) and `p53-tracks.txt` (lines 40, 41) are rewritten in the
  same change. The viewer does not hash: it reads `seq-cat-5` as `highlight.color`'s default (`viewer.ts`), which a
  token redesign (the author's decisions, 4) moves too. A legend cannot list hashed values it has not painted, and says so once (`seq.render.legend.unlisted`,
  info, EC-204).
- **The builder lists what it can see**: `theme.field('label')` on a layer over inline data writes `:values` from that
  column in first-occurrence order, so a small inline table gets distinct colors and the spec says which (EC-133,
  changed: the order is written, not inferred over all rows). Data read by region is never listed by the builder.
- **A `scale` node** (a new kind under `views`, id required, `/views/<id>`), as Vega-Lite's named scales: `type`
  (`linear`, `log`, `sqrt`), `domain` (`[a, b]`, either end `null`: open, EC-202), `palette` (color scales) or `range`
  (size scales, `[from, to]` in the channel's unit), `mid`, `nice` (legend ticks only), `label`, `unit`. An inline
  scale in a `field` form takes the same params (`:range [2 8]`); a radius scale defaults to `sqrt`
  (fixed for a radius unless the scale says `type`), so a mark's area follows its value.
  - **One scale, two tracks**: layers whose themes name one scale share its domain and palette; an open end is the
    union over those layers' data in view, computed when the view settles and held while it moves (EC-56, EC-201). An
    explicit end wins. **An export** unions an open end over what it draws (every line of a text, a ring whole), not
    over a viewer's view, and records the domain it used (`data-seq-scale-*`); **a navigator's miniature** uses the
    view's resolved domain, never its own, so its colors match the view's. Why a node, not a group string: it has an address, so one `set` changes every layer that uses
    it (one undo step), the legend is one per scale (merged by construction), and a view, a story step or a link can
    name it.
  - **A values layer's y** names a scale too (its `scale` param; its `domain` param stays as the inline form). Its
    color may name the same scale or another: a palette's domain is separate from the y domain (EC-205).
  - **Log** (EC-215): a value at or below zero is missing on a log scale, counted once per layer
    (`seq.render.scale.log`, info; the scale's `non-positive`: `missing`, or `clamp` to the smallest positive in view).
  - Values are drawn at their value; the domain bounds the height and the palette, never the value (EC-55).
- **Missing is not zero**: NaN, an absent field, an unmapped position take the missing token, distinct from any ramp's
  low end (ux-foundations section 1; a test checks the distance in oklab for every built-in ramp, both themes).

### Schemes: named, curated, data

- A **scheme** is a registry entry ([extensibility](extensibility.md): "themes, palettes, residue schemes"): a name, a
  label, a subject (letters of given alphabets, numbers, categories), what it needs (a measure: a column's share, a
  reference row), its light colors and its dark variant (given, or derived: each color's oklab lightness mapped into the
  dark theme's band, hue kept, chroma clamped to the gamut), whether it is color-blind safe, its safe substitute if not,
  its legend labels, and its provenance. It expands to a `letters`, `field` or `rules` form. Built-ins are entries of
  the same kind, in the default composition of the viewer and of `createSeqContext`, so the canvas and the export draw
  the same scheme.
- **Built in** (provenance kept per entry, as the current Seq\*'s table: published colors and rules, never GPL code):

  | Subject | Schemes |
  |---|---|
  | protein letters | `classes` (the default: the `seq-aa-*` tokens), `classes-safe` (the safe substitute), `clustal-x`, `zappo`, `taylor`, `hydrophobicity`, `helix`, `strand`, `turn`, `buried`, `percent-identity`, `blosum62`, `espript` |
  | nucleotide letters | `nucleotide-classes` (the default: the `seq-nt-*` tokens), `nucleotide`, `purine-pyrimidine`, `nucleotide-cb` (Okabe–Ito) |
  | numbers | `ramp` (the default: the `seq-ramp-*` tokens), `diverging` (the `seq-div-*` tokens), `viridis`, `cividis`, `magma`, `blues`, `rd-bu`, `plddt` (AlphaFold's four bins, binned) |
  | categories | `categories` (the default: the `seq-cat-*` tokens), `okabe-ito` (core's `lists.ts` has it), `tableau10`, `category10` |
  | secondary structure | `secondary` (the default: `seq-ss-helix`, `seq-ss-strand`, `seq-ss-turn`) |
  | chromosome bands | `giemsa` (the `bands` layer's default, the navigators design: a `stain` field's values as grays derived from `msui-text`, so dark flips them; `acen` in `seq-centromere`) |

- **A scheme applies only where it fits**: the picker offers only fitting ones; a spec naming one that does not fit (an
  alphabet it does not serve, a measure not available: Clustal X, percent identity, BLOSUM62 and ESPript before the rows
  layer's column measures) draws the layer's default and reports `seq.render.theme.scheme` with what it needs.
- **ESPript, fixed** (EC-135): identical where every row of the column has the letter, gaps counted against it
  (the scheme's `:gaps count`; `ignore` the old behavior); similar where the letter's group holds at least `:similarity` (0.7) of the column; cells
  filled for identical, letters inked for similar. Stated as an approximation of ESPript's score by groups.
  `:mode bw` is ESPript's black-and-white figure: identical cells black (ink white by the ink rule), similar cells grey
  (`#c8c8c8`), nothing elsewhere; it reads the same in print and needs no hue.
- **Defaults are to be color-blind safe; today's are not.** The target: every default passes the conformance check in
  both themes. The check: any two entries of a categorical palette at least 0.06 apart in oklab (ΔE) under normal
  vision and under simulated protanopia, deuteranopia and tritanopia (Machado et al. 2009, severity 1); a ramp's
  adjacent stops at least 0.03 apart in oklab lightness, monotone. Today's `seq-cat-*` (light table) fail it: the
  smallest distance is 0.032 under normal vision (`#9a7aa6` and `#a67a95`), 0.016 under deuteranopia, 0.023 under
  protanopia. So the category tokens are redesigned for the check (the author's decisions, 4), in both tables; every committed record that
  draws a category changes once, with the hash's change above. A scheme that is not safe (Zappo, Taylor, Clustal X) is
  offered with that said in the picker, never chosen by default.

### Colors: literals, tokens, light and dark, forced colors

- **A color in a spec** is a hex literal, a token (`token:<name>`), or `light-dark(<a>, <b>)` (the CSS function's
  meaning: the first in the light theme, the second in the dark). A plain literal is drawn as written in both themes:
  the writer's color. Tokens, schemes and the default palettes follow the token table, so a theme switch redraws them
  and only them. Why `light-dark`: a figure's own color that must read in both themes needs both values, and the CSS
  spelling is one a reader already knows.
- **New tokens** (light and dark tables, and `--seq-*` a host may set): `seq-ramp-0` … `seq-ramp-4` (a ramp from near
  the track's color to strong, monotone in lightness), `seq-div-neg`, `seq-div-mid`, `seq-div-pos`, `seq-ss-helix`,
  `seq-ss-strand`, `seq-ss-turn`, `seq-variant`, `seq-ruler`; chrome the other designs draw: `seq-break` (the axis
  design's break ink), `seq-window` (a window's fill and edge: rings and navigators), `seq-grip` (its grips),
  `seq-wash` (outside a navigator's window and its unreachable span), `seq-centromere` (a hue, 3:1 against the surface in
  both themes). The `seq-density-*` pair stays the density strip's. A host recolors any of them through `--seq-*`.
- **Ink follows its fill** for letters and labels on any colored mark (4.5:1). The tables' inks are not pure black and
  white, so a mid-tone fill can leave both under 4.5:1 (a saturated scheme color between them): there the ink is pure
  black or white, whichever is higher (the better of the two is always above 4.5:1). An `ink` theme's color under 3:1
  against its fill is moved in oklab lightness until it reaches 3:1 (fixed: a theme must not cost the letters, ux-foundations section 1). Tested over every color of every built-in scheme, both themes, not only the tokens.
- **The hover outline holds 3:1 against any fill** (EC-138): where `seq-hover` falls under 3:1 against the mark it
  outlines (an arbitrary scheme color), a 1 px halo in the opposite ink is drawn outside it; tested over every scheme
  color. Exports draw no hover.
- **Compatibility**: a string in `layer.color` reads as before, and `token:<name>` keeps its three-step lookup
  (`seq-<name>`, `msui-<name>`, `<name>`). A non-string value (accepted today by the `any` param and reported at draw
  as `seq.render.color`, a warning) becomes a refusal at the check (`seq.lang.theme.form`) unless it is a form; a
  string that is no color is `seq.lang.theme.color` at the check, and `seq.render.color` stays only for a token the
  table lacks. `light-dark()` extends `resolveColor` (`render/tokens.ts`, which picks by the table's `theme`), and
  `ref.color` (documented as hex or token) and the viewer's `highlight.color` items take it too.
- **A reader's color-blind substitution** (`theme.cvdSafe`, a viewer setting, off): each scheme that is not safe is
  drawn as its declared substitute, and a categorical palette as `okabe-ito`; the legend shows what is drawn and says
  it is substituted. Never written into the spec. **Every substitute is itself safe** (decided, 2026-11-19): a scheme
  may not name an unsafe scheme as its substitute (a test over the registry). The protein letter schemes that are not
  safe (`classes`, `clustal-x`, `zappo`, `taylor`) substitute to `classes-safe`, `classes`' groups on Okabe–Ito hues
  (passing the conformance check in both themes). `classes` stays the default (decided: the reader's setting covers
  it; changing the default changes every figure).
- **Letter cells in dark** (decided, 2026-11-19): the dark scheme variants sit at oklab L 0.32–0.61, 1.4–1.7:1
  against the dark page, and cannot go lighter without their ink failing 4.5:1. So in the dark theme a letter cell
  draws a hairline outline (1 device px, `seq-cell-edge`, a new token: 3:1 against the dark page; none in light),
  rather than a lighter page for every figure.
- **Forced colors** (planned, with accessibility): the viewer reads `forced-colors: active` and builds its token table's chrome and
  emphasis tokens from system colors (`Canvas`, `CanvasText`, `Highlight`, `GrayText`), so the selection, the hover and
  the cursor stay visible; data colors are kept (fixed), since the data's meaning is in them; the
  legend part's swatches opt out of forcing (`forced-color-adjust: none`) so they match the canvas. An export never
  uses forced colors.

### Legends: data in the core, a part a host places

- **A legend is derived, never written**: from a layer's resolved channel (its theme, its scale, its domain as now
  resolved) the core makes `Legend` data: `title`, `channel`, the scale's address when shared, `entries` (a `ramp` with
  stops, domain and ticks; `steps` with edges; `swatches` with labels and colors; `sizes` with marks at chosen values),
  whether the domain is open ("in view", with its numbers), `source` (`vocabulary`: every listed value; `painted`: what
  was drawn; EC-204 says which), `more` (the count past the part's `maxItems`, 12), notes (hashed values unlisted). It replaces
  today's `LegendSpec` from `legend?(lay)`.
- **Merged**: legends with the same key (the theme's hash, the scale, the domain, the entries) are one, titled by what
  their layers' labels share, else by the scale's label (EC-132).
- **Published on change only**: the viewer's `legends` state is derived when the view settles, a theme, scale or token
  table changes, or the hover moves onto another entry; never per frame, compared structurally.
- **A part a host places** ([ux-foundations](ux-foundations.md) section 9: "Legends sit with their track or in one legend
  panel, never over data"): `Legend` (`@molstar/ui` components and tokens) takes a list of legends and renders it
  wherever a host mounts it (`mountLegend(element, viewer, { of })`, `of` a track's or a scale's address to place
  that one's legend, every legend when absent): beside a track header, in a panel, in a Mol\* plugin's own
  legend control (a Seq\* legend converts to `@molstar/core/util/legend`'s `TableLegend` or `ScaleLegend`). The default
  mount is a `SeqViewerSpec.ui` entry (`ui.legends`, the same call), its gestures a `BehaviorProvider` (`legend`), as
  the navigators' are. Seq\* never draws a legend on its canvas.
- **The legend follows the hover**: the hovered item's entry is outlined; on a ramp, its value is ticked.
- **Highlight channels in the legend part**: it lists the active highlight channels (a link's, a legend hover's, a
  host's) as rows with their color and a clear, which runs `highlight.clear` (a new command, beside the
  selection's bookmark command).
- **Every entry is reached from the keyboard** (ux-foundations section 8: every gesture has a keyboard way): the part
  is one tab stop, a list with a roving focus; arrows move between entries, Enter or Space selects (as a click),
  Cmd/Ctrl+Enter toggles, the context-menu key opens its actions (the color input, the scheme picker); focus on an
  entry hovers its items as the pointer does. Each entry's accessible name says its value and title ("Domains:
  DNA-binding"); a ramp says its domain and whether it is open.

### Editing a theme: one patch, previewed

- A theme form inside a slot has no node ids (a slot's body is replaced whole, [patches](patches.md)); an edit `set`s
  the channel param or a scale node's param: one patch, one undo step. Why whole: forms are small, and ids inside them
  would be addresses nobody needs.
- **Commands** (`seqx-viewer`, editing): `layer.color` (a form, one layer or several: one patch), `layer.channel` (any
  channel), `scale.set` (domain, type, palette of a scale node), `theme.value-color` (one value's color in the layer's
  palette or its scale's), `theme.color-selection` (planned), `legend.select`, `legend.hover`. A pick hovered in a
  menu is a preview; a click ends it in one step; Esc cancels ([changes](changes.md), previews).
- **On a default palette**, `theme.value-color` writes the layer's theme out: the layer's default becomes an explicit
  `field` form over the same field with a `categorical` palette whose dictionary holds that one value and whose
  `:others` is the default palette, so every other value keeps its color. On a scale's palette it writes the scale.
- **Back to the default**: `layer.color` with no form unsets the param (the layer's default again); `theme.value-color`
  with no color removes the value's entry; `scale.set` with a `null` end reopens it. A layer leaves a shared scale by
  a `layer.channel` that names another scale or an inline one; the scale node stays, drawn by no one and listed by no
  legend until a layer names it again or it is removed.
- **Legend selection follows the grammar's click rule**: clicking the only selected value clears the selection.
- **Story steps**: a theme change is a patch like any other; between two steps a color channel blends in
  `story.color.space` ([stories](stories.md)).

### What the SVG export draws

- **The canvas's picture**: the same draw lists, every token resolved with the export's table and written as a literal
  (`light` headless unless asked; the viewer's shown theme for `viewer.renderSvg`), never forced colors or a reader's
  color-blind substitution unless asked (`cvdSafe` as an export option).
- **Legends off the canvas, in the export only when asked** (settled by ux-foundations section 9 and EC-138: legends
  never over data, no hover in an export; never drawn with a hover outline)
  (`renderSvg({ legend: 'none' | 'below' | 'right' })`, `none` by default): laid out by the
  core's `layoutLegends(legends, box, tokens, text)` into the `labels` pass in a band outside the data box (a ramp as
  one row of the `image` primitive, as density strips are), so the export is a figure without host chrome. The SVG
  records each scale's resolved domain (`data-seq-scale-*`), the union over what the export draws (above).
  **Export choices are the view's** (decided): `legend` (and the navigators' and the ring's window, in their designs)
  are optional spec params of a view (`export: { legend, navigators, window }`), so a figure says how it exports;
  `renderSvg`'s options override them for one call.

### The small layers: what each reads from the theme

| Layer | Reads | Default | Notes |
|---|---|---|---|
| `ruler` | tokens only (`seq-ruler` ticks, `msui-text-muted` labels); a color theme on it is `seq.render.theme.channel` | — | `numbering` (spec param: the numbering it labels by, the section's default numbering when unset; switching schemes, Kabat against IMGT, is a patch of it); insertion labels (100A, 111.1) as minor ticks, labeled only with room; tabular figures, steps of 1, 2, 5 × 10ⁿ, a label skipped when it does not fit (a gap of 48 px, a foundations number in the presets), never shrunk; above, below or absent is its track's place and existence (EC-234); fallback: omitted |
| `secondary` | `color` (by the `kind` field), `size` (band) | `(scheme :name secondary)` | helices and strands as blocks (the layer's `style`: `blocks`; `cartoon` with the text figure, [layouts.md](layouts.md)); the strand's arrowhead on the piece holding its C-terminal end, whatever the direction or the origin (EC-130, fixed); cut ends square; 200 copies draw what one does (EC-67); fallback: `ranges` |
| `variants` | `color`, `size` (radius), `opacity`, `shape` | `token:variant` (the builder writes `(field :name consequence)`) | a stem and a head a site; at `bars` and `density`, a count strip in the density tokens; fallback: `ranges` |
| `values` | `color` (uniform, or a field over `value`), `size` (line thickness), `opacity` (area) | `token:value` | `encoding`: `bars` (built), `line`, `area`; a line changes color halfway between positions and holds over half cells at the ends and at gaps (EC-57); missing breaks the line; the area under the line at the `opacity` channel (0.35 by default); y through its scale |
| `ranges` | `color`, `size` (band), `opacity`, `line`, `arrow` | `token:range` (the builder writes `(field :name category)` where it has the column) | a line style anchored to the content (the dash phase moves with the data, not the screen); an invalid style drawn solid and said once (`seq.render.line-style`); an arrowhead on the piece holding the end it points to, never on a cut end |
| `letters` | `color` (cells), `ink` | `(scheme :name classes)` or `nucleotide-classes` by the system's declared alphabet | ink `auto` follows the fill |
| `bands` (the navigators design) | `color` (by the `stain` field) | `(scheme :name giemsa)` | the `acen` pair as triangles; fallback: `ranges` |

Every default above is fixed by the language (the layer type's description), never by the host or the data: one spec
draws the same in the canvas, the SVG, headless and another viewer (decided). Where a default depends on a column, the
builder writes the `field` form into the spec instead.

### Performance

- **Compiled once, looked up per mark**: a form compiles in tier 1 to an evaluator, cached by its hash; per frame a
  continuous palette is a lookup table (256 colors, fixed) rebuilt only when the tokens' version or the
  domain changes; a categorical field over a dictionary-coded column is a color per dictionary code, made once per tile.
  Why: the layout's hot loop reads an index, never a string or an allocation.
- **Targets** ([performance](performance.md); each names its path and test, a quick representative where it is built):
  - the layout with a field theme at a fixed screen costs at most 1.1× the same layout uniform (10,000 ranges by a
    category; 300,000 values on a continuous log scale), `itWithin` in the heavy lane;
  - draw and hover warm: 0 bytes with themes (the opt-in allocation tier);
  - legends derived in under 1 ms for 100 tracks, on settle only (not per frame).
- **Size**: the author, 2026-11-18: "don't overdo it on the limits (they were sort of self-imposed and estimate), but
  obviously try to minimize the amount of code". So a size estimate is never a reason to cut a feature; the code is kept
  minimal instead: Mol\*'s color code reused, no private copies, nothing without a use. Estimates, minified, to be
  measured:

  | Piece | kB |
  |---|---|
  | palette evaluation (Mol\*'s value palettes without the lists) | 5 |
  | the color lists the built-in schemes name (seven) | 2 |
  | the scheme registry, expansion and residue tables | 5 |
  | oklab conversion and blending (core has `lab`, not oklab: added there) | 1 |
  | the LUT and the evaluator compile (forms, rules, measures) | 4 |
  | the `seq.theme` description and the `scale` kind | 3 |
  | legends (data, merge, `layoutLegends`) | 3 |
  | the small layers (ruler 2, secondary 2.5, variants 3, values line and area 2, ranges' line styles and arrowheads 2.5) | 12 |
  | **total** | **about 35** |

  The headless target is measured when themes land and moved then if the slice needs it. The lists:
  `value-palette.ts` builds `ValuePaletteLists` from every list at module load (about 48 kB of source with `lists.ts`).
  Decided: make it a lazy getter in `@molstar/core/util/color` (the public API kept; MolViewSpec reads it through
  `helpers/colors.ts`), with Seq\* passing explicit colors to the evaluation and naming lists only through its schemes
  module; measured when themes land. The lazy getter saved work but not bytes (`value-palette.ts` imports `lists.ts`
  statically: about 13.6 kB in headless rendering and the viewer), so (decided, 2026-11-19) the evaluation moves to a
  lists-free module in `@molstar/core/util/color` that Seq\* imports, and `value-palette.ts` re-exports it with
  `ValuePaletteLists`: non-breaking for Mol\*'s and MolViewSpec's code, no entry in `breaking-v6-changes.md`.

## The interaction matrix's cells

This design fills the **track or layer** row's edit cell, the **selection** row's edit cell in part, and adds a row,
**a legend entry** (a swatch, a step, a ramp), to [ux-foundations](ux-foundations.md) section 8.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a track or layer** | ([grids-and-rows.md](grids-and-rows.md)) | — | — | ([grids-and-rows.md](grids-and-rows.md)) | — | — | — | `track.hide`, `layer.hidden` (built) | its color and channels from the track menu (`layer.color`, `layer.channel`; each scheme previewed while hovered, one patch); several layers at once in one patch; a shared scale's domain, type or palette (`scale.set`: every layer naming it follows, one step) |
| **a legend entry** | click, or Enter on a focused entry: its items (`legend.select`, the selection becomes `<data>/(item :<field> v)`; the only selected value clicked again clears); a step or a ramp: no, until the query sublanguage filters by value | no: a field's values have no order for categories; steps wait with the value filter | Cmd/Ctrl+click or Cmd/Ctrl+Enter toggles the value's items in or out | no gesture of its own: a value's items can span the whole data; Enter focuses the selection as for any selection | Alt+Left (the view's) | no: a legend has no axis | no: a legend has no axis | the value's items filtered out: `layer.filter` from the entry's menu (one patch; the track's "n filtered" chip clears it); the part is the host's to hide | commands only, the controls the application's (the author, 2026-11-19: changing colors and ranges from a legend "is an application concern"): its color (`theme.value-color`, one patch), a ramp's domain (`scale.set`), the scheme (`layer.channel`), each previewed; the part shows no editing control |
| **the selection** | — | — | — | — | — | — | — | — | color it on a layer: planned (`theme.color-selection`, a rule `(in …)` added at the top, one patch) |
| **a bookmark or link** | — | — | — | — | — | — | — | — | its color: a color as above (`light-dark` allowed), by patch with its label (planned with the label) |
| **the view** | — | — | — | — | — | — | — | — | no theme edit: light, dark, forced colors and the color-blind substitution are the reader's settings, never the spec |

Hover, not a column: hovering a legend entry highlights its items in view (the `legend` channel); hovering an item
outlines its legend entry.

## Configuration: where each choice lives

The author, on this design: "a lot of this is "configurable behavior" and where appropriate, the code should allow the
configuration via spec ports etc.; we already have a lot of patterns for this, so use your own judgement what's
appropriate to keep just as a single feature and what to make configurable". The rule applied (the same in every
design of the axis, themes, layouts, grids and rows, navigators): a choice about what a figure means is a **spec param** (the writer's), its default plain, never
inferred from data; a drawing number is a **`RenderConfig` preset**; a size a spec may ask is a **`SpecBounds`** bound;
how a viewer or its reader behaves is a **`ViewerConfig`** item (the host's at composition, a reader's settings layer
later, [extensibility](extensibility.md)); a named provider is a **registry entry**; one call's choice is an **option**
of that call (`renderSvg`'s for every export choice); a mounted part's look is its **prop**; anything with no second
real use is **fixed**, with its reason. **No new port**: themes need no provider that is one per viewer
([ports](ports.md): named providers are registry entries); tokens keep their existing source (`theme.source`, CSS).

| Choice | Whose | Default | Why there |
|---|---|---|---|
| a channel's form, palette, scale, domain, rules | spec param (`layer` channels, `scale` nodes) | the layer type's default | the figure's meaning |
| a scale's `type`, `mid`, `non-positive` | spec param (scale, or inline) | `linear`; a radius `sqrt`; none (continuous); `missing` | the writer's reading of the numbers (EC-215) |
| ESPript's `:gaps`, `:similarity`, `:mode` | spec param (the scheme's) | `count`, 0.7, `color` | what the figure claims (EC-135) |
| `encoding` of values; `style` of secondary | spec param (layer) | `bars`; `blocks` | the figure's form |
| line style, arrowheads, area opacity | spec param (`line`, `arrow`, `opacity` channels) | `solid`, `none` (EC-232), 0.35 | channels like any other |
| the theme a layer takes when its spec gives none | fixed by the language (the layer type's description) | uniform tokens and the token schemes (`categories`, `ramp`, `classes`, `nucleotide-classes`, `secondary`, `giemsa`) | one spec draws the same everywhere (decided); a host recolors through `--seq-*` tokens; a column-driven theme is written by the builder |
| the reader's color-blind substitution | `ViewerConfig`, a reader's setting (`theme.cvdSafe`) | off | the reader's, never the spec (decided) |
| schemes, measures | registry entries | the built-ins | named providers, extensible ([extensibility](extensibility.md)) |
| nesting depth | `SpecBounds` (`themeDepth`) | 4 | the bounds pattern ([security](security.md)) |
| the ruler's label gap | `RenderConfig` preset | 48 px | the presets pattern (ux-foundations section 0) |
| legends in the export | spec param of the view (`export.legend`, optional), `renderSvg`'s `legend` overriding it | `none` | how the figure exports is the figure's; one call may differ (decided) |
| the color-blind substitution in an export | option of `renderSvg` (`cvdSafe`) | off | the reader's, never the spec |
| hiding a layer's items by value | command `layer.filter`, a patch of `layer.data` | — | a filter hides, a theme never does (decided) |
| the ruler's numbering | spec param (`ruler` layer's `numbering`) | the section's default numbering | which labels the figure shows |
| the legend part's default mount | `SeqViewerSpec.ui` entry (`legends`) | on in the standalone viewer | the host places parts |
| a legend's item count | the `Legend` part's prop (`maxItems`) | 12 | the host places and sizes the part |
| the builder's listing of inline values | option of `theme.field` (`list`) | true | one call's choice |
| interpolation space | fixed: oklab | — | perceptual steps; no second use (the stories' blend has its own item) |
| the lookup table | fixed: 256 colors | — | indistinguishable steps at any width; no second use |
| the hash of unlisted values | fixed: FNV-1a over UTF-8 | — | other readers must reproduce it (EC-203) |
| a scheme's dark variant when it gives none | fixed: derived | — | a scheme that wants another gives it as data |
| offering only fitting schemes | fixed | — | a scheme that cannot apply is noise (EC-136) |
| ink contrast, the hover halo | fixed | — | accessibility rules (ux-foundations section 1, EC-138) |
| legend source, merging | fixed: vocabulary when listed, else painted, always said; merged | — | the truthful key (EC-132, EC-204) |
| data colors under forced colors | fixed: kept (accessibility work may bring evidence) | — | the data's meaning is in them; no second use for `monochrome` |

## Tensions

- **One implementation, two spellings**: the palettes evaluate with Mol\*'s value palettes (a value has one color in
  MolViewSpec and Seq\*) but are spelled in this language's words (`binned` for `discrete`, kebab-case names).
  Decided as recommended: shared evaluation, our spelling, a converter both ways in `seqx-core`.
- **Hash collisions**: eight default categories and a hash give two unlisted values one color often. The builder's
  listing cures small inline data; data read by region keeps the risk, and EC-203 (one color everywhere) is worth it.
  A legend shows a collision as two labels on one swatch. Standing.
- **Literals in both themes**: a plain literal chosen on light may vanish on dark. Recommendation: `light-dark()` and
  tokens for figures meant for both; no automatic change of a writer's literal. Standing.
- **An open shared scale is per viewer**: two viewers of one spec each union what they show, so one value can differ
  in color between them. Recommendation: accepted; a figure fixes its domain; the legend says "in view" with the
  numbers. Standing.
- **Data-in-view pictures differ**: an open domain, a legend's painted entries and a density step depend on what is
  in view, so the screen, a navigator's miniature (it takes the view's domain) and an export (it unions over what it
  draws) can show one value in different colors. Standing: a figure that must match fixes its domains.
- **Data colors under forced colors**: kept colors serve meaning, forced ones serve the platform's contrast promise.
  Recommendation: keep data, force chrome and emphasis. Open until the accessibility work's look.
- **Column schemes arrive late**: Clustal X, percent identity, BLOSUM62 and ESPript need the rows layer's column
  measures, built after themes; the first themes view cannot show them. Recommendation: themes build rules and the
  registry of measures, the rows layer registers the column ones and their tests; the Globins views carry them. Accepted.
- **Meaning only in color**: a greyscale print or a reader who sees no hue loses a categorical key. Line styles, shapes
  and the ESPript's cell-versus-ink split carry some; patterns (hatching) are not designed. Standing (Self-check).

## Self-check

The claim: **this design covers what a user does with themes and legends**.

1. **Assume it holds.** Every channel of every layer this slice draws is a form in a slot; every form evaluates the
   same in the canvas, the SVG and another reader; every color has a light and a dark value; every legend comes from a
   channel; every edit is one patch.
2. **What it rests on, negated.** It rests on *the writer chooses the colors and the reader reads them*: themes in the
   spec, legends as a key. Negated: **the reader chooses, and the legend is acted on**: a reader recolors a figure they
   cannot edit; a reader asks the legend "which are these?" and "make this one different"; a reader cannot see the
   colors at all.
3. **The negated world.** A published figure or a full-state story step refuses edits, so a reader with a color vision
   deficiency has no way to make Zappo readable: the recolor must be a reader's setting, outside the spec, like dark
   and light. The legend is the inverse of the theme (value → color, read back as color → items): a click on a swatch
   is a selection by value, a hover a highlight, a color input an edit; a ramp's ends are a domain. The plural: two
   layers sharing one scale must change together (one step), and two legends that say the same must be one. The
   missing view: the export, a figure with no host chrome, needs its key drawn in it. And a reader who sees no hue at
   all (greyscale print, achromatopsia) gets no categorical key from color. The inverses of an edit: an override reset
   to its default, a layer leaving a shared scale, a selection by legend clicked again. The symmetric case of "color
   shows": "color hides" (a rule at opacity 0), which keeps hit boxes and so is no hiding at all.
4. **Both, held.** *Integrated* into the Decision: the writer–reader split (a rule of its own); the color-blind
   substitution as a reader's setting with a substitute declared per scheme; the legend row of the matrix (select,
   toggle, hover, edit color and domain) and the legend following the hover; scales as nodes so the plural edit is one
   step and the legend one per scale; legends in the export when asked; `light-dark()` for a writer's literal in both
   themes; reset, leaving a scale and the clearing click (Editing a theme); themes never hide, filters do; the legend
   reached by keyboard; ESPript's black-and-white mode as a figure needing no hue; and from the reader's side of
   "another viewer": a layer's default fixed by the language, never by the host or the data in view, so two readers of
   one spec draw it alike (the host recolors through tokens only). *Planned*: forced colors (with accessibility); column-measure schemes and the threshold slider (with the rows
   layer); coloring the selection (a command). *Standing tensions*: meaning only in color (patterns not designed), hash collisions in data read by
   region, an open shared scale differing between viewers.

## Consequences

- **`seqx-core`**:
  - `lang/`: the `seq.theme` sublanguage (the forms, the palettes, the clauses), described as the query sublanguage is,
    the channel slots on `layer` (`color`, `ink`, `size`, `opacity`, `line`, `arrow`, `shape`), the `scale` kind under
    `views`, `values.encoding` and `scale`, the new layer types' kinds (`ruler`, `secondary`, `variants`); the builder's
    `theme.*` and `b.views.scale`; the language version `0.4.0` (a string color stays valid). Issues, cataloged:
    `seq.lang.theme.form` (an unknown form or param, a non-string non-form value), `.color` (a string that is no
    color), `.channel` (a value not of the channel's kind: a color in `size`), `.type` (a comparison across types),
    `.stops` (stops or edges out of order), `.scale` (an address naming no scale), `.depth` (past `themeDepth`);
    `seq.render.theme.channel`, `.scheme`, `.measure`, `.color` (a `direct` value that does not read),
    `seq.render.scale.log`, `seq.render.line-style`, `seq.render.legend.unlisted`; `seq.render.color` narrowed to a
    token the table lacks.
  - `theme/` (an entry, `@molstar/seqx-core/theme`, so another reader evaluates a theme without the layout): compile a
    form to an evaluator; palettes over Mol\*'s value palettes; the FNV-1a hash; the scheme registry kind and the
    built-in schemes (their own module, measured); measures as registry entries; `Legend` data, merging, the converter
    to Mol\*'s legend types and to and from MolViewSpec's palettes.
  - `resolve`, `view`: themes in the resolved view as normalized forms with scale nodes followed (renderer-independent:
    the Nightingale adapter reads them); shared open domains unioned at settle.
  - `render/`: channels read by each layer type; LUTs per tokens' version and domain; the new tokens in both tables;
    `layoutLegends` and the export's `legend` option; the small layers.
  - Tests: property tests (a continuous palette monotone in its ramp's lightness; log monotone; bins closed below; the
    hash's known vectors, stable across tiles and readers; a dictionary code's color equal to its value's); a
    conformance suite over every built-in scheme and both token tables (ink 4.5:1 on every cell color, missing apart
    from every ramp's low end, color-blind distinctness for the safe ones, a dark variant for each); drawing checked
    against the data (each mark's color is its theme's color for its datum, each painted legend entry has a mark, each
    drawn color is in the legend); the size channel never moving a hit box; the allocation tier and the 1.1× bound.
- **`@molstar/core`, non-breaking** (the author: "if there are gaps in mol*'s themes, add them in non-breaking way";
  new exports and optional fields only, each with its test, so Mol\*'s own themes and MolViewSpec can use them; nothing
  of it kept privately in Seq\*):
  - `util/color/value-palette.ts`: `ValuePaletteLists` a lazy getter (same name and type); optional fields on the
    continuous and discrete props and specs: `scale` (`linear`, `log`, `sqrt`, with log's non-positive rule) and `mid`
    (a diverging center); on the categorical: `unlisted: 'hash'` (FNV-1a over UTF-8 into the list, or `others`) beside
    today's missing color, the default unchanged for MolViewSpec;
  - `util/color/spaces/oklab.ts` (new, beside `lab.ts`): conversions, `Color.interpolateOklab`, oklab distance, and the
    dark-variant derivation (lightness mapped into a band, hue kept, chroma clamped);
  - `util/color/cvd.ts` (new): color-vision deficiency simulation (Machado et al. 2009) and a distinctness check over
    a palette (the conformance threshold as a parameter);
  - `util/color/lists.ts` unchanged: every list the built-in schemes name is there already (`okabe-ito` included);
  - `util/color/residue-schemes.ts` (new): the residue and nucleotide dictionaries and groups (Clustal X's groups,
    Zappo, Taylor, helix, strand, turn, buried, nucleotide, purine/pyrimidine, Okabe–Ito nucleotides) as data with
    provenance, which Seq\*'s schemes expand and a Mol\* theme may read;
  - the hydrophobicity scale moves from `@molstar/model` (`ResidueHydrophobicity`, `structure/model/types.ts`) to
    `@molstar/core` beside the residue schemes, `@molstar/model` re-exporting it under its name (non-breaking), so the
    `hydrophobicity` scheme reads Mol\*'s one table;
  - `Color.luminance` and `Color.contrast` exist already: Seq\*'s `render/tokens.ts` drops its own `luminance` for them;
  - `util/legend.ts` unchanged: Seq\*'s legends convert to `TableLegend` and `ScaleLegend` (steps as a table of ranges),
    the converter in `seqx-core`; a new legend kind would widen a union Mol\*'s code switches on.
- **`seqx-viewer`**: the `legends` state (published on change only), the `Legend` part with `mountLegend`, its default
  mount (`SeqViewerSpec.ui.legends`) and its gestures as the `legend` behavior provider, the scheme picker in the track
  menu and the legend's title with previews, the commands above, `theme.cvdSafe`, the legend highlight channel and the
  hover's entry; forced colors (with accessibility).
- **`seqx-examples`**: the p53 themes view (the protein's letters under each fitting scheme; domains by a categorical
  field with its legend; a per-residue confidence as binned `plddt` on the letters and as a continuous line; two
  values tracks naming one scale; variants by consequence, sized by a log-scaled frequency; secondary structure; light
  and dark side by side) and its editing view (a scheme previewed then picked, a shared scale's domain set and undone
  in one step for two tracks), each three ways, with a headless twin.
- **Storybook**: a Design page "Themes and legends" from this design; a schemes gallery (light, dark, simulated
  color-vision deficiencies); the legend part placed three ways (a track header, a panel, the export's band).
- **`_spec` sections**: core "Language" (the theme slot, `scale`), a new "Themes" section, "Tokens", "Layer types"
  (the channels each reads), "Export" (legends), "Issues"; viewer "Legends", "Commands", the Configuration table;
  [features](features.md) rows Themes, Residue schemes, Shared domains, Legends, Size channels and the small layers;
  ux-foundations section 8's new row.
- **Edge cases with outcomes here**: kept 52, 55, 56, 57, 67, 128, 132, 137, 138, 139, 201, 202, 203, 204, 205, 215,
  232, 234, 242 (shape, size, line style, opacity); changed 133 (the order written by the builder, else the hash), 134
  (its FNV-1a over UTF-8 restored; a dictionary palette hashes unlisted values into `:others`, EC-203 deciding), 136
  (fitting schemes and dark variants kept; a pick is a patch, undoable and saved, not UI state cleared by a new
  document); fixed 130, 135.

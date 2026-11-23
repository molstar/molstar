# Extension openings: the descriptor and `ext:` kinds, custom navigators, legends' surface

A draft (XD, 2026-11-19) for the fourth slice's lanes EX, CN and LE. One design in three parts, each cut whole by Gate 0's
answer 2 (CN first, then LE); part 1 stands alone.

## Intent

The author, 2026-11-17: "what we can make configurable, we should"; 2026-11-19 on the genome ring: "Rondo-like,
extension": the fifth slice's showcase is an extension package on the public API, so what it needs from the core lands
here first ([version-1](../../../design/version-1.md): "`ext:` layer kinds … land before it"). On legends, 2026-11-19:
changing colors and ranges from a legend "is an application concern".

The design says what an outside package may ship, how a spec names it, what a reader without it draws, and which of
the library's parts (a navigator's picture, its window, the legend part) a host may replace or place. Built-ins use the
same points ([extensibility](../../../design/extensibility.md): "a point is real only if its first user is the library
itself").

**What exists** (read in the code, 2026-11-19): the registry holds `format-runtime`, `layer-type` and `layout` entries,
checked by name (same object again no change, another an error, all listed) in `composeSeqContext`; a `LayerType` has
a `fallback` and the renderer walks a layer's `fallbacks` chain (`seq.render.layer.type` when none is drawn); a layout
the composition lacks is laid out as tracks (`seq.render.layout`). But the resolve side's `LayerTypes` (the shape each
type reads, `view/view.ts`) is a static table; `layer.type`, `section.layout` and `navigator.content` are closed enums
in the language, so a spec cannot name an extension's type, and an older reader meets a newer built-in type as an enum
error, not EC-235's "unavailable". `repr` already keeps an unknown kind in a foreign namespace as a warning
(`repr.schema.kind.extension`), and `custom.<ns>` fields are carried untouched. Sizings (`heights`, `trackSizing`)
are designed in [layouts](../../../design/layouts.md), not built. The navigator's picture is a fixed switch on
`content` (`pictureOf`, `render/navigator.ts`) that builds spec tracks; `render/window.ts` is reachable through the
`./*` export but not stated public. Legends: `LegendModel`, `legendRows`, `mountLegend(element, viewer, { of,
maxItems })`, the export's `legend` band.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| The current Seq\* | namespaces, registry kinds, conflict checks, a mandatory fallback per extension kind, version upgrades | a fixed registration order with priorities; a "kernel" layer under every extension |
| Mol\* v6 `PluginSpec` | entries as data, atomic, counted, undoable; behaviors with params; one composition root | no descriptor: a plugin's parts are spread by hand (no namespace or version check) |
| JBrowse 2 plugins | one plugin object installing several pluggable element types (adapters, track types, displays, renderers, widgets); a config naming an unknown type fails that track only | registration by `install(pluginManager)` code (here data); runtime plugin loading by URL |
| igv.js `registerTrackClass`, HiGlass plugin tracks | a track type by name | a global registry; registration on import (refused, extensibility.md) |
| `repr` (this codebase) | kinds `ns:kind`; an unknown foreign kind kept, warned, skipped; `custom.<ns>` carried | — |
| UCSC, IGV, the current Seq\* navigation track | an overview strip, an ideogram, a miniature: pictures under one window | a navigator as twelve switches (navigators.md) |
| Vega-Lite legends, matplotlib | `orient`, `legend: null` per encoding, `direction`, `columns` | a legend placed over the plot (ux-foundations section 9) |
| Mol\*'s legend control (`core/util/legend`) | a host shows Seq\*'s legends in its own control (the converter, built) | — |

**Interactions in the precedents**: a user of JBrowse or HiGlass sees a missing plugin's track as an error box (here:
the fallback drawn, the reason said once); a user drags an overview's window (built, any picture); hides a legend
(Vega-Lite by spec, built here as a layer param, below); collapses a legend panel (the host's chrome, refused as a
library gesture: the part is the host's); a JBrowse user adds a plugin's track from a menu (planned: a host's header
entry, track header parts).

## Part 1: the extension descriptor and `ext:` kinds (lane EX)

### Decision

- **An extension is a descriptor, data plus providers**, given at composition beside the registry:

  ```ts
  interface SeqExtension {                 // @molstar/seqx-core/context
    readonly id: string;                   // its namespace: `rondo`
    readonly version: string;              // semver of the package
    readonly requires: string;             // the language versions its contract targets: '0.5'
    readonly doc?: string;
    readonly registry: readonly SeqRegistryEntry[];   // layer types, layouts, sizings, formats, schemes, pictures
    readonly params?: Readonly<Record<string, ParamDefs>>; // its `custom.<id>` params by the type that reads them
    readonly config?: readonly SeqConfigItem<any>[];  // its settings, keys under its namespace
  }
  interface SeqViewerExtension extends SeqExtension {   // @molstar/seqx-viewer
    readonly behaviors?: readonly SeqViewerSpec.Behavior[];
    readonly commands?: …; readonly bindings?: …; readonly headerActions?: …;  // as the built-ins' kinds
  }
  // SeqContextSpec.extensions / SeqViewerSpec.extensions: readonly SeqExtension[]
  ```

  Why a descriptor and not entries spread by hand: the namespace check, the version check and the undo need the
  extension as one thing; why a spec field and not a registry entry: registry entries are pure data with no behaviors
  ([viewer-composition](../../../design/viewer-composition.md)), and an extension carries behaviors. `createSeqContext`
  takes either type and ignores the viewer's fields.
- **What an extension may ship** (this version): layer types, section layouts, sizings, format runtimes, residue
  schemes and measures, navigator pictures (part 2), header actions, commands, bindings, behaviors, config items. Not
  node kinds or data shapes: an extension's layer reads one of the six shapes (`letters`, `ranges`, `values`, `cells`,
  `pairs`, `grid`), its data from the built-in kinds through a format runtime. Why: resolve, the store, queries,
  selection and the second consumer know those shapes; a new shape is a new node kind, opened when an extension needs
  one (Rondo's bin pairs are `pairs`).
- **Names are namespaced by the descriptor's id**: every entry an extension ships is named `<id>:<name>`
  (`rondo:ribbons`), checked at composition (`seq.context.extension.namespace`, an error listing every offender).
  Built-ins keep bare names (`ranges`, `tracks`); `seq` is reserved. `ext:` in the designs is the placeholder for an id,
  not a literal prefix. Why: two extensions cannot collide on a shared `ext:`, a built-in cannot be shadowed, and an
  issue names the package a spec needs.
- **Registration is the registry's, the descriptor whole**: every conflict across all extensions and the registry found
  first, one error, nothing changed; the same descriptor twice counted; another under one id an error; a `requires`
  that does not cover the language's version (0.x: the minor is the major) refused
  (`seq.context.extension.version`). `viewer.extend(ext)` registers one at run time and returns its undo (behaviors
  disposed, entries removed); `retire` undoes it.
- **A spec names an extension type in the built-in slot, with a fallback**:

  ```ts
  { kind: 'layer', id: 'links', params: { type: 'rondo:ribbons', fallback: 'pairs', data: '/data/hic' },
    custom: { rondo: { bundle: 0.85, stems: true } } }
  ```

  `layer.type`, `section.layout`, `section.heights` and `navigator.content` take a built-in name or `<ns>:<name>` (a
  union of the enum and the namespaced pattern); a layer's `fallback` names a built-in type or `none`; the extension's
  own params live in `custom.<id>`, validated by its `params` when it is registered (its issues in its namespace),
  carried untouched otherwise. Why the same `layer` node and not a kind `rondo:ribbons`: its `data`, theme channels,
  `hidden`, id and address stay the language's, checked by every reader, so the fallback draws the same data in the
  same colors, and patches, history, legends and addresses need nothing new.
- **What a reader without it draws** (EC-235): the layer's chain is its type, then its `fallback`, then the built-in
  chain from there (`bands` → `ranges`); the first the composition draws is drawn, its shape checked against the data
  (`seq.view.layer.shape` when it cannot read it). Said once per namespace and view, a warning:
  `seq.view.layer.extension` ("3 layers need `rondo` 1.2: drawn as `pairs`"). `fallback: none` or a chain with nothing
  drawn: not drawn, the same issue. A namespaced type with no `fallback`: `seq.lang.layer.fallback` (warning), read as
  `none`. Layouts fall back to `tracks` (fixed, as built), sizings to `content`, navigator pictures to `axis`, schemes
  to the layer type's default color, all said. **A format runtime has no fallback**: a source in a format the reader
  lacks is unavailable, its layers say so. Why: a reader can draw data more simply, never invent it; a writer who wants
  portability ships the data in a built-in format.
- **The same rule covers a newer built-in** (the symmetric case, EC-235 proper): an unknown bare `layer.type` in a
  spec of a newer minor version is unavailable with its `fallback`, a warning, never the enum's error. So the language
  check reads `type` as a string and reports an unknown value itself (`seq.view.layer.type`, warning), and writers of
  0.5.0 on give every new built-in type its `fallback` in the builder.
- **Versions**: the spec's optional `extensions` param (`{ rondo: '1.2' }`, the least version it was written for),
  written by the builder when it writes an extension's type. A reader with the namespace at a lower major (0.x:
  minor), or without it, draws the fallback and says the version; at a lower minor in the same major draws it and says
  so (info: params it does not know are kept in `custom`). No migrations before 1.0 (version-1: "migrations start at
  1.0"). Why in the spec: the one place a reader learns what a figure needs before drawing it.
- **Registered late, drawn again**: the draw lists' key holds the composition's type table version, so a layer drawn
  as its fallback is drawn as its type once `viewer.extend` lands, and as the fallback again after its undo.
- **The contract an extension imports is one entry**: `@molstar/seqx-core/extension` re-exports what an extension may
  use (`LayerType`, `LayoutContext`, `StoreReader`, `PassWriter`, `Hit`, `Emphasis`, `Legend`, the channels'
  `ThemeContext`, `Placement` and the warp for `scope: 'section'` layers, `SectionLayoutType`, `FormatRuntime`, the
  navigator picture and window types, the descriptor); everything else is internal and may change in a 0.x minor.
  Checked: the example extension imports nothing else (its package's `noRestrictedImports`). Why: freezing every type
  under `./*` would freeze the renderer; one entry is what a review reads before a break.
- **Hits are closed**: an extension's `hit` gives an element, a position or both; it adds no `HitTarget`. Why: the
  viewer's behaviors switch on the targets; an extension's gestures are its own behaviors over its own commands.
- **Enough for the genome ring?** Its package ships `rondo:ribbons` (a `scope: 'section'` layer over `pairs` on the
  `rings` layout, through the warp), a bin-pair format runtime (`rondo:binpairs`, or the core's pairs level, fifth
  slice), a theme measure per block, a story channel only if the core opens one (not here). What it cannot be: the
  square navigator over several systems (part 2) and a ring zoom channel; both are the fifth slice's core work.

### Interaction matrix

| Row | Cells |
|---|---|
| a track or layer (extension type) | every cell as the built-ins': select, extend, add/remove, focus through its `hit`'s element or position; hide (`layer.hide`, `track.hide`); edit `custom.<id>` by `layer.set` (one patch); drawn as its fallback, the fallback's gestures |
| an element of an extension layer | as any element: its key and extent from `hit`, `items` for announcements; no gesture of its own unless its behavior binds one, through the grammar ([ux-foundations](../../../design/ux-foundations.md) section 8) |
| the view | no new cell: the fallback is said in the issues panel, never on the canvas |

### Configuration

| Item | Where | Default | Why; alternatives |
|---|---|---|---|
| extensions | `spec.extensions` (context and viewer) | none | listed at the root; run time `viewer.extend` |
| a layer's fallback | spec param `fallback` | required for namespaced types (absent: `none`, warned) | the writer knows what reads alike; a per-type default from the descriptor (unknown to a reader without it) |
| the spec's needs | spec param `extensions` | written by the builder | a reader says what is missing; none (versions unknown) |
| fallback of layouts, sizings, pictures | fixed: `tracks`, `content`, `axis` | — | every section can be laid out as tracks; one less param a writer gets wrong |

### Tensions

- **Params in `custom` against first-class params**: `custom.rondo.bundle` is clumsier than `bundle` and the builder
  must type it. Recommendation: `custom`, since a reader without the extension must not see unknown params as errors,
  and the extension's builder helper hides it.
- **A fallback written per layer against one per type**: a type's natural fallback is the descriptor's knowledge, and a
  writer may get it wrong. Recommendation: the builder helper the extension ships writes it; the check warns when a
  registered type's declared `fallback` differs from the written one.
- **One entry point freezes a contract under 0.x**: `LayoutContext` has twenty-odd fields, many for built-ins (`chip`,
  `rowColumns`, `bandTop`). Recommendation: export it whole, mark built-in-only fields `@internal`, so an extension sees
  the stable subset in its types.

### Self-check

Claim: "this part covers what an outside package needs to add a way to draw, and what every other reader does with it."
1. Rests on: the extension's data has a built-in shape; the reader can draw something simpler; one composition knows
   the extension at start.
2. Negated: data with no built-in shape; a reader that can draw nothing like it; an extension that comes and goes at
   run time; a newer core meeting an older extension; a newer spec's built-in type meeting an older reader.
3. The negated world surfaced: the symmetric case of a newer built-in (EC-235 proper, the enum's error); run-time
   registration and its undo redrawing the fallback; the version both ways (`requires` refused at composition, the
   spec's `extensions` said at reading); formats, which have no fallback.
4. Integrated: all four above into the Decision. Standing: new data shapes (deferred with node kinds); an extension's
   story channels.

### Consequences

- **`seqx-core`**: the language 0.5.0's union for `type`, `layout`, `heights`, `content`, the `fallback` and
  `extensions` params, `seq.lang.layer.fallback`; `SeqExtension`, `spec.extensions`, the namespace and version checks;
  resolve reads the composition's layer kinds (shape and fallback per type) instead of the static `LayerTypes`, and
  `ResolvedLayer` carries `custom`; `seq.view.layer.extension`; the `extension` export. Tests: an extension type drawn,
  its fallback without it (canvas and SVG alike), the issue once a namespace, the shape mismatch, a newer bare type,
  registration atomic over two extensions, the version refused; property: a spec with an extension prints, patches and
  diffs identically with and without it. `_spec` "Context", "Language", "Render: Layer types".
- **`seqx-viewer`**: `SeqViewerExtension`, `viewer.extend` with its undo, the redraw on the type table's version.
- **Storybook**: `Features/Extensibility` with an example extension (a small layer type with params in `custom`, a
  scheme, a command and a binding), toggled at run time to show the fallback; its package's restricted imports.
- Sizes: the minimal viewer not grown (the extension field is a list).

## Part 2: custom navigators (lane CN; cut first)

### Decision

- **A navigator's picture is a registry entry** (`navigator-picture`), named and namespaced as part 1; the built-ins
  (`axis`, `tracks`, `bands`, `rows`) are entries written to the same contract, listed by the default spec.
  `navigator.content` names one; its params in `custom.<id>`; a reader without it draws `axis`, said.
- **A picture is spec tracks**, drawn by the navigator's own renderer through the layout interface as now:

  ```ts
  interface NavigatorPicture {
    readonly name: string;
    readonly axes: readonly ('positions' | 'rows' | 'both')[];   // what it can cover
    readonly size: number | readonly [number, number];           // its default thickness, or a miniature's box
    picture(ctx: {
      model: StaticModel; view: ResolvedView; navigator: ResolvedNavigator;
      bound: ResolvedSection; boundNode: Node; height: number; zoom: boolean; config: RenderConfig;
    }): { readonly tracks: readonly Node[]; readonly rows?: string } | undefined;
    readout?(position: number, ctx: …): string | undefined;      // the hover readout; the position's label when absent
  }
  ```

  Why tracks and not a draw function: a miniature drawn by other code drifts from the view it summarizes
  (navigators.md); as tracks the picture reads at `view` priority, colors by the view's held domains, draws any
  registered layer type (an extension's too), exports, and needs no draw code. `zoom: true` asks for the zoom strip's
  picture; `undefined` falls back to the axis there.
- **The window stays the navigator's**, whatever the picture: its geometry, floor, parts, gestures, keys, idle rule,
  zoom strip and announcements. Why: navigators.md's "one part"; a picture that moved the window would fork them.
- **The window is public at two levels**:
  - **drawn by the library**: `render/window.ts`'s `floorExtent`, `navigatorPart`, `drawNavigatorWindow`,
    `NavWindowFlags` and `NavigatorWindow` (the frame: true and drawn ends, span, axes, positions, rows), in the
    `extension` entry and stated in the `_spec`, for a picture or a layout that draws a window of its own (a ring's);
  - **drawn by the host**: `NavigatorModel` (headless) is the contract for a navigator a host draws itself (DOM,
    WebGL): `frame.windows`, `values`, `readout`, `cursor`, and the input methods (`resize`, `hit`, `move`, `down`,
    `up`, `click`, `doubleClick`, `wheel`, `key`, `cancel`, `leave`, `setFocus`, `focusNext`). Why: the model already
    holds every rule; a host replacing the canvas keeps them by calling it.
- **Refused: a navigator over several systems** (Rondo's square, a genome's chromosome grid): navigators.md refuses
  switching systems; a picture draws the navigator's one system. A navigator kind over systems is new geometry (a
  window per system pair), the fifth slice's to design.

### Interaction matrix

| Row | Cells |
|---|---|
| a navigator (custom picture) | every cell as built (navigators.md's table): click centres, drag pans, an end zooms, outside sweeps, double-click fits, wheel, Esc, keys, back; hover a readout (the picture's `readout`); hide (`hidden`); edit: no, a navigator has no content to edit |
| a host-drawn navigator | the same cells through `NavigatorModel`'s methods; the host maps its events, the rules are the model's |

### Configuration

| Item | Where | Default | Why; alternatives |
|---|---|---|---|
| the picture | spec param `content` | `axis` (a rail: none declared) | built-ins as entries; a host's own |
| a missing picture's fallback | fixed: `axis` | — | every system has an axis |
| the built-in pictures | `DefaultSeqViewerSpec().registry` | all four | a minimal spec lists fewer |

### Tensions

- **Spec tracks cannot draw what no layer type draws** (a picture of a different kind of mark): an extension then ships
  the layer type too. Recommendation: accepted; it keeps one drawing path.
- **`NavigatorModel` as a contract** freezes eighteen methods. Recommendation: public, since the gestures are already
  bound by the grammar and tested headless; any rename after is a breaking change recorded.

### Self-check

Claim: "a host can make any navigator it needs." Rests on: one system, one window per bound section, a picture
expressible as tracks. Negated: several systems; a picture no layer draws; a host that will not use the canvas. The
negation surfaced the host-drawn level (integrated: `NavigatorModel` public) and the layer type an extension ships
(integrated); several systems stays refused, a standing tension handed to the fifth slice.

### Consequences

- **`seqx-core`**: `navigator-picture` entries; `pictureOf`'s switch becomes the four built-in entries; `content`'s
  union; the window types in the `extension` entry. Tests: a custom picture drawn and its window driven headless; the
  fallback to `axis`; a picture whose tracks use an extension layer; the built-ins' tests unchanged (the refactor's
  check). `_spec` "Render: Navigators".
- **`seqx-viewer`**: `NavigatorModel` stated public in the `_spec` "Composition, Navigators"; a test drives a
  host-drawn navigator through the methods alone.
- **Storybook**: `Features/Navigators` gains a custom picture (a coverage strip) and a host-drawn DOM navigator.

## Part 3: legends' library surface (lane LE; cut second)

### Decision

- **Placement is the host's** (as built): `mountLegend(element, viewer, options)` and `legendRows(viewer, options)`
  for a host's own components, the Mol\* converter for a Mol\* control; never on the canvas (ux-foundations section 9).
  New options:

  | Option | Values | Default | Why |
  |---|---|---|---|
  | `of` | an address or a list (a track's, a layer's, a scale's) | every legend | a legend beside its track, or a panel of some |
  | `omit` | a list of addresses | none | the panel minus the ones a host shows elsewhere |
  | `form` | `full`, `compact` | `full` | below |
  | `orient` | `vertical` (groups stacked), `horizontal` (groups in a row, wrapping) | `vertical` | a side panel, or a band under the canvas |
  | `maxItems` | number | 12 | as built |
  | `channels` | `reserve`, `collapse` | `reserve` | as decided (themes.md) |

- **Hiding has two owners**: the figure's, a layer param `legend` (`auto`, the default; `none`), honoured by the part,
  `legendRows` and the export (a constant-color layer whose key says nothing); the host's, `of` and `omit`. **No
  session command hides a legend**: a user hides the panel through the host's chrome. Why: a legend is derived from the
  spec; what a figure keys is the figure's, where it sits is the host's, and a third, session owner would make
  "why is there no key" unanswerable.
- **The compact form** is one line per legend: its title, then swatches without labels (each labelled for the
  accessible name and on hover), a ramp as a short bar with its domain's ends, steps as joined cells, "+n" past
  `maxItems`. Same entries, roving focus, gestures and hover outline as the full form; only the drawing differs. Why: a
  key beside a track header or in a toolbar has a line, not a panel (track header parts' planned widget slot).
- **Editing stays the application's** (decided): commands (`theme.value-color`, `scale.set`, `layer.channel`), no
  control in the part.
- **The export** keeps `none`, `below`, `right` and the full form. Why: a figure's key is read on paper, labels needed.

### Interaction matrix

| Row | Cells |
|---|---|
| a legend entry | as themes.md's row, in both forms (compact: the label on hover and in the accessible name) |
| a legend | hide: `legend: none` on the layer (one patch, one undo step) or the host's `omit`; others no, a legend has no axis |

### Configuration

The options above; the layer param `legend` (`auto`; `none`), spec, since the figure owns its key.

### Tensions

- **A layer param for a viewer's part**: `legend: none` is read by a host's part. Recommendation: accepted, it is the
  figure's key (the export reads it too), as `export.legend` is.

### Self-check

Claim: "a host can show the legends it wants where it wants." Rests on: a host has room for a panel; every legend is
worth keying; the host decides. Negated: no room (the compact form, integrated); a legend not worth keying (the layer
param, integrated); the writer decides, not the host (the param, integrated); a user decides (refused, standing: the
host's chrome).

### Consequences

- **`seqx-core`**: the `legend` layer param (0.5.0), read by `layoutLegends` and the legends' derivation.
- **`seqx-viewer`**: `mountLegend`'s and `legendRows`' options, the compact form. Tests: `of` as a list, `omit`, the
  layer param in the part and the export, the compact form's keys and accessible names, `orient`.
- **Storybook**: the legend part placed four ways (a track header line compact, a side panel, a band under the canvas
  horizontal, the export).

## Gaps for the author

1. Namespaces are the extension's id (`rondo:ribbons`), `ext:` a placeholder, built-ins bare? Recommend **yes**.
2. An extension's params in `custom.<id>` rather than namespaced keys in `params`? Recommend **yes**.
3. The spec's `extensions` param (what a figure needs, by version)? Recommend **yes, optional, written by the builder**.
4. A namespaced layer type without `fallback`: a warning read as `none`, not an error? Recommend **yes** (the rest of
   the spec opens).
5. One `@molstar/seqx-core/extension` entry as the only contract an extension may import? Recommend **yes**.
6. `NavigatorModel` public as the host-drawn navigator's contract? Recommend **yes** (CN).
7. A navigator over several systems (Rondo's square) stays refused here, designed in the fifth slice? Recommend **yes**.
8. Legends: no session command to hide one; the layer's `legend: none` and the host's `omit` only? Recommend **yes**.
9. The CLI loading an extension by module path? Recommend **no**: a CLI host composes its own context (no code loaded
   from a spec or a flag).

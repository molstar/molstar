/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout's settings (_spec/README.md "Render", Configuration): every number is a setting with its recommended
// default; the density steps and their thresholds, shared by every layer type.

import { DefaultMaxStrings } from './draw-lists.js';
import type { SeqTokenName } from './tokens.js';

/** What a cell shows as it narrows, at fixed thresholds per position, without moving anything. */
export type DensityStep = 'letters' | 'cells' | 'dots' | 'bars' | 'density';

export const DensitySteps: readonly DensityStep[] = ['letters', 'cells', 'dots', 'bars', 'density'];

/**
 * The lowest pixels per position of each step, in the unit `DensityUnits` gives it: letters from 10 CSS px, cells
 * from 3 device px, dots from 1, bars from 1/64.
 */
export interface DensityThresholds {
  readonly letters: number;
  readonly cells: number;
  readonly dots: number;
  readonly bars: number;
}

/** The unit of each step's threshold: CSS pixels (what a reader sees: letters) or device pixels (what bounds lists). */
export type DensityUnits = { readonly [K in keyof DensityThresholds]: 'css' | 'device' };

export interface RenderConfig {
  /** The page's inset around everything, in CSS px. */
  readonly inset: number;
  /**
   * The composed extensions' settings by key (`rondo.stemPx`), as the composition's `config` gave them: what an
   * extension's layer types read through `settingOf` (_spec/README.md "Context", Extensions, Settings). Absent: none.
   */
  readonly settings?: Readonly<Record<string, unknown>>;
  /**
   * The one label column: its width, and the gap between it and the data. With `labelColumn.measure` (the default) the
   * width is the widest track label's in the view (a direction mark fades at its end), between `labelColumn.min` and
   * `labelColumn.max` (80 and 160 px), and a label past it ends in an ellipsis; without, `labelWidth` (112) is the width
   * (_spec/README.md "Render", Configuration: decided from the defaults audit).
   */
  readonly labelWidth: number;
  readonly labelGap: number;
  readonly labelColumn: { readonly measure: boolean; readonly min: number; readonly max: number };
  /** Gaps between the tracks of a section (a group) and between sections. */
  readonly trackGap: number;
  /** The gap before and after a group of tracks (_spec/README.md "Render", Track groups), in place of `trackGap`. */
  readonly groupGap: number;
  /**
   * A folded group's summary (_spec/README.md "Render", Track groups): a values member's strip `stripPx` high, at most
   * `maxStrips` of them (the rest counted in the header).
   */
  readonly group: { readonly summary: { readonly stripPx: number; readonly maxStrips: number } };
  readonly sectionGap: number;
  /**
   * A track's padding when its spec gives none (_spec/README.md "Render", Configuration): px inside its well above and
   * below its data (0, 0). The well is the track's height plus its padding; its data's band is inset by it.
   */
  readonly trackPadding: { readonly top: number; readonly bottom: number };
  /**
   * A track's well (_spec/README.md "Render", Shapes): a squircle by the squircle rule over the well's height, its data
   * clipped to it (`squircle`, the default), or a square box (`square`).
   */
  readonly trackShape: 'squircle' | 'square';
  /** A track's height when the spec gives none: the largest of its layers' types'. */
  readonly heights: {
    readonly letters: number;
    readonly ranges: number;
    readonly values: number;
    readonly ruler: number;
    readonly secondary: number;
    readonly variants: number;
    /** A frame's title row (its outline covers the tracks it names). */
    readonly frame: number;
    /** An alignment's rows band when it gives no height (its rows scroll inside), and a logo. */
    readonly rows: number;
    readonly logo: number;
    /** A pairs layer's arcs band (in a grid its dots take the plot). */
    readonly pairs: number;
    /** A grid layer's track when a grid section's `aspect` is `none` (else the plot's height is the aspect's). */
    readonly grid: number;
  };
  /**
   * A lane of ranges (14 px: a 12 px mark and its 2 px gap) and the padding above and below a band's lanes (4 px): a
   * track with no height is as tall as its lanes, up to `track.maxPx` (_spec/README.md "Render", Track heights).
   */
  readonly laneHeight: number;
  readonly laneGap: number;
  readonly lanePadding: number;
  /**
   * A track of lanes past its content height (_spec/README.md "Render", Track sizing): `maxPx`, the cap (240 px: 16
   * lanes and padding, a px number: no view height to take a share of), and `overflow` past it: its lanes `scroll`
   * inside the track (the default), `compact` (squeezed to the band) or the rest in a `density` row.
   */
  readonly track: {
    readonly maxPx: number;
    readonly overflow: 'scroll' | 'compact' | 'density';
    /**
     * A track with no items in what it covers (its content known whole, no lanes): `note` (the default) a one-line row
     * of `emptyPx` naming it with a quiet "no items", `collapse` the row alone, `keep` its well as any track's. A
     * track's own `empty` overrides it; a track still waiting for its data is pending, never empty.
     */
    readonly empty: 'note' | 'collapse' | 'keep';
    /** The height of an empty track's row (px). */
    readonly emptyPx: number;
  };
  /**
   * An export's lanes (`all`: every lane of a capped track at its content height; `view`: the band as on screen with
   * its cut edges), and the most a track's px an export draws (4,000), the rest said in its header.
   */
  readonly export: { readonly lanes: 'all' | 'view'; readonly maxTrackPx: number };
  readonly density: DensityThresholds;
  readonly densityUnits: DensityUnits;
  /** Letter sizes, a short scale (CSS px); the largest that fits the cell is used. */
  readonly letterSizes: readonly number[];
  /** The mono face's advance per em, so a letter's size is known from its cell with no measuring. */
  readonly monoAdvance: number;
  /** The label font (a family role, `sans` or `mono`: a backend maps it to its stack) and size. */
  readonly labelFont: string;
  /** Items' labels inside their marks. */
  readonly itemFont: string;
  /**
   * Items' labels inside their marks: one that fits neither whole inside nor beside goes inside cut with an ellipsis
   * while at least `minChars` letters are kept (4; a label of so many or fewer only whole), else is dropped.
   */
  readonly labels: { readonly minChars: number };
  /**
   * The letters layer's mark of a modified position (`modified: mark`, _spec/README.md "Render", Layer types): a dot
   * `size` px at the cell's top right (`shape: 'dot'`) or a bar along its bottom (`underline`), in the letter's ink; a
   * bar over a component of several positions, a tick at a strip's top below the cells.
   */
  readonly letters: {
    readonly modifiedMark: { readonly shape: 'dot' | 'underline'; readonly size: number };
    /**
     * The largest letter size, CSS px of the mono face: a letter is never larger than the chrome's label (cap heights
     * equal: 0, the default, is the label font's size times `capRatio`, so it follows `labelFont`). A host that wants
     * letters to follow the cell (teaching, slides) sets it past the scale's largest size.
     */
    readonly maxPx: number;
    /** The label face's cap height over the mono face's, per em: measured by the host (the canvas); 1 when it has none. */
    readonly capRatio: number;
  };
  /** How far past its max width a text fades out, CSS px. */
  readonly fade: number;
  /** Squircle corners: radius a third of the height, capped; none below `squircleMin` in either dimension. */
  readonly squircleRatio: number;
  readonly squircleMax: number;
  readonly squircleMin: number;
  /** Outline widths: hover, selection, cursor (thicker). */
  readonly outline: { readonly hover: number; readonly selection: number; readonly cursor: number };
  /**
   * The emphasis's look: where the selection's tint is drawn (`over` the data, translucent, the default; `under` it, in
   * the background pass: hidden by opaque cells, and a selection change redraws the base), its opacity, a highlight
   * channel's tint opacity, and how strong what is projected onto a linked section is against its own look (a factor
   * of the opacity: 0.5).
   */
  readonly emphasis: {
    readonly selectionPlace: 'over' | 'under';
    /**
     * The selection's, a highlight's and the projected hover's shape on a section: one band a range spanning its tracks
     * and the gaps between them (`column`, the default), or a box a track (`tracks`).
     */
    readonly selectionShape: 'column' | 'tracks';
    readonly selectionTint: number;
    /** The selection's outline: the accent at this opacity (0.6), `outline.selection` px wide (1). */
    readonly selectionOutline: number;
    /** A hovered member's row band where a layer draws one (the matrix's): its tint's opacity (0.1). */
    readonly hoverTint: number;
    readonly highlightTint: number;
    /** Stripes a band is split into where channels overlap (4); the rest share the last one. */
    readonly highlightStripes: number;
    /** The height of a `band`-style channel's bar, px (3). */
    readonly highlightBand: number;
    readonly projected: number;
    /** How far a hover's outline stands outside its mark (1 px): concentric, never touching what is under it. */
    readonly hoverOffset: number;
    /** The keyboard cursor's ink, a token (`msui-text`: neutral, never the selection's accent). */
    readonly cursorInk: SeqTokenName;
  };
  /**
   * The axis's marks: a reversed axis says so in its section's first label (`directionMark`: `label`, an arrow and
   * "reversed", or "along" what it reads along; `none` for a host whose own chrome says it). Fixed: a break is a
   * hairline in the break ink (`seq-break`) across the section's band with a notch above it, zero units wide, hit
   * within 3 CSS px of its notch and on its hidden items' ticks.
   */
  readonly axis: { readonly directionMark: 'label' | 'none' };
  /**
   * A header's line in the label column (_spec/README.md "Render", Track headers): `slots` action columns of `slotPx`
   * and the "⋯" at the column's end (2, 20 px: an icon of 14), none on a track under `minHeightPx` (16: its states
   * marks); a section's line `sectionPx` high (20); the actions drawn on `hover` (the track hovered, the cursor on its
   * header), `always` or `never` (the menu only); the states in an export (`include`) or not (`none`).
   */
  readonly header: {
    readonly slots: number;
    readonly slotPx: number;
    readonly minHeightPx: number;
    readonly sectionPx: number;
    /** A group's header line (open), and its tracks' labels' indent (_spec/README.md "Render", Track groups). */
    readonly groupPx: number;
    readonly groupIndentPx: number;
    readonly actions: { readonly show: 'hover' | 'always' | 'never' };
    readonly states: { readonly export: 'include' | 'none' };
  };
  /** A section's header line: drawn when it has a state or a slot (`auto`), for every section, or never. */
  readonly section: {
    readonly header: 'auto' | 'always' | 'none';
    /**
     * A section's box with no `height` of its own (_spec/README.md "Render", Section scroll): its content's height up to
     * `maxPx` (600 px), past it its content (a text's lines, a tracks section's tracks below its ruler) scrolls inside
     * with the one rail; an export shows all of it.
     */
    readonly maxPx: number;
  };
  /** What a want adds on each side of the view, as a fraction of it (the store's `margin` priority). */
  readonly margin: number;
  /** A follower's leader pieces carried as their image's hull, or the shortest arc on a cyclic follower. */
  readonly followArc: 'shortest' | 'hull';
  /** Tracks one template expands into at most; more is `seq.render.template.limit` and the rest are not laid out. */
  readonly maxTemplateTracks: number;
  /**
   * Positions a numbering whose labels are a source's may cover at most (2^18, 262,144: a protein's or a small
   * genome's): past it the labels are neither asked for nor held (`seq.render.numbering.limit`). Why: a source's
   * labels are held whole, a string a position.
   */
  readonly maxNumberingPositions: number;
  /**
   * Dynamic nodes: with `holdShown`, a frame's layout lists the nodes it shows (a placed instance, the node of an item
   * system a shown layer goes through, and what those name), for its owner to hold in the store (`ViewLayout.shown`).
   */
  readonly nodes: { readonly holdShown: boolean };
  /**
   * Ranges: items past the lanes collapse into a density row at the bottom (`extra-lanes`), or the whole layer does;
   * the connector joining an item's pieces across shown positions (_spec/README.md "Render", Layer types): `width` CSS
   * px (1), `ink` the item's own colour (`item`) or the muted ink (`muted`), at `alpha` of it (0.6).
   */
  readonly ranges: {
    readonly overflow: 'extra-lanes' | 'layer';
    /** The most lanes a layer packs (1,000): past them the last lane is the density row (layout work is bounded). */
    readonly maxLanes: number;
    readonly connector: { readonly width: number; readonly ink: 'item' | 'muted'; readonly alpha: number };
  };
  /**
   * Values: below a device pixel a position, each column's extent (`minmax`, negative extremes show) or its maximum;
   * how many times a device pixel's positions a level's bin may be (`coarser`, 1 by default: the coarsest bin no wider than a
   * pixel's positions, exact to a device pixel; 2, up to two pixels'); whether data on the section's own system is placed directly (`direct`;
   * false: through the general mapping, the same lists, slower).
   */
  readonly values: { readonly aggregate: 'minmax' | 'max'; readonly coarser: number; readonly direct: boolean };
  /**
   * The ruler: the least px between two ticks (`labelGap`, 48: the foundations' number; the step of 1, 2, 5 × 10ⁿ is the
   * least that reaches it). Labels also clear each other by an eighth of it; one that does not fit is skipped, never shrunk.
   */
  readonly ruler: { readonly labelGap: number };
  /**
   * Blocks in a row and on a ring (_spec/README.md "Render", Blocks): the gaps in view share at most `gapShare` of the
   * length (0.25, each at most the layout's gap: 12 px, 1°); consecutive blocks each narrower than `runPx` (2 CSS px) are
   * one run (a band, no gaps inside, labelled when it has room; 0: never); inside a run a block is read and drawn while
   * at least `readPx` device px wide (1; 0 reads every block); a gap between two assemblies' groups is `groupGap` block
   * gaps (3, taken first, fixed).
   */
  readonly blocks: {
    readonly gapShare: number;
    readonly runPx: number;
    readonly readPx: number;
    readonly groupGap: number;
    /**
     * A stacked block's track with no items for that block: `collapse` (the default) leaves it out of that block's row
     * (the first block's row and a track of letters or any fixed height are kept), `keep` draws its well.
     */
    readonly emptyTracks: 'collapse' | 'keep';
  };
  /**
   * Rings (_spec/README.md "Render", Rings layout): the hole a share of the radius at least (0.35: room for the title,
   * chords readable), the diameter's cap (720 px: a map, not a wall), the gap between blocks (1°: a boundary, Circos),
   * letters `upright` or along the ring (`tangent`), the hole's title (`title`: the name and length; `none`), the
   * outside labels' gap from the rings and their leaders' bend (8 and 12 px: leaders fan out over the ruler), the
   * grouping of a dense fan (more than 6 labels within 12° are one label, "13 items"), and when
   * the window is drawn (`followed`: only while a section follows the ring; `always`, idle when none does).
   */
  readonly rings: {
    readonly hole: number;
    /** A lanes track's radial width past which it compacts (48 px): an annulus cannot scroll. */
    readonly trackMaxPx: number;
    readonly maxDiameter: number;
    readonly blockGap: number;
    /**
     * The share of the circle the gaps between blocks may take together (0.25: positions keep three quarters at any
     * count, each gap at most `blockGap`), and the gap between two assemblies' groups (4°, taken first, fixed).
     */
    readonly gapShare: number;
    readonly groupGap: number;
    readonly letters: 'upright' | 'tangent';
    readonly center: 'title' | 'none';
    readonly labelGap: number;
    readonly leaderBend: number;
    /**
     * More than `over` outside labels within `degrees` are one grouped label ("13 items"); 0 degrees: never grouped. A
     * hover opens a group of up to `fanMax` as a fan, pivoted wider than its own arc; a larger one is listed in its card.
     */
    readonly cluster: { readonly over: number; readonly degrees: number; readonly fanMax: number };
    readonly window: 'followed' | 'always';
    /** The runs an image's rows draw on a ring a frame (65,536): past it rows are taken a step apart (`seq.render.ring.image`). */
    readonly imageRuns: number;
  };
  /**
   * The strings (texts and fonts) a renderer's draw lists keep interned across frames: a frame drawn again interns
   * nothing new; past this many, a full redraw starts the table again.
   */
  readonly internedStrings: number;
  /**
   * The shared image cache's bytes (64 MB): every `image` tile a layer draws (grid, rows, navigators) in one LRU, this
   * over a tile's bytes its capacity (_spec/README.md "Render", The shared reduction). Why one: two caches could not lend
   * each other room.
   */
  readonly imageBytes: number;
  /**
   * The rows layer (_spec/README.md "Render", The rows layer): the least row height a label, a number and an icon are
   * drawn at (`labelMinPx`, 8: a label needs a row), and a pinned row kept in its place in the order too (`pinInPlace`,
   * true: the tree stays a tree). A row's height (`RowHeight`, 20) and a tree hit's reach (`TreeHitPx`, 8) are fixed.
   * Row columns (Row columns): a column's least width before it is hidden (`columns.floorPx`, 32: a short number or a
   * chip), the gap before each (`columns.gapPx`, 6), and the data's least width they yield to (`columns.dataMinPx`,
   * 200: the data keeps a third of a 600 px view beside a tree, a label and two columns), and the title row a rows
   * layer with columns takes above its rows (`columns.titlePx`, 16: a line of the items' font), so its titles never
   * reach the track above. While the view moves a frame builds at most `buildBudget` cells of the shared reduction's
   * new level (2^18: about 4 ms), the level drawn last drawn meanwhile (a settled frame builds the rest), as the
   * grid's `buildBudget`.
   */
  readonly rows: {
    readonly labelMinPx: number;
    readonly pinInPlace: boolean;
    /** A rows band with no `height` is as tall as its rows up to this many px (480: an alignment wants more than features). */
    readonly maxPx: number;
    readonly buildBudget: number;
    readonly columns: {
      readonly floorPx: number;
      readonly gapPx: number;
      readonly dataMinPx: number;
      readonly titlePx: number;
    };
  };
  /**
   * Navigators and the window they share with a ring's (_spec/README.md "Render", Navigators): the window's least drawn
   * width (`minPx`, at most half the span), its least hit width by pointer and by touch, its grips' width (tabs just
   * outside each end, at most a quarter of the window inside it); a floored window drawn `centred` on the truth or as a
   * scrollbar's `thumb`; a miniature's tracks at `trackScale` of their height, at least `trackMinPx`; the marks drawn over
   * the picture; the ruler rail grown to `grownPx` while used, its grab zone on an alignment's ruler its lower
   * `grabShare`; a narrow frame's order of giving way (A narrow frame): a vertical overview narrowed toward
   * `overview.minPx` (48) to keep the plot `overview.minPlotPx` (160) wide, hidden below it, then the label column and
   * the tree, then the frame scrolls.
   */
  readonly navigator: {
    readonly window: {
      readonly minPx: number;
      readonly hitPx: number;
      readonly touchHitPx: number;
      readonly handlePx: number;
      readonly floor: 'centred' | 'thumb';
      /** One stroke weight for every line of a navigator's part, CSS px (1): lit is said by the fill and the grips. */
      readonly strokePx: number;
      /** The quiet alpha of every "where" mark and an idle frame (0.4): the control stays the loudest. */
      readonly idleAlpha: number;
    };
    /** The zoom strip's funnel (two lines from the context's mark to the control row): `muted` or `none`. */
    readonly zoomStrip: { readonly funnel: 'muted' | 'none' };
    readonly trackScale: number;
    readonly trackMinPx: number;
    readonly marks: readonly ('selection' | 'highlights' | 'cursor')[];
    readonly rail: {
      readonly grownPx: number;
      readonly grabShare: number;
      /** The rail's tint: `neutral` (the text ink at 55%, the default) or `accent` (the accent at 55%). */
      readonly tint: RailTint;
    };
    readonly overview: { readonly minPlotPx: number; readonly minPx: number };
    /**
     * The square (two axes over blocks × blocks): past `maxCells` block (or run) pairs (4,096) its lattice's lines are
     * not drawn, said once; its picture is always drawn.
     */
    readonly square: { readonly maxCells: number };
  };
  /**
   * The pairs layer (_spec/README.md "Render", The pairs layer): the most arcs drawn (`maxArcs`, 4,096: bounded by the
   * screen; past it one arc a bin of the ends' device columns, widest first, the rest a density strip of ends) and how
   * near a point is to a curve to hit it (`hitPx`, 4 CSS px: a curve is thin). Read by region: bins at the coarsest
   * level no wider than `binPx` device pixels' positions (2), items asked once the bin pairs in view sum to at most
   * `maxArcs` and kept while the items in view stay at most `itemsHysteresis` × `maxArcs` (1.25: no flip at the
   * threshold).
   */
  readonly pairs: {
    readonly maxArcs: number;
    readonly hitPx: number;
    readonly binPx: number;
    readonly itemsHysteresis: number;
    /** A bin pair's stroke, CSS px, on a square-root scale of its count over the counts in view: 1 to 6. */
    readonly binStroke: { readonly min: number; readonly max: number };
  };
  /**
   * A bridge (_spec/README.md "Render", Bridges): its band's height when its spec gives none (`height`, 96 CSS px: room
   * for a crossed ribbon to read; the band's resize grip is the viewer's, not built).
   */
  readonly bridges: { readonly height: number };
  /**
   * A bridge's ribbons (_spec/README.md "Render", Bridges): filled at `fillOpacity` (0.35: overlaps read through), a bin
   * by its aligned length at `binOpacity` (0.15 to 0.8, square root over the bins in view), each end at least `minPx`
   * device px wide (1: a short run on a long axis is still seen), a stub `stubShare` of the band's height (0.25).
   */
  readonly ribbons: {
    readonly fillOpacity: number;
    readonly binOpacity: { readonly min: number; readonly max: number };
    readonly minPx: number;
    readonly stubShare: number;
  };
  /**
   * The grid layer (_spec/README.md "Render", The grid layer): a cell drawn as a rect with its outline from `rectMinPx`
   * (12 CSS px: rects only where few), image tiles of the shared reduction below; at most `buildBudget` samples (65,536:
   * a tile) built a frame while the view moves, the coarser level drawn meanwhile (a settled frame builds the rest).
   * Read by region: the finest level whose bin is at least `binPx` device pixels (1) on both axes, a ring of
   * `marginTiles` tiles (1) round the view at margin priority.
   */
  readonly grid: {
    readonly rectMinPx: number;
    readonly buildBudget: number;
    readonly binPx: number;
    readonly marginTiles: number;
  };
  /**
   * The text layout (_spec/README.md "Render", Text layout): its letter size (`scale`, px a position: the zoom), the gap
   * between groups (a share of a cell), between lines and between a line's bands (px), whether a band no layer used on
   * a line collapses there, the most lines a section lays out, whether a piece's edge starts a line, the gap between a line and its numbers, a track's height by its
   * layers' type when the spec gives none (thin: a ranges track reads as underlines, a values track as a strip).
   */
  readonly text: {
    readonly scale: number;
    readonly groupGap: number;
    readonly lineGap: number;
    readonly bandGap: number;
    readonly collapseEmpty: boolean;
    readonly maxLines: number;
    readonly pieceBreak: 'none' | 'line';
    readonly numberGap: number;
    readonly heights: RenderConfig['heights'];
  };
}

/** The rail's tint: a neutral grey or a faint accent. */
export type RailTint = 'neutral' | 'accent';

export const DefaultRenderConfig: RenderConfig = {
  inset: 8,
  labelWidth: 112,
  labelGap: 8,
  labelColumn: { measure: true, min: 80, max: 160 },
  trackGap: 4,
  groupGap: 8,
  group: { summary: { stripPx: 4, maxStrips: 4 } },
  sectionGap: 20,
  trackPadding: { top: 0, bottom: 0 },
  trackShape: 'squircle',
  heights: {
    letters: 20,
    ranges: 36,
    values: 32,
    ruler: 22,
    secondary: 16,
    variants: 32,
    frame: 16,
    rows: 240,
    logo: 48,
    pairs: 64,
    grid: 240,
  },
  laneHeight: 14,
  laneGap: 2,
  lanePadding: 4,
  track: { maxPx: 240, overflow: 'scroll', empty: 'note', emptyPx: 18 },
  export: { lanes: 'all', maxTrackPx: 4000 },
  density: { letters: 10, cells: 3, dots: 1, bars: 1 / 64 },
  densityUnits: { letters: 'css', cells: 'device', dots: 'device', bars: 'device' },
  letterSizes: [9, 10, 11, 12, 13, 14],
  monoAdvance: 0.6,
  labelFont: '12px sans',
  itemFont: '10px sans',
  labels: { minChars: 4 },
  letters: { modifiedMark: { shape: 'dot', size: 3 }, maxPx: 0, capRatio: 1 },
  fade: 12,
  squircleRatio: 1 / 3,
  squircleMax: 8,
  squircleMin: 6,
  outline: { hover: 1.5, selection: 1, cursor: 1 },
  emphasis: {
    selectionPlace: 'over',
    selectionShape: 'column',
    selectionTint: 0.18,
    selectionOutline: 0.6,
    hoverTint: 0.1,
    highlightTint: 0.3,
    highlightStripes: 4,
    highlightBand: 3,
    projected: 0.5,
    hoverOffset: 1,
    cursorInk: 'msui-text',
  },
  axis: { directionMark: 'label' },
  header: {
    slots: 2,
    slotPx: 20,
    minHeightPx: 16,
    sectionPx: 20,
    groupPx: 20,
    groupIndentPx: 12,
    actions: { show: 'hover' },
    states: { export: 'include' },
  },
  section: { header: 'auto', maxPx: 600 },
  margin: 0.5,
  followArc: 'shortest',
  maxTemplateTracks: 200,
  maxNumberingPositions: 1 << 18,
  nodes: { holdShown: true },
  ranges: { overflow: 'extra-lanes', maxLanes: 1000, connector: { width: 1, ink: 'item', alpha: 0.6 } },
  values: { aggregate: 'minmax', coarser: 1, direct: true },
  ruler: { labelGap: 48 },
  blocks: { gapShare: 0.25, runPx: 2, readPx: 1, groupGap: 3, emptyTracks: 'collapse' },
  rings: {
    hole: 0.35,
    trackMaxPx: 48,
    maxDiameter: 720,
    blockGap: 1,
    gapShare: 0.25,
    groupGap: 4,
    letters: 'upright',
    center: 'title',
    labelGap: 8,
    leaderBend: 12,
    cluster: { over: 6, degrees: 12, fanMax: 12 },
    imageRuns: 65_536,
    window: 'followed',
  },
  internedStrings: DefaultMaxStrings,
  imageBytes: 64 * 1024 * 1024,
  rows: {
    labelMinPx: 8,
    pinInPlace: true,
    maxPx: 480,
    buildBudget: 1 << 18,
    columns: { floorPx: 32, gapPx: 6, dataMinPx: 200, titlePx: 16 },
  },
  navigator: {
    window: { minPx: 16, hitPx: 12, touchHitPx: 24, handlePx: 5, floor: 'centred', strokePx: 1, idleAlpha: 0.4 },
    zoomStrip: { funnel: 'muted' },
    trackScale: 0.25,
    trackMinPx: 2,
    marks: ['selection', 'highlights', 'cursor'],
    rail: { grownPx: 6, grabShare: 0.4, tint: 'neutral' },
    overview: { minPlotPx: 160, minPx: 48 },
    square: { maxCells: 4096 },
  },
  pairs: { maxArcs: 4096, hitPx: 4, binPx: 2, itemsHysteresis: 1.25, binStroke: { min: 1, max: 6 } },
  bridges: { height: 96 },
  ribbons: { fillOpacity: 0.35, binOpacity: { min: 0.15, max: 0.8 }, minPx: 1, stubShare: 0.25 },
  grid: { rectMinPx: 12, buildBudget: 65_536, binPx: 1, marginTiles: 1 },
  text: {
    scale: 12,
    groupGap: 0.5,
    lineGap: 10,
    bandGap: 3,
    collapseEmpty: true,
    maxLines: 100_000,
    pieceBreak: 'none',
    numberGap: 8,
    heights: {
      letters: 18,
      ranges: 4,
      values: 3,
      ruler: 14,
      secondary: 6,
      variants: 10,
      frame: 14,
      rows: 18,
      logo: 24,
      pairs: 24,
      grid: 120,
    },
  },
};

/** Fields of a config to override: the nested ones whole or field by field. */
export type RenderConfigOverrides = {
  readonly [K in keyof RenderConfig]?: RenderConfig[K] extends readonly unknown[]
    ? RenderConfig[K]
    : RenderConfig[K] extends object
      ? Partial<RenderConfig[K]>
      : RenderConfig[K];
};

/** A config with some fields overridden (the nested ones whole or field by field). */
export function renderConfig(over: RenderConfigOverrides = {}): RenderConfig {
  return {
    ...DefaultRenderConfig,
    ...(over as Partial<RenderConfig>),
    nodes: { ...DefaultRenderConfig.nodes, ...over.nodes },
    axis: { ...DefaultRenderConfig.axis, ...over.axis },
    header: {
      ...DefaultRenderConfig.header,
      ...over.header,
      actions: { ...DefaultRenderConfig.header.actions, ...over.header?.actions },
      states: { ...DefaultRenderConfig.header.states, ...over.header?.states },
    },
    section: { ...DefaultRenderConfig.section, ...over.section },
    // a width given and no column rule: that width, fixed
    labelColumn: {
      ...DefaultRenderConfig.labelColumn,
      ...(over.labelWidth !== undefined && over.labelColumn === undefined ? { measure: false } : {}),
      ...over.labelColumn,
    },
    trackPadding: { ...DefaultRenderConfig.trackPadding, ...over.trackPadding },
    heights: { ...DefaultRenderConfig.heights, ...over.heights },
    density: { ...DefaultRenderConfig.density, ...over.density },
    densityUnits: { ...DefaultRenderConfig.densityUnits, ...over.densityUnits },
    ranges: {
      ...DefaultRenderConfig.ranges,
      ...over.ranges,
      connector: { ...DefaultRenderConfig.ranges.connector, ...over.ranges?.connector },
    },
    letters: {
      modifiedMark: { ...DefaultRenderConfig.letters.modifiedMark, ...over.letters?.modifiedMark },
      maxPx: over.letters?.maxPx ?? DefaultRenderConfig.letters.maxPx,
      capRatio: over.letters?.capRatio ?? DefaultRenderConfig.letters.capRatio,
    },
    values: { ...DefaultRenderConfig.values, ...over.values },
    ruler: { ...DefaultRenderConfig.ruler, ...over.ruler },
    blocks: { ...DefaultRenderConfig.blocks, ...over.blocks },
    rings: { ...DefaultRenderConfig.rings, ...over.rings },
    track: { ...DefaultRenderConfig.track, ...over.track },
    export: { ...DefaultRenderConfig.export, ...over.export },
    outline: { ...DefaultRenderConfig.outline, ...over.outline },
    emphasis: { ...DefaultRenderConfig.emphasis, ...over.emphasis },
    text: { ...DefaultRenderConfig.text, ...over.text },
    rows: {
      ...DefaultRenderConfig.rows,
      ...over.rows,
      columns: { ...DefaultRenderConfig.rows.columns, ...over.rows?.columns },
    },
    navigator: {
      ...DefaultRenderConfig.navigator,
      ...over.navigator,
      window: { ...DefaultRenderConfig.navigator.window, ...over.navigator?.window },
      zoomStrip: { ...DefaultRenderConfig.navigator.zoomStrip, ...over.navigator?.zoomStrip },
      rail: { ...DefaultRenderConfig.navigator.rail, ...over.navigator?.rail },
      overview: { ...DefaultRenderConfig.navigator.overview, ...over.navigator?.overview },
      square: { ...DefaultRenderConfig.navigator.square, ...over.navigator?.square },
    },
    pairs: { ...DefaultRenderConfig.pairs, ...over.pairs },
    grid: { ...DefaultRenderConfig.grid, ...over.grid },
  };
}

/**
 * The step a cell of `cssPx` CSS pixels per position shows at a device pixel ratio, each threshold in its unit. The
 * letters step never starts below the width the smallest letter size needs, so it never shows cells without letters.
 */
export function densityStep(cssPx: number, dpr: number, config: RenderConfig): DensityStep {
  const t = config.density;
  const u = config.densityUnits;
  const at = (k: keyof DensityThresholds) => (u[k] === 'css' ? cssPx : cssPx * dpr);
  const minLetters = (config.letterSizes[0] ?? 0) * config.monoAdvance + 2;
  const letters = u.letters === 'css' ? Math.max(t.letters, minLetters) : Math.max(t.letters, minLetters * dpr);
  if (at('letters') >= letters) return 'letters';
  if (at('cells') >= t.cells) return 'cells';
  if (at('dots') >= t.dots) return 'dots';
  if (at('bars') >= t.bars) return 'bars';
  return 'density';
}

const Caps = new WeakMap<RenderConfig, number>();

/**
 * The largest letter size: `letters.maxPx`, else the label font's size times `letters.capRatio` (the label face's cap
 * height over the mono face's, per em: 1 where the host measured none), so a letter's cap height equals a label's.
 * Computed once per config.
 */
export function letterCap(config: RenderConfig): number {
  let cap = Caps.get(config);
  if (cap === undefined) {
    cap = config.letters.maxPx > 0 ? config.letters.maxPx : fontSize(config.labelFont) * config.letters.capRatio;
    Caps.set(config, cap);
  }
  return cap;
}

/**
 * One size for a row of letters: the cap itself where the cell (px) and `room` (px, the height it may take) allow it,
 * else the largest of the scale that fits, never past the cap. The smallest size when none does (the step starts
 * where it fits a cell, so a short row still shows its letters). `capped: false` is a text figure's, read as text: the
 * scale alone.
 */
export function letterSize(config: RenderConfig, cell: number, room: number, capped = true): number {
  const cap = capped ? letterCap(config) : Infinity;
  const limit = Math.min((cell - 2) / config.monoAdvance, room);
  if (cap <= limit) return cap;
  let size = 0;
  for (const s of config.letterSizes) if (s <= cap && s <= limit) size = s;
  return size || Math.min(config.letterSizes[0] ?? 9, cap);
}

/** The size of a font's shorthand in px (12 when it names none). */
export function fontSize(font: string): number {
  const m = /(\d+(?:\.\d+)?)px/.exec(font);
  return m ? Number(m[1]) : 12;
}

/** Whether a step draws a mark per position (per cell), not per device column. */
export function perPosition(step: DensityStep): boolean {
  return step === 'letters' || step === 'cells' || step === 'dots';
}

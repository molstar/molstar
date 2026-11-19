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
  /** The one label column: its width, and the gap between it and the data. */
  readonly labelWidth: number;
  readonly labelGap: number;
  /** Gaps between the tracks of a section (a group) and between sections. */
  readonly trackGap: number;
  readonly sectionGap: number;
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
  };
  /** A lane of ranges, and the gap within it. */
  readonly laneHeight: number;
  readonly laneGap: number;
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
    readonly highlightTint: number;
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
  /** What a want adds on each side of the view, as a fraction of it (the store's `margin` priority). */
  readonly margin: number;
  /** A follower's leader pieces carried as their image's hull, or the shortest arc on a cyclic follower. */
  readonly followArc: 'shortest' | 'hull';
  /** Tracks one template expands into at most; more is `seq.render.template.limit` and the rest are not laid out. */
  readonly maxTemplateTracks: number;
  /**
   * Dynamic nodes: with `holdShown`, a frame's layout lists the nodes it shows (a placed instance, the node of an item
   * system a shown layer goes through, and what those name), for its owner to hold in the store (`ViewLayout.shown`).
   */
  readonly nodes: { readonly holdShown: boolean };
  /** Ranges: items past the lanes collapse into a density row at the bottom (`extra-lanes`), or the whole layer does. */
  readonly ranges: { readonly overflow: 'extra-lanes' | 'layer' };
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
   * Rings (_spec/README.md "Render", Rings layout): the hole a share of the radius at least (0.35: room for the title,
   * chords readable), the diameter's cap (720 px: a map, not a wall), the gap between blocks (1°: a boundary, Circos),
   * letters `upright` or along the ring (`tangent`), the hole's title (`title`: the name and length; `none`), the
   * outside labels' gap from the rings and their leaders' bend (8 and 12 px: leaders fan out over the ruler), and when
   * the window is drawn (`followed`: only while a section follows the ring; `always`, idle when none does).
   */
  readonly rings: {
    readonly hole: number;
    readonly maxDiameter: number;
    readonly blockGap: number;
    readonly letters: 'upright' | 'tangent';
    readonly center: 'title' | 'none';
    readonly labelGap: number;
    readonly leaderBend: number;
    readonly window: 'followed' | 'always';
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
   * true: the tree stays a tree). A row's height (`RowHeight`, 20) and a tree hit's reach (`TreeHitPx`, 6) are fixed.
   */
  readonly rows: {
    readonly labelMinPx: number;
    readonly pinInPlace: boolean;
  };
  /**
   * Navigators and the window they share with a ring's (_spec/README.md "Render", Navigators): the window's least drawn
   * width (`minPx`, at most half the span), its least hit width by pointer and by touch, its grips' width (tabs just
   * outside each end, at most a quarter of the window inside it); a floored window drawn `centred` on the truth or as a
   * scrollbar's `thumb`; a miniature's tracks at `trackScale` of their height, at least `trackMinPx`; the marks drawn over
   * the picture; the ruler rail grown to `grownPx` while used, its grab zone on an alignment's ruler its lower
   * `grabShare`.
   */
  readonly navigator: {
    readonly window: {
      readonly minPx: number;
      readonly hitPx: number;
      readonly touchHitPx: number;
      readonly handlePx: number;
      readonly floor: 'centred' | 'thumb';
    };
    readonly trackScale: number;
    readonly trackMinPx: number;
    readonly marks: readonly ('selection' | 'highlights' | 'cursor')[];
    readonly rail: { readonly grownPx: number; readonly grabShare: number };
  };
  /**
   * The text layout (_spec/README.md "Render", Text layout): its letter size (`scale`, px a position: the zoom), the gap
   * between groups (a share of a cell), between lines and between a line's bands (px), whether a band no layer used on
   * a line collapses there, the most lines a section lays out, the most px a section in a page takes before it scrolls
   * inside, whether a piece's edge starts a line, the gap between a line and its numbers, a track's height by its
   * layers' type when the spec gives none (thin: a ranges track reads as underlines, a values track as a strip).
   */
  readonly text: {
    readonly scale: number;
    readonly groupGap: number;
    readonly lineGap: number;
    readonly bandGap: number;
    readonly collapseEmpty: boolean;
    readonly maxLines: number;
    readonly maxHeight: number;
    readonly pieceBreak: 'none' | 'line';
    readonly numberGap: number;
    readonly heights: RenderConfig['heights'];
  };
}

export const DefaultRenderConfig: RenderConfig = {
  inset: 8,
  labelWidth: 112,
  labelGap: 8,
  trackGap: 4,
  sectionGap: 8,
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
  },
  laneHeight: 12,
  laneGap: 2,
  density: { letters: 10, cells: 3, dots: 1, bars: 1 / 64 },
  densityUnits: { letters: 'css', cells: 'device', dots: 'device', bars: 'device' },
  letterSizes: [9, 10, 11, 12, 13, 14],
  monoAdvance: 0.6,
  labelFont: '12px sans',
  itemFont: '10px sans',
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
    highlightTint: 0.3,
    projected: 0.5,
    hoverOffset: 1,
    cursorInk: 'msui-text',
  },
  axis: { directionMark: 'label' },
  margin: 0.5,
  followArc: 'shortest',
  maxTemplateTracks: 200,
  nodes: { holdShown: true },
  ranges: { overflow: 'extra-lanes' },
  values: { aggregate: 'minmax', coarser: 1, direct: true },
  ruler: { labelGap: 48 },
  rings: {
    hole: 0.35,
    maxDiameter: 720,
    blockGap: 1,
    letters: 'upright',
    center: 'title',
    labelGap: 8,
    leaderBend: 12,
    window: 'followed',
  },
  internedStrings: DefaultMaxStrings,
  imageBytes: 64 * 1024 * 1024,
  rows: { labelMinPx: 8, pinInPlace: true },
  navigator: {
    window: { minPx: 16, hitPx: 12, touchHitPx: 24, handlePx: 5, floor: 'centred' },
    trackScale: 0.25,
    trackMinPx: 2,
    marks: ['selection', 'highlights', 'cursor'],
    rail: { grownPx: 6, grabShare: 0.4 },
  },
  text: {
    scale: 12,
    groupGap: 0.5,
    lineGap: 8,
    bandGap: 2,
    collapseEmpty: true,
    maxLines: 100_000,
    maxHeight: 600,
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
    heights: { ...DefaultRenderConfig.heights, ...over.heights },
    density: { ...DefaultRenderConfig.density, ...over.density },
    densityUnits: { ...DefaultRenderConfig.densityUnits, ...over.densityUnits },
    ranges: { ...DefaultRenderConfig.ranges, ...over.ranges },
    values: { ...DefaultRenderConfig.values, ...over.values },
    ruler: { ...DefaultRenderConfig.ruler, ...over.ruler },
    rings: { ...DefaultRenderConfig.rings, ...over.rings },
    outline: { ...DefaultRenderConfig.outline, ...over.outline },
    emphasis: { ...DefaultRenderConfig.emphasis, ...over.emphasis },
    text: { ...DefaultRenderConfig.text, ...over.text },
    rows: { ...DefaultRenderConfig.rows, ...over.rows },
    navigator: {
      ...DefaultRenderConfig.navigator,
      ...over.navigator,
      window: { ...DefaultRenderConfig.navigator.window, ...over.navigator?.window },
      rail: { ...DefaultRenderConfig.navigator.rail, ...over.navigator?.rail },
    },
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

/** The size of a font's shorthand in px (12 when it names none). */
export function fontSize(font: string): number {
  const m = /(\d+(?:\.\d+)?)px/.exec(font);
  return m ? Number(m[1]) : 12;
}

/** Whether a step draws a mark per position (per cell), not per device column. */
export function perPosition(step: DensityStep): boolean {
  return step === 'letters' || step === 'cells' || step === 'dots';
}

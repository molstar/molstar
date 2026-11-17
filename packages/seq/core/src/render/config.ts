/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout's settings (_spec/README.md "Render", Configuration): every number is a setting with its recommended
// default; the density steps and their thresholds, shared by every layer type.

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
  readonly heights: { readonly letters: number; readonly ranges: number; readonly values: number };
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
  /** What a want adds on each side of the view, as a fraction of it (the store's `margin` priority). */
  readonly margin: number;
  /** A follower's leader pieces carried as their image's hull, or the shortest arc on a cyclic follower. */
  readonly followArc: 'shortest' | 'hull';
  /** Tracks one template expands into at most; more is `seq.render.template.limit` and the rest are not laid out. */
  readonly maxTemplateTracks: number;
  /** Ranges: items past the lanes collapse into a density row at the bottom (`extra-lanes`), or the whole layer does. */
  readonly ranges: { readonly overflow: 'extra-lanes' | 'layer' };
  /** Values below a device pixel a position: each column's extent (`minmax`, negative extremes show) or its maximum. */
  readonly values: { readonly aggregate: 'minmax' | 'max' };
}

export const DefaultRenderConfig: RenderConfig = {
  inset: 8,
  labelWidth: 112,
  labelGap: 8,
  trackGap: 4,
  sectionGap: 8,
  heights: { letters: 20, ranges: 36, values: 32 },
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
  outline: { hover: 1.5, selection: 1.5, cursor: 2.5 },
  margin: 0.5,
  followArc: 'shortest',
  maxTemplateTracks: 200,
  ranges: { overflow: 'extra-lanes' },
  values: { aggregate: 'minmax' },
};

/** A config with some fields overridden (the nested ones whole or field by field). */
export function renderConfig(over: Partial<RenderConfig> = {}): RenderConfig {
  return {
    ...DefaultRenderConfig,
    ...over,
    heights: { ...DefaultRenderConfig.heights, ...over.heights },
    density: { ...DefaultRenderConfig.density, ...over.density },
    densityUnits: { ...DefaultRenderConfig.densityUnits, ...over.densityUnits },
    ranges: { ...DefaultRenderConfig.ranges, ...over.ranges },
    values: { ...DefaultRenderConfig.values, ...over.values },
    outline: { ...DefaultRenderConfig.outline, ...over.outline },
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

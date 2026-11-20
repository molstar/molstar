/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RenderConfig } from './config.js';

// A narrow frame (_spec/README.md "Render", A narrow frame): what gives way, in order, when the plot (the data's width)
// runs short, keyed on the plot's width, never the frame's: the vertical overview narrows toward its floor and then
// hides (the rows rail stays), then the label column shrinks to its floor and the tree collapses, then the plot keeps
// its floor and the frame scrolls. The order is by what is cheapest to lose: a miniature before names, names before data.

/** The columns a view lays out with in a frame of its width, and what the frame needs. */
export interface NarrowColumns {
  /** The label column's width (its own, else shrunk toward its floor). */
  readonly labelWidth: number;
  /** The tree's width (its own, else 0: collapsed). */
  readonly tree: number;
  /** The plot's width with every column at its own width (what a navigator beside the frame reads). */
  readonly plotPx: number;
  /** Everything but the plot at the columns' own widths, CSS px. */
  readonly columnsPx: number;
  /** The least frame width that keeps the plot at its floor, the columns at theirs: a narrower frame scrolls. */
  readonly minWidth: number;
}

/**
 * The label column and the tree for a frame `width` px wide whose columns take `fixed` px (the label column
 * `labelWidth` of it, the tree `tree`; row columns as granted at this width): unchanged while the plot keeps
 * `navigator.overview.minPlotPx`; else the label column gives up to its floor (`labelColumn.min`, or its own width if
 * less), then the tree all of its. `columnsPx`: the columns at their own widths, whatever the frame (`fixed` when
 * absent): what `plotPx` and a navigator beside the frame read, so a width they decide never feeds back.
 */
export function narrowColumns(
  width: number,
  fixed: number,
  labelWidth: number,
  tree: number,
  config: RenderConfig,
  columnsPx = fixed,
): NarrowColumns {
  const minPlot = config.navigator.overview.minPlotPx;
  const plotPx = width - columnsPx;
  const floor = Math.min(labelWidth, config.labelColumn.min);
  const minWidth = fixed - (labelWidth - floor) - tree + minPlot;
  let lack = minPlot - (width - fixed);
  if (!(lack > 0)) return { labelWidth, tree, plotPx, columnsPx, minWidth };
  const label = Math.max(floor, Math.floor(labelWidth - lack));
  lack -= labelWidth - label;
  return { labelWidth: label, tree: lack > 0 ? 0 : tree, plotPx, columnsPx, minWidth };
}

/**
 * A vertical overview's width beside a frame whose plot would be `plotPx` wide without it: its own `size` while the
 * plot keeps `navigator.overview.minPlotPx` beside it; narrowed toward `navigator.overview.minPx` to keep it; 0
 * (hidden: the rows rail stays the way to move) when even at its floor the plot would be under it.
 */
export function overviewWidth(plotPx: number, size: number, config: RenderConfig): number {
  const { minPlotPx, minPx } = config.navigator.overview;
  if (!(minPlotPx > 0)) return size;
  const floor = Math.min(size, minPx);
  if (plotPx - floor < minPlotPx) return 0;
  return Math.max(floor, Math.min(size, Math.floor(plotPx - minPlotPx)));
}

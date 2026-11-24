/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The example story's grids. A grid's column is `minmax(0, 1fr)`, never the implicit `auto`: an `auto` column grows to
// its widest item's content, and the canvas (measured from its container) would then widen the page past its right edge,
// the vertical overview beside it clipped. Every item has `min-width: 0` for the same reason.

const Column = 'grid-template-columns: minmax(0, 1fr);';

/** A titled section's box. */
export const Box = `display: grid; ${Column} gap: var(--msui-space-2); min-width: 0; align-content: start;`;

/** A find box's result block: lines kept as said, a long one wrapped inside the box (never cut at its right edge). */
export const FoundBlock =
  'margin: 0; padding: var(--msui-space-3); overflow: auto; white-space: pre-wrap; overflow-wrap: anywhere; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';

/** The page: its one column the page's width. */
export const PageGrid = `display: grid; ${Column} gap: var(--msui-space-5); padding: var(--msui-space-5); min-height: 100vh; min-width: 0;`;

/** The canvas and the SVG: side by side where each has 420 px, else (or `stacked`) one under the other. */
export function pairGrid(stacked: boolean): string {
  const columns = stacked ? 'minmax(0, 1fr)' : 'repeat(auto-fit, minmax(min(420px, 100%), 1fr))';
  return `display: grid; grid-template-columns: ${columns}; gap: var(--msui-space-5); min-width: 0;`;
}

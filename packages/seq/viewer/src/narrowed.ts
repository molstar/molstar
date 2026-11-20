/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rowsRailTracks } from '@molstar/seqx-core/render';
import type { ResolvedNavigator, ResolvedView } from '@molstar/seqx-core/view';
import type { ViewerConfig } from './config.js';

// The navigators a narrow frame hid (core's "Render", A narrow frame; _spec/README.md "Composition", Navigators): a
// vertical overview at 0 px for want of width is hidden as a hidden one is, so under `navigator.rowsRail:
// unlessOverview` its rows rail comes back: the way to move is never lost. Kept by viewer (the models say it, the
// viewer's rails read it), with a version the rails' key holds.

interface Narrowed {
  readonly ids: Set<string>;
  version: number;
}

const ByViewer = new WeakMap<object, Narrowed>();

function of(viewer: object): Narrowed {
  let n = ByViewer.get(viewer);
  if (!n) ByViewer.set(viewer, (n = { ids: new Set(), version: 0 }));
  return n;
}

/** Says whether a narrow frame hides navigator `id` (its width 0 for want of room); true when that changed. */
export function setNarrowed(viewer: object, id: string, hidden: boolean): boolean {
  const n = of(viewer);
  if (n.ids.has(id) === hidden) return false;
  if (hidden) n.ids.add(id);
  else n.ids.delete(id);
  n.version++;
  return true;
}

/** The version of the navigators a narrow frame hides: it moves when one is hidden or shown again. */
export function narrowedVersion(viewer: object): number {
  return ByViewer.get(viewer)?.version ?? 0;
}

/**
 * The rows tracks that draw a rows rail (`rowsRailTracks` by `navigator.rowsRail`), a navigator a narrow frame hides
 * counted as hidden.
 */
export function rowsRailsOf(
  viewer: object,
  view: ResolvedView,
  hosted: readonly ResolvedNavigator[],
  config: ViewerConfig,
): string[] {
  const ids = ByViewer.get(viewer)?.ids;
  const out = (list: readonly ResolvedNavigator[]) => list.filter((n) => !ids?.has(n.id) && !ids?.has(n.address));
  const shown = ids?.size ? { ...view, navigators: out(view.navigators) } : view;
  return rowsRailTracks(shown, ids?.size ? out(hosted) : hosted, config.navigator.rowsRail);
}

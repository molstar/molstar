/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GridLay, LayerLayout, RowsLay, SectionLayout, ViewLayout } from '@molstar/seqx-core/render';
import { trackOfLayer } from './rows.js';

// A section's layer of two axes (_spec/README.md "Behaviors", Cursor): its rows layer or its grid, which the keyboard
// cursor's row is on. Internal: not in the behaviors' entry.

/**
 * The layer of two axes of a section in the last frame: the one showing `ySystem` (a rows layer's rows system, a
 * grid's y), else (no `ySystem`) its rows layer; undefined when it has none.
 */
export function twoAxesLayer(section: SectionLayout, ySystem?: string): LayerLayout | undefined {
  const is = (l: LayerLayout) =>
    l.type.name === 'rows'
      ? ySystem === undefined || (l.lay as RowsLay).rowsSystem === ySystem
      : l.type.name === 'grid' && ySystem !== undefined && (l.lay as GridLay).ySystem === ySystem;
  for (const t of section.tracks) for (const l of t.layers) if (is(l)) return l;
  for (const ls of section.lays) for (const l of ls) if (is(l)) return l;
  for (const l of section.sectionLays) if (is(l)) return l;
  return undefined;
}

/** A section's rows layer as the last frame laid it out, and its track (where its band's state is kept). */
export interface RowsOf {
  readonly lay: RowsLay;
  readonly track: string;
}

/** The rows layer of the section at `address` (one a section: the cursor has one row order), from the last frame. */
export function rowsOf(layout: ViewLayout | undefined, address: string): RowsOf | undefined {
  const s = layout?.sections.find((x) => x.section.address === address);
  const l = s && twoAxesLayer(s);
  return l ? { lay: l.lay as RowsLay, track: trackOfLayer(l.key) } : undefined;
}

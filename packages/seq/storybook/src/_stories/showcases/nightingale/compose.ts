/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import {
  readSection,
  type SectionRead,
  type ThemeName,
  themeColors,
  toFeatures,
  toHighlight,
  toLineData,
} from './nightingale-adapter.js';

// The Nightingale showcase's composition, DOM-free: the story's root and its headless twin call these, so the twin
// feeds Nightingale's shapes from the viewer as the story does (the viewer's resolved view and store, its hover and
// selection), and the root only sets them on the elements.

/** The positions the sequence and features story shows, 1-based inclusive: Nightingale draws letters from about 8 px a residue. */
export const Window = { start: 60, end: 130 } as const;

/** What Nightingale's elements are set from: the section as the viewer holds it, in a theme's colors. */
export interface NightingaleData {
  readonly read: SectionRead;
  readonly features: ReturnType<typeof toFeatures>;
  readonly lines: ReturnType<typeof toLineData>;
  readonly colors: ReturnType<typeof themeColors>;
}

/** Reads the viewer's first section (its resolved view, its model, its store) for Nightingale. */
export async function nightingaleData(viewer: SeqViewer, theme: ThemeName): Promise<NightingaleData> {
  const model = viewer.view.model.value;
  const view = viewer.view.resolved.value;
  if (!view) throw new Error('the viewer shows no resolved view');
  const read = await readSection(model, view, viewer.store.store, theme);
  const colors = themeColors(theme);
  return {
    read,
    features: toFeatures(read.ranges),
    lines: toLineData(read.values, 'Conservation', read.valueColor, 'none'),
    colors,
  };
}

/**
 * The viewer's hover and selection on the section's system as Nightingale's `highlight`: selection first, then the
 * hover. The selection's `shown` (what is drawn: a clicked domain is an element, its extent only in `shown`).
 */
export function highlightOf(viewer: SeqViewer, system: string): string {
  return [toHighlight(viewer.selection.value.value.shown, system), toHighlight(viewer.hover.location(), system)]
    .filter(Boolean)
    .join(',');
}

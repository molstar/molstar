/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay, ViewLayout } from '@molstar/seqx-core/render';
import { trackOfLayer } from '@molstar/seqx-viewer/behaviors/two-axes';
import type { ViewModel } from '@molstar/seqx-viewer/models';

// A track's lanes scrolling inside it (core's _spec/README.md "Render", Track heights): the wheel's vertical delta over
// a band whose lanes overflow it scrolls them (the band's `bandTop`, in lanes, session state); at the top or bottom, or
// with every lane in view, it is given back (the section, then the page: scroll chaining).

/** The scrolling ranges band of a track in a layout (its lanes overflow the band), if any. */
export function lanesBandOf(layout: ViewLayout | undefined, track: string): RangesLay | undefined {
  for (const s of layout?.sections ?? [])
    for (const band of s.lays)
      for (const l of band ?? []) {
        if (l.type.name !== 'ranges' || trackOfLayer(l.key) !== track) continue;
        const lay = l.lay as RangesLay;
        if (lay.scrolls && lay.scrollMax > 0) return lay;
      }
  return undefined;
}

/**
 * Scrolls a track's lanes by `dy` px: true when they moved. A band at its end, or a track with no scroll, takes
 * nothing (false: the wheel goes on to the section and the page).
 */
export function scrollLanes(view: ViewModel, layout: ViewLayout | undefined, track: string, dy: number): boolean {
  const lay = lanesBandOf(layout, track);
  if (!lay) return false;
  const top = Math.max(0, Math.min(lay.scrollMax, lay.top + dy / lay.laneHeight));
  if (top === lay.top) return false;
  view.setBand(track, { bandTop: top });
  return true;
}

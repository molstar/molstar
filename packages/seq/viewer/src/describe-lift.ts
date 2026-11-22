/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An item's images through a map read by region (_spec/README.md "Hover card": a liftover's places): the card of an
// item with several (a duplicated stretch, EC-20) lists them, so the reader knows the projected hover lights each.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { spanText } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { type RegionProjector, regionMaps, type ResolvedView } from '@molstar/seqx-core/view';

/** The most images a card names: past it, the count alone. */
const MaxNamed = 4;

/**
 * `[start, end)` on `system` carried through the maps read by region onto every block the view shows: its images as
 * `2 places: chr1 62,001–68,000, chr2 62,001–68,000` when there are at least two, else undefined (one image is the
 * projected hover itself, nothing to list).
 */
export function liftPlaces(
  model: StaticModel,
  view: ResolvedView,
  regions: RegionProjector,
  system: string,
  start: number,
  end: number,
): string | undefined {
  if (!regionMaps(model.spec).length) return undefined;
  const named: string[] = [];
  const seen = new Set<string>([system]);
  const buf: number[] = [];
  let n = 0;
  for (const s of view.sections)
    for (const to of s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system]) {
      if (seen.has(to)) continue;
      seen.add(to);
      buf.length = 0;
      if (!regions.project(model, system, [start, end], to, buf) || !buf.length) continue;
      const r = HalfOpenRanges.normalize(buf);
      const name = model.systems.get(to)?.name ?? to;
      for (let k = 0; k + 1 < r.length; k += 2, n++)
        if (named.length < MaxNamed) named.push(`${name} ${spanText(r[k], r[k + 1])}`);
    }
  if (n < 2) return undefined;
  return `${n} places: ${named.join(', ')}${n > named.length ? ', …' : ''}`;
}

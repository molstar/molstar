/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { partOf, rangesLocation } from '@molstar/seqx-core/query';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { Snap } from './selection.js';
import type { ViewModel } from './view.js';

// A nucleotide's codon for a drag (EC-149; _spec/README.md "Models", Selection): the span of the codon through the
// view's codon maps and the holes inside it (the intron a split codon spans), by position. Called on every pointer
// move, so a unit is computed once per position and kept (`null`: none) while the model and the resolved view are the
// same: a warm move allocates nothing. The kept units are bounded (`MaxUnits`, then dropped).

const MaxUnits = 4096;

/** `snap` for `SelectionOptions`, over the view (read late: it is made after the selection). */
export function snapOfCodons(view: () => ViewModel): (system: string, position: number) => Snap | undefined {
  let model: unknown;
  let resolved: ResolvedView | undefined;
  let systems: string[] = [];
  const units = new Map<string, Map<number, Snap | null>>();
  return (system, position) => {
    const v = view();
    const m = v.model.value;
    const r = v.resolved.value;
    if (m !== model || r !== resolved) {
      model = m;
      resolved = r;
      systems = r?.sections.map((s) => s.axis.system) ?? [];
      units.clear();
    }
    let by = units.get(system);
    if (!by) units.set(system, (by = new Map()));
    const kept = by.get(position);
    if (kept !== undefined) return kept ?? undefined;
    if (by.size >= MaxUnits) by.clear();
    const closed = m.close(rangesLocation(system, [position, position + 1]), systems);
    const ranges = partOf(closed, system)?.ranges;
    let unit: Snap | null = null;
    if (ranges?.length && !(ranges.length === 2 && ranges[1] - ranges[0] === 1)) {
      // the pieces are the codon's own; what lies between them is a hole (an intron), never selected
      const holes: number[] = [];
      for (let i = 1; i + 1 < ranges.length; i += 2) holes.push(ranges[i], ranges[i + 1]);
      unit = { start: ranges[0], end: ranges[ranges.length - 1], ...(holes.length ? { holes } : {}) };
    }
    by.set(position, unit);
    return unit ?? undefined;
  };
}

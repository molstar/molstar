/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GrapeSection, PeachSection } from '@molstar/seqx-examples/across/views';
import type { Blocks } from '@molstar/seqx-examples/across/genomes';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { expect } from 'vitest';

// What the across-sections tests share: a system's name, whether a section shows a place, and the drawn-once check.

export const bridges = (v: SeqViewer) => v.frame!.layout.bridges?.list ?? [];

export const named = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)?.name;

/** Whether a section's window shows a position of a chromosome (a block's place over its window). */
export function shows(v: SeqViewer, section: string, chr: string, p: number): boolean {
  const s = v.frame!.layout.sections.find((x) => x.section.address === section)!;
  const strips = s.placement.strips;
  for (let i = 0; i < strips.n; i++) {
    const axis = strips.axes[i];
    if (named(v, axis.system) !== chr) continue;
    const vis = axis.visible([]);
    for (let j = 0; j + 1 < vis.length; j += 2) if (p >= vis[j] && p < vis[j + 1]) return true;
  }
  return false;
}

/**
 * Every block anchored in either side's window drawn exactly once (its item, or a bin pair of its strand holding both
 * its starts), none twice (EC-209).
 */
export function expectDrawnOnce(v: SeqViewer, bl: Blocks) {
  const ribbons = bridges(v)[0].ribbons();
  const count = new Array<number>(bl.a.length).fill(0);
  for (const r of ribbons) {
    const a = named(v, r.aSystem);
    const b = named(v, r.bSystem);
    for (let j = 0; j < bl.a.length; j++) {
      if (bl.a[j] !== a || bl.b[j] !== b || r.reversed !== bl.strand[j] < 0) continue;
      const as = bl.aStart[j];
      const bs = bl.bStart[j];
      if (
        r.bin === undefined
          ? r.aStart === as && r.bStart === bs
          : as >= r.aStart && as < r.aEnd && bs >= r.bStart && bs < r.bEnd
      )
        count[j]++;
    }
  }
  const wrong: string[] = [];
  for (let j = 0; j < bl.a.length; j++) {
    const shown = shows(v, PeachSection, bl.a[j], bl.aStart[j]) || shows(v, GrapeSection, bl.b[j], bl.bStart[j]);
    if (count[j] > 1 || (shown && count[j] !== 1))
      wrong.push(`${bl.a[j]}:${bl.aStart[j]} ${bl.b[j]}:${bl.bStart[j]} ×${count[j]}`);
  }
  expect(wrong.slice(0, 5)).toEqual([]);
}


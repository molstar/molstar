/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Layouts story's composition, shared by the story and its headless twin: a tracks section past its cap (a ruler,
// a lanes band of 40 lanes, 30 gene tracks), at the default cap (`section.maxPx`, 600), a host's 300 and a spec
// `height` of 200.

import type { Node } from '@molstar/core/repr/tree';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** The section's address. */
export const Genes = '/views/v/genes';

/** A cap of the story: its title, what it says, the host's `section.maxPx` or the spec's `height`. */
export interface CapEntry {
  readonly id: string;
  readonly title: string;
  readonly maxPx?: number;
  readonly height?: number;
}

export const Caps: readonly CapEntry[] = [
  { id: 'default', title: 'The default cap (`section.maxPx`: 600 px)' },
  { id: 'host', title: 'A host’s `section.maxPx`: 300 px', maxPx: 300 },
  { id: 'spec', title: 'The section’s own `height`: 200 px', height: 200 },
];

/** The spec: a ruler, a lanes band of 40 overlapping transcripts and 30 genes of one feature each. */
export function capSpec(cap: CapEntry): Node {
  const b = createSeqSpec({ title: 'A tracks section past its cap' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTTGCA'.repeat(250) });
  const tall = b.data.ranges('transcripts', {
    target: s,
    rows: {
      start: Array.from({ length: 40 }, (_, i) => 40 * i),
      end: Array.from({ length: 40 }, (_, i) => 900 + 25 * i),
      id: Array.from({ length: 40 }, (_, i) => `T-${201 + i}`),
    },
  });
  b.views.view('v', (v) => {
    const genes = v.section('genes', { of: s, ...(cap.height !== undefined ? { height: cap.height } : {}) });
    genes.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    genes.track('transcripts', { label: 'Transcripts' }).layer('l', { type: 'ranges', data: tall });
    for (let i = 0; i < 30; i++) {
      const start = (i * 61) % 1400;
      const d = b.data.ranges(`g${i}`, {
        target: s,
        rows: { start: [start], end: [start + 300 + ((i * 37) % 300)], id: [`GENE${i + 1}`] },
      });
      genes.track(`g${i}`, { label: `Gene ${i + 1}` }).layer('l', { type: 'ranges', data: d });
    }
  });
  return b.build();
}

/** The host's config of a cap (none for the default and a spec's `height`). */
export function capConfig(cap: CapEntry) {
  return cap.maxPx !== undefined ? [[SeqCoreConfig.render, { section: { maxPx: cap.maxPx } }] as const] : [];
}

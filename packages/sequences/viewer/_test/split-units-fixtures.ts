/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A split track whose rows sit on units of their own (HR1): a reference protein, a chain sequence covering its residues
// 10–49 (a mapping) and three copies of it as units A, B and C, each observed in part; the chains' observed stretches
// as one `ranges` annotation on the units (a `target` column), drawn on the reference's block split by chain. Shared by
// the tests, the heap guard and the Storybook story; no I/O.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

export const REFERENCE = 'MSTEQKLISEEDLNGAVKRLHEFAVPWQDTRYCMNPGSTEELKVAAQYIDGRLHHPSEK';
/** The chain: the reference's residues 10–49 (0-based [10, 50)). */
export const CHAIN_AT = 10;
export const CHAIN = REFERENCE.slice(CHAIN_AT, 50);
/** The units and their observed stretches on the chain (0-based, stop exclusive). */
export const UNITS = [
  { ref: 'A', observed: [[0, 40]] },
  {
    ref: 'B',
    observed: [
      [2, 18],
      [24, 38],
    ],
  },
  { ref: 'C', observed: [[5, 35]] },
] as const;

/**
 * The document. `highlights`: chain B's residues 4–7 (a highlight: B's row and the reference's tracks only) and the
 * reference's 30–33 (an accent: every row).
 */
export function splitUnitsDoc(o: { highlights?: boolean; title?: string } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title ?? 'Per-row marks' });
  const ref = b.data.sequence({ ref: 'ref', alphabet: 'protein', value: REFERENCE, label: 'Reference' } as never);
  b.data.sequence({ ref: 'chain', alphabet: 'protein', value: CHAIN } as never);
  for (const u of UNITS)
    b.data.unit({ ref: u.ref, sequence: 'chain', observed: u.observed, label: `chain ${u.ref}` } as never);
  b.data.mapping({
    ref: 'chain-ref',
    from: 'chain',
    to: 'ref',
    method: 'segments',
    segments: [[0, CHAIN_AT, CHAIN.length]],
  } as never);
  const rows = UNITS.flatMap((u) => u.observed.map(([s, t]) => ({ unit: u.ref, s, t })));
  const chains = b.data.annotation({
    ref: 'chains',
    geometry: 'ranges',
    label: 'Observed',
    value: {
      target: rows.map((r) => r.unit),
      start: rows.map((r) => r.s),
      stop: rows.map((r) => r.t),
      id: rows.map((r) => `chain ${r.unit}`),
      chain: rows.map((r) => r.unit),
    },
  } as never);
  const marks = [
    b.data.annotation({ ref: 'on-b', geometry: 'ranges', target: 'B', value: { start: [4], stop: [8] } } as never),
    b.data.annotation({ ref: 'on-ref', geometry: 'ranges', target: ref, value: { start: [30], stop: [34] } } as never),
  ];
  const v = b.view({ ref: 'v', label: 'V' });
  const sec = v.section({
    ref: 'v.main',
    ...(o.highlights
      ? {
          highlights: [
            { data: marks[0], style: 'highlight' },
            { data: marks[1], style: 'accent' },
          ],
        }
      : {}),
  } as never);
  const block = sec.block({ sequence: ref, label: 'Reference' } as never);
  block.track({ ref: 't-letters', label: 'Reference' } as never).letters();
  block.track({ ref: 't-chains', label: 'Chains', expanded: true } as never).ranges({
    data: chains,
    display: 'union',
    split: 'chain',
  } as never);
  return b.getDocument();
}

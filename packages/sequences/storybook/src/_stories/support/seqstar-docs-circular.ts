/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The circular genome's document (circular-maps.stories.tsx), and what the plasmid map (case 16, the Viewer cases
// story) does with its buttons: rotating the ring's origin and zooming its linear detail. Kept out of the story files so
// the smoke test can open and drive each on a headless viewer. Synthetic sequences: no fixture, no network.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { Patch, SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { SeqViewerContext } from '@molstar/seq-viewer';

/** A deterministic synthetic DNA of `n` nucleotides (a linear congruential walk over ACGT). */
function synthetic(n: number, seed: number): string {
  let x = seed;
  let out = '';
  for (let i = 0; i < n; i++) {
    x = (x * 1103515245 + 12345) & 0x7fffffff;
    out += 'ACGT'[(x >> 16) & 3];
  }
  return out;
}

/** Case 16's ring section, and its block (generated: the view's first), whose `origin` a rotation patches. */
export const PLASMID_RING = 'pmap.ring';
const PLASMID_BLOCK = 'pmap.block';
/** pSyn's length. */
export const PSYN_LENGTH = 2686;
/** The ring started from `origin` (0-based): "Rotate the origin". */
export const rotateRing = (origin: number): Patch => ({
  format: 'seq_view_spec_patch',
  ops: [{ op: 'set', ref: PLASMID_BLOCK, params: { origin } }],
});

/**
 * The plasmid map's linear detail zoomed to the multiple cloning site (360 … 490), animated: the ring's window (a
 * sector with grips, CG3) shows it. False: no detail.
 */
export function zoomDetail(ctx: SeqViewerContext): boolean {
  const ref = ctx.viewports.sections.value.find((r) => r !== PLASMID_RING);
  const vp = ref === undefined ? null : ctx.viewports.get(ref);
  vp?.zoomTo(0, 360, 490, { animate: true });
  return vp !== null;
}

/**
 * A synthetic 120,000-nucleotide circular genome: 160 genes on both strands and its GC skew as bars, through
 * `b.preset.circularGenome` (two strand rings, a values ring, and a linear detail).
 */
export function circularGenomeDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'A circular genome' });
  const n = 120_000;
  const g = b.data.sequence({ ref: 'chr', label: 'chr', alphabet: 'dna', cyclic: true, value: synthetic(n, 7) });
  const genes = Array.from({ length: 160 }, (_, i) => {
    const start = Math.floor((i * n) / 160) + ((i * 37) % 200);
    return [start, start + 300 + ((i * 131) % 400), i % 3 === 0 ? '-' : '+', `g${i + 1}`] as const;
  });
  const features = b.data.annotation({
    ref: 'genes',
    geometry: 'ranges',
    target: g,
    value: {
      start: genes.map((x) => x[0]),
      stop: genes.map((x) => x[1]),
      strand: genes.map((x) => x[2]),
      label: genes.map((x) => x[3]),
      category: genes.map((x) => (x[2] === '+' ? 'forward' : 'reverse')),
    },
  } as never);
  const bins = 600;
  const skew = b.data.annotation({
    ref: 'skew',
    label: 'GC skew',
    geometry: 'values',
    target: g,
    value: {
      start: Array.from({ length: bins }, (_, i) => (i * n) / bins),
      stop: Array.from({ length: bins }, (_, i) => ((i + 1) * n) / bins),
      // A skew that turns at the origin and its opposite, with a little noise.
      value: Array.from({ length: bins }, (_, i) =>
        Number((Math.cos((2 * Math.PI * i) / bins) * 0.08 + 0.02 * Math.sin(i * 1.7)).toFixed(3)),
      ),
    },
  } as never);
  const v = b.view({ ref: 'genome', label: 'Genome' });
  b.preset.circularGenome(v, g, { ref: 'genome.ring', label: 'chr', features, values: [skew] });
  return b.getDocument();
}

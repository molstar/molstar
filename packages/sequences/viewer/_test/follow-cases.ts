/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the `follow` behavior's guards (LK1; seqspec:viewer#performance; follow.heap.ts) and their quick twin
// (follow.quick.ts): a leader with two followers (a protein through the codon step and the gene reversed), driven at
// input rate by a drag and by a Ctrl wheel zoom; each step moves all three. Fractional coordinates, as an engine's;
// every default behavior subscribed. The heap guard and the quick twin share them.

import { createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { expect } from 'vitest';
import { Inputs, Modifier, PointerType } from '../src/index.js';
import { driven } from './input-fixtures.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;
let sink = 0;

/** Case 14's gene, a second copy of its two exons' worth of codons, and three sections: gene, protein, gene reversed. */
function doc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const coding = 'ATGGCTTCCAAGGAACTGACCCCACGTGACGGCCATTGGTACTTCAACAGCCAGGTGGAG'.repeat(20);
  const gene = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: `GAATTCGCCACC${coding}TAAGCGGCC` });
  const protein = translate(b.data, { ref: 'protein', from: gene, cds: [[12, 12 + coding.length]] });
  const v = b.view({ ref: 'v', label: 'V' });
  v.section({ ref: 'v.dna' }).block({ sequence: gene }).track().letters();
  v.section({ ref: 'v.protein', follow: 'v.dna' }).block({ sequence: protein }).track().letters();
  v.section({ ref: 'v.rev', follow: 'v.dna', direction: -1 }).block({ sequence: gene }).track().letters();
  return b.getDocument();
}

/** The three viewports of the document, opened. */
async function opened() {
  const t = await driven(doc(), { labels: false });
  const [dna, protein, rev] = ['v.dna', 'v.protein', 'v.rev'].map((r) => t.vp(r));
  return { t, dna, protein, rev };
}

/** A drag on the leader, started: its step moves the pointer; `finish` releases it. */
export async function dragging() {
  const { t, dna, protein, rev } = await opened();
  dna.zoomTo(0, 300, 360);
  const y = t.layout('v.dna').scrolledTop + 8.5;
  const e = Inputs.pointer();
  e.section = t.ctx.layout.sectionIndex('v.dna');
  t.pointer('down', 'v.dna', 500.5, y);
  t.pointer('move', 'v.dna', 490.5, y);
  expect(t.ctx.ui.drag.value?.owner).toBe('dragPan');
  e.type = PointerType.move;
  e.buttons = 1;
  e.button = -1;
  e.y = y;
  const step = (i: number) => {
    // Whole px: the leader's scroll snaps to device px, so a fractional step would drift it to an end over the run.
    e.x = 490.5 + (i & 15) * 3;
    t.ctx.input.pointer(e);
  };
  return { step, dna, protein, rev, finish: () => t.pointer('up', 'v.dna', 500.5, y) };
}

/** A Ctrl wheel zoom on the leader at x 400.5, in and out by turns: its step and the three viewports. */
export async function wheeling() {
  const { t, dna, protein, rev } = await opened();
  dna.zoomTo(0, 300, 400);
  const w = Inputs.wheel();
  w.section = t.ctx.layout.sectionIndex('v.dna');
  w.y = t.layout('v.dna').scrolledTop + 8.5;
  w.modifiers = Modifier.ctrl;
  w.deltaX = 0;
  const step = (i: number) => {
    w.x = 400.5 + (i & 7);
    w.deltaY = i & 1 ? -3.25 : 3.25;
    sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
  };
  return { step, dna, protein, rev };
}

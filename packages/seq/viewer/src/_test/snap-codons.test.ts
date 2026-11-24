/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { snapOfCodons } from '../models/snap-codons.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { proteinOnGene } from './fixture.js';
import { type Frames, installFrames } from './frames.js';

// A drag over nucleotides snaps to whole codons, a codon split by an intron its pieces (EC-149; _spec/README.md
// "Models", Selection): the fixture's exons [2, 7) and [12, 19) code the codons (2 3 4), (5 6 12), (13 14 15),
// (16 17 18); the intron [7, 12) is no codon's.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(proteinOnGene());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const at = (v: SeqViewer, position: number): [number, number] => {
  const s = v.frame!.layout.sections[1];
  const box = s.tracks[0].box;
  return [(s.axis!.x(position) + s.axis!.x(position + 1)) / 2, box.y + box.height / 2];
};
const gene = (v: SeqViewer) => [
  ...(v.selection.value.value.shown.parts.find((p) => p.system === '/data/g')?.ranges ?? []),
];

describe('a drag over a split codon', () => {
  it('EC-149: from a split codon’s anchor, the selection is its three nucleotides and the dragged-over ones, never the intron', async () => {
    const v = await make();
    const [x, y] = at(v, 6);
    v.input.down(x, y);
    // alone: the codon (5 6 12), its pieces
    v.input.move(at(v, 5)[0], y);
    expect(gene(v)).toEqual([5, 7, 12, 13]);
    // over the next codon (13 14 15): the pieces and it whole, the intron out
    v.input.move(at(v, 14)[0], y);
    expect(gene(v)).toEqual([5, 7, 12, 16]);
    v.input.up();
    expect(gene(v)).toEqual([5, 7, 12, 16]);
  });

  it('EC-149: ending on a split codon takes its pieces too, the intron out', async () => {
    const v = await make();
    const [x, y] = at(v, 14);
    v.input.down(x, y);
    v.input.move(at(v, 12)[0], y);
    v.input.up();
    expect(gene(v)).toEqual([5, 7, 12, 16]);
  });
});

describe('the snap per pointer move', () => {
  it('EC-149: a position is computed once while the model and the view stand; a warm move does no work and allocates nothing', async () => {
    const v = await make();
    const closed = vi.spyOn(v.view.model.value, 'close');
    const snap = snapOfCodons(() => v.view);
    const first = snap('/data/g', 6);
    expect(first).toEqual({ start: 5, end: 13, holes: [7, 12] });
    expect(closed).toHaveBeenCalledTimes(1);
    // warm: the same unit, by identity, with no further closure
    for (let i = 0; i < 1000; i++) expect(snap('/data/g', 6)).toBe(first);
    expect(closed).toHaveBeenCalledTimes(1);
    // a protein's position has none: kept as none
    expect(snap('/data/p', 1)).toBeUndefined();
    const n = closed.mock.calls.length;
    for (let i = 0; i < 1000; i++) expect(snap('/data/p', 1)).toBeUndefined();
    expect(closed).toHaveBeenCalledTimes(n);
  });
});

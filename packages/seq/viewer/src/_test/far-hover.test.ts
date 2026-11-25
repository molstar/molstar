/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { FakeCanvas } from '../canvas/_test/recording.js';
import { CanvasBackend } from '../canvas/backend.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The hover at a far zoom (_spec/README.md "Emphasis"): a cell narrower than a device pixel still has its outline, so a
// pointer moving across a section fitted to thousands of positions hovers continuously, on the leader and the follower.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

/** A protein of 1,200 residues on a gene of 4,000 nucleotides (one exon from 200), the gene following the protein. */
function longProteinOnGene() {
  const b = createSeqSpec({ title: 'far' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: `M${'EA'.repeat(600)}`.slice(0, 1200) });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(1000) });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [200], end: [3800] }, strand: '-' });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    const gene = v.section('gene', { of: g, follow: protein });
    gene.track('letters').layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('the hover at a far zoom', () => {
  it('outlines the hovered cell at every pointer move across each section, however narrow the cell', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    const overlay = new FakeCanvas();
    viewer.useBackend(new CanvasBackend(new FakeCanvas(), overlay, new FakeCanvas()));
    await viewer.load(longProteinOnGene());
    viewer.resize(600);
    await viewer.settle();
    frames.flush();
    const layout = viewer.frame!.layout;
    for (const s of layout.sections) {
      const box = s.tracks[0].box;
      const y = box.y + box.height / 2;
      let hovered = 0;
      let outlined = 0;
      for (let x = box.x + 0.5; x < box.x + box.width; x += 1.37) {
        overlay.context.take();
        viewer.input.move(x, y);
        frames.flush();
        if (!viewer.hover.active) continue;
        hovered++;
        // the hovered cell's own outline: one on this section's band (the other section's is the projected hover)
        const own = overlay.context.calls('strokeRect').filter((c) => {
          const top = Number(c.split(/[(,]/)[2]);
          return top >= s.y - 4 && top < s.y + s.height;
        });
        if (own.length > 0) outlined++;
      }
      expect(hovered).toBeGreaterThan(300);
      expect(outlined).toBe(hovered);
    }
  });
});

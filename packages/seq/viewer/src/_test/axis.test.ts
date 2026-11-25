/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ViewModel on the axis (_spec/README.md "ViewModel", Follow): a follower that reads along its leader (`follow`,
// the default) moves the same way on screen; on a coding axis (pieces `coding`, no pad) exactly three units a residue;
// a follower on a reversed collapsed axis moved by its leader and leading it back.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import { type FixtureOptions, GeneSection, ProteinSection, proteinOnGene } from './fixture.js';

function placed(options: FixtureOptions) {
  const view = new ViewModel(new SpecModel(proteinOnGene(options)), viewerConfig());
  view.setWidth(ProteinSection, 20);
  view.setWidth(GeneSection, 20);
  const range = (section: string) => {
    const v = view.get(section)!;
    return [v.start, v.end];
  };
  return { view, range };
}

describe('a follower reading along its leader', () => {
  it('on the minus strand its axis is reversed, and a pan moves both the same way on screen', () => {
    const { view, range } = placed({ strand: '-' });
    expect(view.resolved.value?.sections[1].axis.dir).toBe(-1);
    view.setRange(ProteinSection, 1, 3);
    const g = range(GeneSection);
    // the protein's content to the left (a pan right): the gene's too
    view.pan(ProteinSection, 0.5);
    expect(range(GeneSection)[0]).toBeGreaterThan(g[0]);
    const p = range(ProteinSection);
    view.pan(GeneSection, 1);
    expect(range(ProteinSection)[0]).toBeGreaterThan(p[0]);
  });

  it('on a coding axis, exactly three units a residue, both strands, either way', () => {
    for (const strand of ['+', '-'] as const) {
      const { view, range } = placed({ strand, pieces: 'coding' });
      expect(view.resolved.value?.sections[1].axis.ranges).toEqual(strand === '+' ? [2, 7, 12, 19] : [12, 19, 2, 7]);
      const rng = new PCG(strand === '+' ? 3 : 5);
      for (let i = 0; i < 40; i++) {
        const a = rng.float() * 3;
        const b = a + 0.7 + rng.float() * (4 - 0.7 - a);
        if (b > 4) continue;
        view.setRange(ProteinSection, a, b);
        const [s, e] = range(GeneSection);
        expect(s).toBeCloseTo(3 * a, 9);
        expect(e).toBeCloseTo(3 * b, 9);
      }
      // lead then follow returns the window
      view.setRange(GeneSection, 3, 9);
      expect(range(ProteinSection)[0]).toBeCloseTo(1, 9);
      expect(range(ProteinSection)[1]).toBeCloseTo(3, 9);
      expect(range(GeneSection)[0]).toBeCloseTo(3, 9);
      expect(range(GeneSection)[1]).toBeCloseTo(9, 9);
      // a pan of the gene by three units is a residue of the protein, the same way
      view.pan(GeneSection, 1.5);
      expect(range(ProteinSection)[0]).toBeCloseTo(1.5, 9);
      expect(range(GeneSection)[0]).toBeCloseTo(4.5, 9);
    }
  });

  it('a genomic follower on the minus strand moves the other way (spelled by its writer)', () => {
    const { view, range } = placed({ strand: '-', direction: 'genomic' });
    view.setRange(ProteinSection, 1, 3);
    const g = range(GeneSection);
    view.pan(ProteinSection, 0.5);
    expect(range(GeneSection)[0]).toBeLessThan(g[0]);
  });

  it('a protein following its minus-strand gene reads along it (reversed): moved by it, and leading it back', () => {
    const { view, range } = placed({ strand: '-', reverse: true });
    // the protein follows the gene (genomic, whole): its own axis reversed
    expect(view.resolved.value?.sections[1].axis.dir).toBe(-1);
    view.setRange(GeneSection, 6, 18);
    const p = range(ProteinSection);
    view.pan(ProteinSection, 0.25);
    const g = range(GeneSection);
    expect(g[0]).toBeGreaterThan(6);
    // back where it was: the leader returns
    view.pan(ProteinSection, -0.25);
    expect(range(ProteinSection)[0]).toBeCloseTo(p[0], 6);
    expect(range(GeneSection)[0]).toBeCloseTo(6, 6);
  });
});

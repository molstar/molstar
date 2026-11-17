/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import {
  type FixtureOptions,
  GeneSection,
  PlasmidSection,
  ProteinSection,
  proteinOnGene,
  proteinOnPlasmid,
} from './fixture.js';

function setup(options: FixtureOptions = {}, config = viewerConfig()) {
  const spec = new SpecModel(proteinOnGene(options));
  const view = new ViewModel(spec, config);
  const range = (section: string) => {
    const v = view.get(section);
    return v && [v.start, v.end];
  };
  return { spec, view, range };
}

describe('ViewModel: following both ways and over a cyclic origin', () => {
  it('the protein follows the gene at scale 1/3, continuously within a codon, on both strands', () => {
    for (const strand of ['+', '-'] as const) {
      const { view, range } = setup({ strand, reverse: true });
      view.setWidth(GeneSection, 20);
      view.setWidth(ProteinSection, 20);
      // the untranslated ends map nowhere: the hull of what does
      expect(range(ProteinSection)).toEqual([0, 4]);
      // nucleotides 5, 6, 12: residue 1 (plus) or 2 (minus)
      view.setRange(GeneSection, 5, 13);
      expect(range(ProteinSection)).toEqual(strand === '+' ? [1, 2] : [2, 3]);
      view.setRange(GeneSection, 6.5, 14);
      const [s, e] = range(ProteinSection)!;
      if (strand === '+') {
        expect(s).toBeCloseTo(1.5, 10);
        expect(e).toBeCloseTo(7 / 3, 10);
      } else {
        expect(s).toBeCloseTo(5 / 3, 10);
        expect(e).toBeCloseTo(2.5, 10);
      }
    }
  });

  it('a plasmid following a protein whose coding sequence crosses the origin keeps the wrap as two pieces', () => {
    const view = new ViewModel(new SpecModel(proteinOnPlasmid()));
    view.setWidth(ProteinSection, 20);
    view.setWidth(PlasmidSection, 20);
    expect(view.get(PlasmidSection)).toMatchObject({ start: 24, end: 36, pieces: [24, 30, 0, 6] });
    view.setRange(ProteinSection, 0, 1);
    expect(view.get(PlasmidSection)).toMatchObject({ start: 24, end: 27 });
    expect(view.get(PlasmidSection)?.pieces).toBeUndefined();
    view.setRange(ProteinSection, 1.5, 3);
    expect(view.get(PlasmidSection)).toMatchObject({ start: 28.5, end: 33, pieces: [28.5, 30, 0, 3] });

    const hull = new ViewModel(new SpecModel(proteinOnPlasmid()), viewerConfig({ follow: { cyclicArc: 'hull' } }));
    hull.setWidth(ProteinSection, 20);
    hull.setWidth(PlasmidSection, 20);
    expect(hull.get(PlasmidSection)).toMatchObject({ start: 0, end: 30 });
  });

  it('a protein following a cyclic plasmid: a leader range over the coding origin gives the hull', () => {
    const view = new ViewModel(new SpecModel(proteinOnPlasmid({ reverse: true })));
    view.setWidth(PlasmidSection, 20);
    view.setWidth(ProteinSection, 20);
    expect([view.get(ProteinSection)?.start, view.get(ProteinSection)?.end]).toEqual([0, 4]);
    // nucleotides 24–27: residue 0 and a third of residue 1
    view.setRange(PlasmidSection, 20, 28);
    const v = view.get(ProteinSection)!;
    expect(v.start).toBe(0);
    expect(v.end).toBeCloseTo(4 / 3, 10);
  });

  it('a follower is clamped by the most pixels a position may take', () => {
    const { view, range } = setup();
    view.setWidth(ProteinSection, 20);
    view.setWidth(GeneSection, 400);
    view.setRange(ProteinSection, 1, 2);
    // [5, 13) is 8 positions; 400 px allow 12.5 at least
    expect(range(GeneSection)).toEqual([2.75, 15.25]);
  });

  it('a listener of the leader reads its follower already moved (one batch)', () => {
    const { view, range } = setup();
    view.setWidth(ProteinSection, 20);
    view.setWidth(GeneSection, 20);
    const seen: (number[] | undefined)[] = [];
    const off = view.viewport(ProteinSection).listen(() => seen.push(range(GeneSection)));
    view.setRange(ProteinSection, 1, 3);
    expect(seen).toEqual([[5, 16]]);
    off();
  });
});

describe('ViewModel: updates, intents and views', () => {
  it('without grouping, a changed axis still keeps its center and zoom (EC-111)', () => {
    const { spec, view, range } = setup({}, viewerConfig({ viewport: { groupEqualAxes: false } }));
    view.setWidth(ProteinSection, 20);
    view.setRange(ProteinSection, 1, 2);
    spec.apply(proteinOnGene({ range: [1, 4] }));
    expect(range(ProteinSection)).toEqual([0, 1]);
  });

  it('an intent made while the length is pending is carried when the length arrives', () => {
    const pending = (cyclic: boolean) => {
      const b = createSeqSpec();
      const src = b.data.source('src', { format: 'seq', uri: 'g.seq' });
      const g = b.data.sequence('g', { alphabet: 'dna', source: src, cyclic });
      b.views.view('main', (v) => {
        v.section('a', { of: g });
        if (cyclic) v.section('b', { of: g, range: [8, 2] });
      });
      return b.build();
    };
    const known = () => {
      const b = createSeqSpec();
      const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(6) });
      b.views.view('main', (v) => void v.section('a', { of: g }));
      return b.build();
    };
    const spec = new SpecModel(pending(false));
    const view = new ViewModel(spec);
    view.setWidth('/views/main/a', 20);
    view.setRange('/views/main/a', 3, 9);
    expect(view.get('/views/main/a')).toMatchObject({ start: 3, end: 9, pending: true });
    spec.apply(known());
    expect(view.get('/views/main/a')).toMatchObject({ start: 3, end: 9, pending: false });
    // two pending sections of one system group only when their own ranges are equal
    const cyc = new ViewModel(new SpecModel(pending(true)));
    expect(cyc.groupOf('/views/main/a')).toEqual(['/views/main/a']);
  });

  it('a detached follower stays detached across an update of another part of the spec', () => {
    const { spec, view, range } = setup({}, viewerConfig({ follow: { followerMoves: 'detach' } }));
    view.setWidth(ProteinSection, 20);
    view.setWidth(GeneSection, 20);
    view.setRange(ProteinSection, 1, 3);
    view.pan(GeneSection, 2);
    spec.apply(proteinOnGene({ title: 'retitled' }));
    expect(view.isFollowing(GeneSection)).toBe(false);
    expect(range(GeneSection)).toEqual([7, 18]);
  });

  it('a move that is no move does not detach a follower', () => {
    const { view } = setup({}, viewerConfig({ follow: { followerMoves: 'detach' } }));
    expect(view.zoom(GeneSection, Number.NaN)).toBe(false);
    expect(view.isFollowing(GeneSection)).toBe(true);
  });

  it('show switches views: the other view’s sections get viewports, the first one’s go', () => {
    const { view } = setup({ other: true });
    view.show('/views/other');
    expect(view.get('/views/other/gene')).toMatchObject({ start: 0, end: 24 });
    expect(view.get(ProteinSection)).toBeUndefined();
    view.show(undefined);
    expect(view.get(ProteinSection)).toMatchObject({ start: 0, end: 4 });
  });
});

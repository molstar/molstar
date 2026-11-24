/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import {
  type FixtureOptions,
  Gene2Section,
  GeneSection,
  ProteinSection,
  proteinOnGene,
  TwinSection,
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

describe('ViewModel: scroll groups and follow (EC-113, EC-114)', () => {
  it('sections on equal axes share one viewport; without grouping each has its own', () => {
    const { view, range } = setup({ twin: true });
    view.setWidth(ProteinSection, 20);
    view.setWidth(TwinSection, 20);
    expect(view.groupOf(ProteinSection)).toEqual([ProteinSection, TwinSection]);
    view.setRange(ProteinSection, 1, 3);
    expect(range(TwinSection)).toEqual([1, 3]);
    view.pan(TwinSection, 0.5);
    expect(range(ProteinSection)).toEqual([1.5, 3.5]);

    const solo = setup({ twin: true }, viewerConfig({ viewport: { groupEqualAxes: false } }));
    solo.view.setWidth(ProteinSection, 20);
    solo.view.setWidth(TwinSection, 20);
    solo.view.setRange(ProteinSection, 1, 3);
    expect(solo.range(TwinSection)).toEqual([0, 4]);
  });

  it('the gene follows the protein through the codon map at scale 3, over the intron, on both strands', () => {
    for (const strand of ['+', '-'] as const) {
      const { view, range } = setup({ strand, direction: 'genomic' });
      view.setWidth(ProteinSection, 20);
      view.setWidth(GeneSection, 20);
      expect(view.get(GeneSection)?.leader).toBe(ProteinSection);
      // fitted: the whole coding hull, from the first exon's start to the last's end
      expect(range(GeneSection)).toEqual([2, 19]);
      // residues 1 and 2 are 5, 6, 12 and 13–15 (plus) or 13–15 and 12, 6, 5 (minus): the same hull
      view.setRange(ProteinSection, 1, 3);
      expect(range(GeneSection)).toEqual([5, 16]);
      // within a run it moves smoothly: half a residue is a codon's half, at the strand's end
      view.setRange(ProteinSection, 1.5, 3);
      expect(range(GeneSection)).toEqual(strand === '+' ? [6.5, 16] : [5, 14.5]);
    }
  });

  it('a move on a follower is ignored, or detaches it, when configured so', () => {
    const { view, range } = setup({}, viewerConfig({ follow: { followerMoves: 'ignore' } }));
    view.setWidth(ProteinSection, 20);
    view.setWidth(GeneSection, 20);
    view.setRange(ProteinSection, 1, 3);
    expect(view.pan(GeneSection, 2)).toBe(false);
    expect(range(ProteinSection)).toEqual([1, 3]);
    expect(range(GeneSection)).toEqual([5, 16]);

    const d = setup({}, viewerConfig({ follow: { followerMoves: 'detach' } }));
    d.view.setWidth(ProteinSection, 20);
    d.view.setWidth(GeneSection, 20);
    d.view.setRange(ProteinSection, 1, 3);
    expect(d.view.pan(GeneSection, 2)).toBe(true);
    expect(d.range(GeneSection)).toEqual([7, 18]);
    expect(d.range(ProteinSection)).toEqual([1, 3]);
    expect(d.view.isFollowing(GeneSection)).toBe(false);
    d.view.setRange(ProteinSection, 0, 2);
    expect(d.range(GeneSection)).toEqual([7, 18]);
    // a fit attaches it again
    d.view.fit(GeneSection);
    expect(d.view.isFollowing(GeneSection)).toBe(true);
    expect(d.range(GeneSection)).toEqual([2, 13]);
  });

  it('follows chain: a follower of a follower moves with the first leader', () => {
    const { view, range } = setup({ chain: true });
    for (const s of [ProteinSection, GeneSection, Gene2Section]) view.setWidth(s, 20);
    view.setRange(ProteinSection, 1, 3);
    expect(range(Gene2Section)).toEqual([5, 16]);
    expect(view.get(Gene2Section)?.leader).toBe(GeneSection);
  });
});

describe('ViewModel: viewports across updates and widths (EC-111, EC-112)', () => {
  it('a move before the section has a width is kept as an intent, then placed', () => {
    const { view, range } = setup();
    view.setRange(ProteinSection, 1, 2);
    expect(view.get(ProteinSection)).toMatchObject({ start: 1, end: 2, pending: true });
    view.zoom(ProteinSection, 2);
    expect(range(ProteinSection)).toEqual([1.25, 1.75]);
    // placed at 20 px: at most 32 px a position, so at least 0.625 positions, about the center
    view.setWidth(ProteinSection, 20);
    expect(view.get(ProteinSection)).toMatchObject({ start: 1.1875, end: 1.8125, pending: false });
  });

  it('a section hidden for lack of width keeps its zoom; a resize keeps it about the center', () => {
    const { view, range } = setup();
    view.setWidth(ProteinSection, 20);
    view.setRange(ProteinSection, 1, 2);
    view.setWidth(ProteinSection, 0);
    expect(view.get(ProteinSection)).toMatchObject({ start: 1, end: 2, pending: true });
    view.setWidth(ProteinSection, 20);
    expect(range(ProteinSection)).toEqual([1, 2]);
    view.setWidth(ProteinSection, 40);
    expect(range(ProteinSection)).toEqual([0.5, 2.5]);
    // a fitted viewport stays fitted
    view.fit(ProteinSection);
    view.setWidth(ProteinSection, 10);
    expect(range(ProteinSection)).toEqual([0, 4]);

    const kept = setup({}, viewerConfig({ viewport: { onResize: 'keep-range' } }));
    kept.view.setWidth(ProteinSection, 20);
    kept.view.setRange(ProteinSection, 1, 3);
    kept.view.setWidth(ProteinSection, 40);
    expect(kept.range(ProteinSection)).toEqual([1, 3]);
  });

  it('an equal axis keeps its viewport across a patch; a changed one keeps center and zoom when shown, else fits', () => {
    const { spec, view, range } = setup();
    view.setWidth(ProteinSection, 20);
    view.setRange(ProteinSection, 1, 2);
    spec.apply(proteinOnGene({ title: 'retitled' }));
    expect(range(ProteinSection)).toEqual([1, 2]);
    // the axis now shows [1, 4): the center 1.5 is offset 0.5
    spec.apply(proteinOnGene({ range: [1, 4] }));
    expect(range(ProteinSection)).toEqual([0, 1]);
    // centered on position 0.4, which [1, 4) does not show: fits
    spec.apply(proteinOnGene());
    view.setRange(ProteinSection, 0, 0.8);
    spec.apply(proteinOnGene({ range: [1, 4] }));
    expect(range(ProteinSection)).toEqual([0, 3]);
    // the follower moves with its leader's new viewport
    view.setWidth(GeneSection, 20);
    expect(range(GeneSection)).toEqual([5, 19]);
  });

  it('a section’s viewport is told only to its readers, and not for an equal one', () => {
    const { view } = setup({ twin: true });
    view.setWidth(ProteinSection, 20);
    let told = 0;
    const off = view.viewport(GeneSection).listen(() => told++);
    view.setTop(ProteinSection, 10);
    expect(told).toBe(0);
    view.setRange(ProteinSection, 1, 3);
    expect(told).toBe(1);
    view.setRange(ProteinSection, 1, 3);
    expect(told).toBe(1);
    off();
  });
});

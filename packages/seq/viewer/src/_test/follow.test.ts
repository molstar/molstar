/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { viewerConfig } from '../config.js';
import { windowBehavior } from '../behaviors/window.js';
import { SpecModel } from '../models/spec.js';
import { followerRange, leaderRange } from '../models/view-follow.js';
import { ViewModel } from '../models/view.js';
import {
  type FixtureOptions,
  Gene2Section,
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
      const { view, range } = setup({ strand, reverse: true, direction: 'genomic' });
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

describe('the follow functions write into the caller’s record', () => {
  it('a range is the record asked for, so two views never share one', () => {
    const { view } = setup();
    view.setWidth(GeneSection, 20);
    view.setWidth(ProteinSection, 20);
    const section = view.resolved.value!.sections.find((s) => s.follow)!;
    const gene = view.get(section.follow!.leader)!;
    const a = { start: 0, end: 0 };
    const b = { start: 0, end: 0 };
    expect(leaderRange(section, gene, 5, 13, 1, 1, undefined, 'shortest', a)).toBe(a);
    expect(leaderRange(section, gene, 8, 14, 1, 1, undefined, 'shortest', b)).toBe(b);
    expect(a).not.toEqual(b);
    const p = { start: 0, end: 0, pieces: undefined };
    const q = { start: 0, end: 0, pieces: undefined };
    expect(followerRange(section, gene, undefined, 20, 'shortest', 40, p)).toBe(p);
    expect(followerRange(section, { start: 1, end: 2 }, undefined, 20, 'shortest', 40, q)).toBe(q);
    expect(p).not.toEqual(q);
  });
});

describe('a follower’s window dragged by its body (a pan)', () => {
  it('2x the axis past either end: the window keeps its width and stops at the end', () => {
    const { view, range } = setup();
    view.setWidth(GeneSection, 20);
    view.setWidth(ProteinSection, 20);
    const w = windowBehavior(view, viewerConfig());
    view.setRange(ProteinSection, 1, 3);
    const at = range(ProteinSection)!;
    const width = at[1] - at[0];
    w.begin(ProteinSection, 'body', 2);
    w.dragTo(2 + 200);
    const hi = range(ProteinSection)!;
    expect(hi[1] - hi[0]).toBeCloseTo(width, 6);
    w.dragTo(2 - 200);
    const lo = range(ProteinSection)!;
    expect(lo[1] - lo[0]).toBeCloseTo(width, 6);
    expect(lo[0]).toBeGreaterThanOrEqual(0);
    w.end();
  });
});

describe('ViewModel: a follower leads (follow.followerMoves: lead)', () => {
  const close = (r: number[] | undefined, expected: number[]) => {
    expect(r).toBeDefined();
    expect(r![0]).toBeCloseTo(expected[0], 9);
    expect(r![1]).toBeCloseTo(expected[1], 9);
  };
  const placed = (options: FixtureOptions = {}) => {
    const s = setup(options);
    s.view.setWidth(ProteinSection, 20);
    s.view.setWidth(GeneSection, 20);
    return s;
  };

  it('a pan on the gene moves the protein a third as far, and the gene as asked', () => {
    const { view, range } = placed();
    view.setRange(ProteinSection, 1, 3);
    expect(view.canMove(GeneSection)).toBe(true);
    expect(view.pan(GeneSection, 0.5)).toBe(true);
    close(range(GeneSection), [5.5, 16.5]);
    close(range(ProteinSection), [7 / 6, 19 / 6]);
    expect(view.isFollowing(GeneSection)).toBe(true);
  });

  it('on the minus strand the protein moves the other way, and the gene still as asked', () => {
    const { view, range } = placed({ strand: '-', direction: 'genomic' });
    view.setRange(ProteinSection, 1.5, 3);
    close(range(GeneSection), [5, 14.5]);
    // the gene's content to the left: the protein's to the right
    view.pan(GeneSection, -0.5);
    close(range(GeneSection), [4.5, 14]);
    close(range(ProteinSection), [5 / 3, 19 / 6]);
  });

  it('a zoom on the gene about a point keeps it in place, on both strands', () => {
    for (const strand of ['+', '-'] as const) {
      const { view, range } = placed({ strand, direction: 'genomic' });
      view.setRange(ProteinSection, 0, 4);
      close(range(GeneSection), [2, 19]);
      expect(view.zoom(GeneSection, 1.25, 17)).toBe(true);
      close(range(GeneSection), [5, 18.6]);
      expect(range(ProteinSection)![1] - range(ProteinSection)![0]).toBeLessThan(4);
    }
  });

  it('a codon split by the intron: each of its nucleotides is an edge the follower can show', () => {
    // minus strand: residue 2 is 12, 6 and 5, place by place, so a start moved onto 5 or 6 stays where it was asked
    const { view, range } = placed({ strand: '-', direction: 'genomic' });
    view.setRange(ProteinSection, 1.5, 3);
    view.pan(GeneSection, 0.5);
    expect(range(GeneSection)![0]).toBeCloseTo(5.5, 6);
    expect(range(ProteinSection)![1]).toBeCloseTo(8.5 / 3, 6);
    view.pan(GeneSection, -0.5);
    expect(range(GeneSection)![0]).toBeCloseTo(5, 6);
  });

  it('a pan over an intron moves into it as far as the leader can show, then steps over it, both ways', () => {
    const { view, range } = placed();
    view.setRange(ProteinSection, 1, 3);
    // [5, 16] by 2: the start in the intron [7, 12), no place's image: it steps over to the next place, 12 (residue
    // 1's third nucleotide), the end where it was asked
    view.pan(GeneSection, 2);
    close(range(GeneSection), [12, 18]);
    // on by 2: the start would be in the intron past that image, so it steps over to residue 2's first nucleotide
    view.pan(GeneSection, 2);
    const [s, e] = range(GeneSection)!;
    expect(s).toBeGreaterThanOrEqual(12);
    expect(range(ProteinSection)![0]).toBeGreaterThanOrEqual(2);
    expect(e).toBeGreaterThan(18);
    // back again: it moves left, never sticks
    view.pan(GeneSection, -2);
    expect(range(GeneSection)![0]).toBeLessThan(s);
    view.pan(GeneSection, -4);
    expect(range(GeneSection)![0]).toBeLessThan(7);
  });

  it('small steps across the intron never stick, either way', () => {
    const { view, range } = placed();
    view.setRange(ProteinSection, 1, 2);
    let last = range(GeneSection)![0];
    for (let i = 0; i < 30; i++) {
      view.pan(GeneSection, 0.4);
      const now = range(GeneSection)![0];
      // at the protein's end the pan stops: the last step taken is where the way back starts
      if (range(ProteinSection)![1] >= 4) {
        last = Math.max(last, now);
        break;
      }
      expect(now).toBeGreaterThan(last);
      last = now;
    }
    expect(last).toBeGreaterThan(12);
    for (let i = 0; i < 30 && range(ProteinSection)![0] > 0; i++) {
      view.pan(GeneSection, -0.4);
      const now = range(GeneSection)![0];
      expect(now).toBeLessThan(last);
      last = now;
    }
    expect(last).toBeLessThan(7);
  });

  it('a protein following its gene leads too: a pan on it moves the gene three times as far', () => {
    const { view, range } = placed({ reverse: true });
    view.setRange(GeneSection, 5, 13);
    close(range(ProteinSection), [1, 2]);
    view.pan(ProteinSection, 0.5);
    close(range(ProteinSection), [1.5, 2.5]);
    close(range(GeneSection), [6.5, 14.5]);
  });

  it('a chain leads up to its first leader', () => {
    const { view, range } = setup({ chain: true });
    for (const s of [ProteinSection, GeneSection, Gene2Section]) view.setWidth(s, 20);
    view.setRange(ProteinSection, 1, 3);
    view.pan(Gene2Section, 0.5);
    close(range(Gene2Section), [5.5, 16.5]);
    close(range(ProteinSection), [7 / 6, 19 / 6]);
  });

  it('a move with invalid arguments moves nothing; a pan that is no pan is one when the leader can move', () => {
    const { view, range } = placed();
    view.setRange(ProteinSection, 1, 3);
    expect(view.zoom(GeneSection, Number.NaN)).toBe(false);
    expect(view.pan(GeneSection, Number.POSITIVE_INFINITY)).toBe(false);
    expect(view.pan(GeneSection, 0)).toBe(true);
    expect(range(ProteinSection)).toEqual([1, 3]);
  });

  it('a pan on a follower its own clamp widened moves the leader, which keeps its width (protein follows gene)', () => {
    const { view, range } = placed({ reverse: true });
    view.setRange(GeneSection, 6, 6.1);
    // 20 px at most 32 px a position: the gene shows 0.625 nucleotides, the protein a third of that, widened
    const [g0, g1] = range(GeneSection)!;
    expect(g1 - g0).toBeCloseTo(0.625, 9);
    view.pan(ProteinSection, 0.01);
    close(range(GeneSection), [g0 + 0.03, g1 + 0.03]);
  });

  it('a pan on a follower its own clamp widened moves the leader, which keeps its width (gene follows protein)', () => {
    const { view, range } = setup();
    view.setWidth(ProteinSection, 20);
    view.setWidth(GeneSection, 400);
    view.setRange(ProteinSection, 1, 2);
    // [5, 13) widened to 12.5 nucleotides by the gene's own clamp
    expect(range(GeneSection)).toEqual([2.75, 15.25]);
    view.pan(GeneSection, 0.3);
    close(range(ProteinSection), [1.1, 2.1]);
  });

  it('a zoom out on a follower past its axis reaches the leader’s whole range, never stalls', () => {
    const { view, range } = placed();
    view.setRange(ProteinSection, 1, 3);
    for (let i = 0; i < 6; i++) view.zoom(GeneSection, 0.5);
    expect(range(ProteinSection)).toEqual([0, 4]);
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

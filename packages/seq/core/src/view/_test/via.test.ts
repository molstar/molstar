/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's `via` in the resolved view (_spec/README.md "Resolved view", Via): two sections of one chain, each naming
// one of two maps joining the sequence and the chain's molecule, show it each its own way (a figure comparing SIFTS
// with an alignment); a section naming none is ambiguous (EC-24). Its follow, its `direction: { along }`, its layers'
// data and the links onto it take its `via`.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { projectRanges, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { resolveView } from '../view.js';

/** `up` (20) onto `ent` (16) by SIFTS (from 2) and by an alignment (from 3); unit A of ent; sections over them. */
function comparison(): Node {
  const b = createSeqSpec({ title: 'SIFTS against an alignment' });
  const up = b.data.sequence('up', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFS' });
  const ent = b.data.sequence('ent', { alphabet: 'protein', letters: 'EPQSDPSVEPPLSQET' });
  let a = '';
  b.data.complex('c', {}, (c) => {
    a = c.unit('A', { of: ent, chain: 'A' });
  });
  const sifts = b.data.map('sifts', {
    from: up,
    to: ent,
    method: 'structure',
    segments: { from: [2], to: [0], length: [16], provenance: ['sifts'] },
  });
  const aln = b.data.map('aln', {
    from: up,
    to: ent,
    method: 'segments',
    segments: { from: [3], to: [0], length: [15] },
  });
  const site = b.data.ranges('site', { target: up, rows: { start: [10], end: [12] } });
  b.views.view('v', (v) => {
    const lead = v.section('lead', { of: up });
    v.section('by-sifts', { of: a, follow: lead, via: sifts }, (s) =>
      s.track('t').layer('l', { type: 'ranges', data: site }),
    );
    v.section('by-aln', { of: a, follow: lead, via: [aln], direction: { along: sifts } }, (s) =>
      s.track('t').layer('l', { type: 'ranges', data: site }),
    );
    v.section('neither', { of: a, follow: lead });
  });
  return b.build();
}

describe('a section’s via (EC-24)', () => {
  const model = resolveSpec(comparison());
  const view = resolveView(model, '/views/v')!;
  const section = (id: string) => view.sections.find((s) => s.address === `/views/v/${id}`)!;

  it('follows its leader through the map it names; a section naming none is ambiguous', () => {
    expect(model.issues).toEqual([]);
    expect(section('by-sifts').via).toEqual(['/data/sifts']);
    expect(Array.from(projectRanges(section('by-sifts').follow!.path, [10, 11]))).toEqual([8, 9]);
    expect(Array.from(projectRanges(section('by-aln').follow!.path, [10, 11]))).toEqual([7, 8]);
    expect(section('neither').follow).toBeUndefined();
    expect(section('neither').via).toBeUndefined();
    expect(view.issues.map((i) => [i.code, i.address])).toEqual([['seq.view.follow.unplaced', '/views/v/neither']]);
  });

  it('places its layers’ data through its map', () => {
    const placed = (id: string) => {
      const l = section(id).tracks[0].layers[0];
      return l.placed && l.path ? Array.from(projectRanges(l.path, [10, 12])) : undefined;
    };
    expect(placed('by-sifts')).toEqual([8, 10]);
    expect(placed('by-aln')).toEqual([7, 9]);
  });

  it('what is projected onto it or from it takes its via: one hover lights different residues in the two sections', () => {
    const link = (from: string, to: string) =>
      view.links.find((l) => l.from === `/views/v/${from}` && l.to === `/views/v/${to}`);
    expect(Array.from(projectRanges(link('lead', 'by-sifts')!.path, [10, 11]))).toEqual([8, 9]);
    expect(Array.from(projectRanges(link('lead', 'by-aln')!.path, [10, 11]))).toEqual([7, 8]);
    // from it too (the author, 2026-11-20): onto the lead (no via of its own) from the chain, through its map
    expect(Array.from(projectRanges(link('by-sifts', 'lead')!.path, [8, 9]))).toEqual([10, 11]);
    expect(Array.from(projectRanges(link('by-aln', 'lead')!.path, [7, 8]))).toEqual([10, 11]);
    expect(link('lead', 'neither')).toBeUndefined();
    // between sections of one system: the identity
    expect(link('by-sifts', 'by-aln')?.path.edges).toEqual([]);
  });

  it('warns of a via naming a map that joins no pair on any of its paths, naming the map', () => {
    const b = createSeqSpec();
    const up = b.data.sequence('up', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    const ent = b.data.sequence('ent', { alphabet: 'protein', letters: 'EPQSDPSV' });
    const other = b.data.sequence('other', { alphabet: 'protein', letters: 'MEEP' });
    const m = b.data.map('m', { from: up, to: ent, method: 'offset', by: -2 });
    b.views.view('v', (v) => {
      const lead = v.section('lead', { of: up });
      v.section('used', { of: ent, follow: lead, via: m });
      v.section('elsewhere', { of: other, via: m });
    });
    const v = resolveView(resolveSpec(b.build()), '/views/v')!;
    const unused = v.issues.filter((i) => i.code === 'seq.resolve.via.unused');
    expect(unused.map((i) => [i.address, i.severity, i.places?.map((p) => p.address)])).toEqual([
      ['/views/v/elsewhere', 'warning', ['/data/m']],
    ]);
    expect(unused[0].message).toContain('/data/m');
  });

  it('reads `direction: { along }` through its via where the map does not reach its system itself', () => {
    expect(section('by-aln').axis.dir).toBe(1);
    expect(section('by-aln').axis.along).toBe('/data/sifts');
    expect(view.issues.some((i) => i.code === 'seq.view.direction.unresolved')).toBe(false);
  });
});

describe('two sections naming rival maps of one pair', () => {
  it('the link between them takes the target’s map and says so (seq.view.via.conflict); it never vanishes', () => {
    const b = createSeqSpec();
    const up = b.data.sequence('up', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFS' });
    const ent = b.data.sequence('ent', { alphabet: 'protein', letters: 'EPQSDPSVEPPLSQET' });
    const sifts = b.data.map('sifts', { from: up, to: ent, method: 'offset', by: -2 });
    const aln = b.data.map('aln', { from: up, to: ent, method: 'offset', by: -3 });
    const site = b.data.ranges('site', { target: up, rows: { start: [10], end: [12] } });
    b.views.view('v', (v) => {
      v.section('up', { of: up, via: sifts }, (s) => s.track('t').layer('l', { type: 'ranges', data: site }));
      v.section('ent', { of: ent, via: aln }, (s) => s.track('t').layer('l', { type: 'ranges', data: site }));
    });
    const view = resolveView(resolveSpec(b.build()), '/views/v')!;
    const link = (from: string, to: string) =>
      view.links.find((l) => l.from === `/views/v/${from}` && l.to === `/views/v/${to}`);
    // onto `ent`: its map (aln, by -3); onto `up`: its map (sifts, from ent by +2)
    expect(Array.from(projectRanges(link('up', 'ent')!.path, [10, 11]))).toEqual([7, 8]);
    expect(Array.from(projectRanges(link('ent', 'up')!.path, [8, 9]))).toEqual([10, 11]);
    expect(view.issues.filter((i) => i.code === 'seq.view.via.conflict').map((i) => i.address)).toEqual([
      '/views/v/ent',
      '/views/v/up',
    ]);
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view of the p53 example and of small specs: sections and axes, tracks, layers with their data, paths,
// themes and fallbacks, follow links with their scales, and what does not resolve.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { createResolver, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { EXONS, P53_LETTERS, p53Spec } from '../../lang/_test/p53.js';
import { p53GeometrySpec } from '../../resolve/_test/p53-geometry.js';
import { followRange, resolveView, resolveViews } from '../view.js';

const model = resolveSpec(p53Spec());
const main = resolveView(model, '/views/main')!;
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

describe('resolved view: the p53 example', () => {
  it('has its sections, each an axis on a system', () => {
    expect(main.address).toBe('/views/main');
    expect(main.layout).toBe('tracks');
    expect(main.issues).toEqual([]);
    expect(main.sections.map((s) => [s.address, s.axis])).toEqual([
      [
        '/views/main/protein',
        { system: '/data/p53', length: 30, cyclic: false, ranges: [0, 30], dir: 1, pieces: 'whole' },
      ],
      // the gene's length is not known yet: its axis is pending
      ['/views/main/gene', { system: '/data/tp53', cyclic: false, dir: -1, pieces: 'whole' }],
    ]);
    expect(resolveViews(model)).toEqual([main]);
    expect(resolveView(model, '/views/main')).toBe(main);
    expect(resolveView(model, '/views/none')).toBeUndefined();
  });

  it('has its tracks and layers, each layer’s data on its section’s system', () => {
    const [letters, domains] = main.sections[0].tracks;
    expect([letters.address, letters.label, domains.label]).toEqual([
      '/views/main/protein/letters',
      'letters',
      'Domains',
    ]);
    expect(letters.layers).toEqual([
      {
        index: 0,
        address: '/views/main/protein/letters/letters',
        type: 'letters',
        fallbacks: ['letters'],
        params: {},
        theme: {},
        data: { kind: 'address', address: '/data/p53', system: '/data/p53', shape: 'letters' },
        placed: true,
      },
    ]);
    expect(domains.layers[0].data).toEqual({
      kind: 'address',
      address: '/data/domains',
      system: '/data/p53',
      shape: 'ranges',
    });
    expect(domains.layers[0].path).toBeUndefined();
    expect(main.sections[1].tracks).toEqual([]);
  });

  it('resolves the gene’s follow to the codon path, three nucleotides a residue', () => {
    const follow = main.sections[1].follow!;
    expect([follow.leader, follow.scale, follow.dir, follow.path.from, follow.path.to, follow.path.step]).toEqual([
      '/views/main/protein',
      3,
      -1,
      '/data/p53',
      '/data/tp53',
      1,
    ]);
    expect(followRange(main.sections[1], [0, 1])).toEqual([197, 200]);
    expect(followRange(main.sections[1], [0, 30])).toEqual([100, 200]);
    expect(followRange(main.sections[0], [0, 1])).toBeUndefined();
  });

  it('moves a follower smoothly: the continuous image of a leader position, against it on the minus strand', () => {
    const follow = main.sections[1].follow!;
    expect(follow.point(0)).toBe(200);
    expect(follow.point(0.5)).toBe(198.5);
    expect(follow.point(1)).toBe(197);
    expect(follow.point(31)).toBeUndefined();
  });

  it('precomposes the paths between its sections', () => {
    expect(main.links.map((l) => [l.from, l.to, l.path.step])).toEqual([
      ['/views/main/protein', '/views/main/gene', 1],
      ['/views/main/gene', '/views/main/protein', -1],
    ]);
    expect(main.links[0].path).toBe(main.sections[1].follow!.path);
  });
});

describe('resolved view: a hidden track', () => {
  it('is kept in the spec and left out of the view, its row collapsed; its layers keep their own `hidden`', () => {
    const spec = (hidden: boolean) => {
      const b = createSeqSpec();
      const p = b.data.sequence('p', { alphabet: 'protein', letters: P53_LETTERS });
      b.views.view('main', (v) => {
        const s = v.section('s', { of: p });
        s.track('a', { hidden }).layer('l', { type: 'letters', hidden: true });
        s.track('b').layer('l', { type: 'letters' });
      });
      return b.build();
    };
    const shown = resolveView(resolveSpec(spec(false)), '/views/main')!;
    expect(shown.sections[0].tracks.map((t) => t.address)).toEqual(['/views/main/s/a', '/views/main/s/b']);
    expect(shown.sections[0].tracks[0].layers[0].params).toEqual({ hidden: true });
    const hidden = resolveView(resolveSpec(spec(true)), '/views/main')!;
    expect(hidden.sections[0].tracks.map((t) => t.address)).toEqual(['/views/main/s/b']);
  });
});

describe('resolved view: layers on other systems, queries, themes', () => {
  /** p53 and its gene inline; a view with the protein's domains drawn on the gene and on a lone sequence. */
  function spec(): Node {
    const b = createSeqSpec();
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
    const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: 'A'.repeat(200) });
    const lone = b.data.sequence('lone', { alphabet: 'protein', letters: 'MKW' });
    b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: EXONS, strand: '+' });
    const domains = b.data.ranges('domains', {
      target: p53,
      rows: { start: new Int32Array([0]), end: new Int32Array([14]) },
    });
    b.views.view('v', (v) => {
      const g = v.section('gene', { of: gene, range: [90, 200] });
      g.track('t').layer('d', { type: 'ranges', data: domains, color: { kind: 'uniform', value: 'red' } });
      g.track('q').layer('ranges', {
        type: 'ranges',
        data: '/data/domains/(items :overlapping /data/p53/(range :from 1 :to 5))',
      });
      const l = v.section('lone', { of: lone });
      l.track('t').layer('ranges', { type: 'ranges', data: domains });
      v.section('follower', { of: lone, follow: g });
    });
    return b.build();
  }
  const m = resolveSpec(spec());
  const v = resolveView(m, '/views/v')!;

  it('carries a layer’s data along the static path to its section, with its theme as a value', () => {
    const d = v.sections[0].tracks[0].layers[0];
    expect([d.address, d.placed, d.path?.step, d.fallbacks]).toEqual(['/views/v/gene/t/d', true, 1, ['ranges']]);
    expect(d.theme).toEqual({ color: { kind: 'uniform', value: 'red' } });
    expect(d.params).toEqual({ color: { kind: 'uniform', value: 'red' } });
  });

  it('keeps a query as a query, on the system its first node is on', () => {
    const q = v.sections[0].tracks[1].layers[0];
    expect(q.data).toEqual({
      kind: 'query',
      address: '/data/domains/(items :overlapping /data/p53/(range :from 1 :to 5))',
      system: '/data/p53',
      shape: 'ranges',
    });
    expect(q.placed && q.path?.to).toBe('/data/tp53');
  });

  it('reports data that does not reach its section and a follow with no path, and draws neither', () => {
    const lone = v.sections[1].tracks[0].layers[0];
    expect([lone.placed, lone.path]).toEqual([false, undefined]);
    expect(v.sections[2].follow).toBeUndefined();
    // the follows are resolved before the axes (a follower's direction reads its path), so its issue comes first
    expect(codes(v.issues)).toEqual(['seq.view.follow.unplaced', 'seq.view.layer.unplaced']);
    expect(v.issues.map((i) => i.address)).toEqual(['/views/v/follower', '/views/v/lone/t/ranges']);
    expect(v.issues[1].hint).toContain('No path of maps joins the systems');
  });

  it('keeps a section’s range', () => {
    expect(v.sections[0].axis).toEqual({
      system: '/data/tp53',
      length: 200,
      cyclic: false,
      ranges: [90, 200],
      dir: 1,
      pieces: 'whole',
    });
  });

  it('reports a layer whose data names nothing (an unchecked spec)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MK' });
    b.views.view('v', (view) =>
      view.section('s', { of: s }).track('t').layer('ranges', { type: 'ranges', data: '/data/none' }),
    );
    const r = resolveView(resolveSpec(b.build()), '/views/v')!;
    expect(codes(r.issues)).toEqual(['seq.view.layer.data']);
    expect(r.sections[0].tracks[0].layers[0].placed).toBe(false);
  });
});

describe('resolved view: cyclic sections', () => {
  const b = createSeqSpec();
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGTAC' });
  const pl = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
  b.data.map('m', { from: a, to: pl, method: 'offset', by: 7 });
  b.views.view('v', (v) => {
    const lead = v.section('a', { of: a });
    v.section('pl', { of: pl, range: [8, 2], follow: lead });
  });
  const v = resolveView(resolveSpec(b.build()), '/views/v')!;

  it('a range that wraps the origin is two ranges (EC-9, EC-10)', () => {
    // in display order: from its start to the origin, then on
    expect(v.sections[1].axis.ranges).toEqual([8, 10, 0, 2]);
  });

  it('a follower on a cyclic system takes the shortest arc, in the axis’s shape (EC-12)', () => {
    expect(v.sections[1].follow?.scale).toBe(1);
    expect(followRange(v.sections[1], [1, 5])).toEqual([8, 10, 0, 2]);
    expect(followRange(v.sections[1], [0, 2])).toEqual([7, 9]);
    expect(followRange(v.sections[1], [1, 5], { arc: 'hull' })).toEqual([0, 10]);
  });

  it('a reversed cyclic follower takes its arc in its own display order: the piece below the origin first', () => {
    const c = createSeqSpec();
    const a = c.data.sequence('a', { alphabet: 'dna', letters: 'ACGTAC' });
    const pl = c.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
    c.data.map('m', { from: a, to: pl, method: 'offset', by: 7 });
    c.views.view('v', (view) => {
      const lead = view.section('a', { of: a });
      view.section('pl', { of: pl, range: [8, 2], follow: lead, direction: 'reversed' });
    });
    const f = resolveView(resolveSpec(c.build()), '/views/v')!.sections[1];
    expect(f.axis.ranges).toEqual([0, 2, 8, 10]);
    expect(followRange(f, [1, 5])).toEqual([0, 2, 8, 10]);
    expect(followRange(f, [0, 2])).toEqual([7, 9]);
  });

  it('a leader range over its cyclic origin is split before it is carried', () => {
    const c = createSeqSpec();
    const pl = c.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
    const copy = c.data.sequence('copy', { alphabet: 'dna', letters: 'ACGTACGTAC' });
    c.data.map('m', { from: pl, to: copy, method: 'offset', by: 0 });
    c.views.view('v', (view) => {
      const lead = view.section('pl', { of: pl });
      view.section('copy', { of: copy, follow: lead });
    });
    const f = resolveView(resolveSpec(c.build()), '/views/v')!.sections[1];
    // [8, 2) on the plasmid: 8, 9, 0, 1; their hull on the copy, which is not cyclic
    expect(followRange(f, [8, 2])).toEqual([0, 10]);
    expect(followRange(f, [8, 10, 0, 2])).toEqual([0, 10]);
    expect(followRange(f, [3, 5])).toEqual([3, 5]);
  });
});

describe('resolved view: p53 on chr17', () => {
  it('draws the protein’s domains on the gene through the codon path, following the protein', () => {
    const m = resolveSpec(p53GeometrySpec());
    const v = resolveView(m, '/views/main')!;
    expect(v.issues).toEqual([]);
    const gene = v.sections[1];
    expect(gene.axis).toEqual({ system: '/data/chr17', cyclic: false, dir: -1, pieces: 'whole' });
    expect(gene.tracks[0].layers[0].path?.step).toBe(1);
    expect(followRange(gene, [93, 292])).toEqual([7673743, 7676089]);
    expect(gene.follow?.dir).toBe(-1);
  });

  it('an edit through a resolver rebuilds the view and keeps its paths', () => {
    const resolver = createResolver();
    const spec = p53GeometrySpec();
    const a = resolveView(resolver.resolve(spec), '/views/main')!;
    // a track's label changes: the data is the same
    const views = spec.children!.find((c) => c.kind === 'views')!;
    const relabel = (n: Node): Node =>
      n.kind === 'track'
        ? { ...n, params: { ...n.params, label: 'Domains' } }
        : { ...n, children: n.children?.map(relabel) };
    const next: Node = { ...spec, children: spec.children!.map((c) => (c === views ? relabel(c) : c)) };
    const b = resolveView(resolver.resolve(next), '/views/main')!;
    expect(b).not.toBe(a);
    expect(b.sections[1].tracks[0].label).toBe('Domains');
    expect(b.sections[1].follow!.path).toBe(a.sections[1].follow!.path);
    expect(b.sections[1].tracks[0].layers[0].path).toBe(a.sections[1].tracks[0].layers[0].path);
  });
});

describe('resolved view: what does not resolve', () => {
  function view(build: (b: ReturnType<typeof createSeqSpec>) => void) {
    const b = createSeqSpec();
    build(b);
    return resolveView(resolveSpec(b.build()), '/views/v')!;
  }

  it('a section of no system, a layer of an unknown type, a layer whose type does not draw its data', () => {
    const v = view((b) => {
      const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MK' });
      const vals = b.data.values('vals', { target: s, values: new Float32Array([1, 2]) });
      b.views.view('v', (vb) => {
        vb.section('none', { of: '/data/missing' }).track('t').layer('ranges', { type: 'ranges', data: vals });
        const t = vb.section('s', { of: s }).track('t');
        t.layer('ranges', { type: 'ranges', data: vals });
        t.layer('arcs', { type: 'arcs' as 'ranges', data: vals });
      });
    });
    expect(codes(v.issues)).toEqual(['seq.view.section.system', 'seq.view.layer.shape', 'seq.view.layer.type']);
    expect(v.sections[0].tracks[0].layers[0].placed).toBe(false);
    expect(v.sections[1].tracks[0].layers.map((l) => l.placed)).toEqual([false, false]);
  });
});

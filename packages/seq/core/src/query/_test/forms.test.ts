/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type FormResolver, parseAddress, printAddress, selectAddress } from '@molstar/core/repr/tree';
import { describe, expect, it, vi } from 'vitest';
import { type Evaluation, evaluateQuery, QueryFormKinds, queryResolvers } from '../forms.js';
import { QueryIssues } from '../issues.js';
import { fromNamed } from '../result.js';
import { type Location, rangesLocation } from '../location.js';
import type { QueryModel } from '../model.js';
import { Annotations, Model, Spec } from './fixture.js';

const run = (address: string, limit = 100, model: QueryModel = Model) => evaluateQuery(Spec, model, address, { limit });

/** A location's ranges per system, or the error codes. */
function loc(address: string): [string, number[]][] | string[] {
  const r = run(address);
  if (r.result?.type !== 'location') return codes(r);
  return r.result.location.parts.map((p): [string, number[]] => [p.system, [...p.ranges]]);
}

function keys(address: string): (string | number)[] | string[] {
  const r = run(address);
  if (r.result?.type !== 'elements') return codes(r);
  return r.result.elements.keys.map((k) => k.key);
}

const codes = (r: Evaluation) => r.issues.map((i) => i.code);
const p53 = '/data/p53';

describe('query forms: location forms', () => {
  it('range: 1-based labels by default, half-open indices inside (EC-1)', () => {
    expect(loc('/data/p53/(range :from 94 :to 292)')).toEqual([[p53, [93, 292]]]);
    expect(loc('/data/p53/(range :from 94 :to 292 :numbering index)')).toEqual([[p53, [94, 293]]]);
    expect(loc('/data/p53/(range :from 1 :to 1)')).toEqual([[p53, [0, 1]]]);
    // a missing end is the system's
    expect(loc('/data/p53/(range :from 390)')).toEqual([[p53, [389, 393]]]);
    expect(loc('/data/p53/(range :to 3)')).toEqual([[p53, [0, 3]]]);
    expect(loc('/data/empty/(range)')).toEqual([]);
  });

  it('range: in a declared numbering, one that skips zero (EC-3)', () => {
    expect(loc('/data/offset/(range :from -1 :to 1 :numbering skip)')).toEqual([['/data/offset', [2, 4]]]);
    // 0 is not a label: a start snaps to 1
    expect(loc('/data/offset/(range :from 0 :to 3 :numbering skip)')).toEqual([['/data/offset', [3, 6]]]);
  });

  it('range: a missing end snaps inward, with an info; nothing left is a warning (EC-6)', () => {
    const r = run('/data/ab/(range :from 84 :to 200 :numbering kabat)');
    expect(r.result?.type === 'location' && [...r.result.location.parts[0].ranges]).toEqual([8, 12]);
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([['seq.query.label.snapped', 'info']]);
    const e = run('/data/ab/(range :from 84 :to 84 :numbering kabat)');
    expect(e.result?.type === 'location' && e.result.location.parts).toEqual([]);
    expect(codes(e)).toEqual(['seq.query.label.snapped', 'seq.query.range.empty']);
    expect(loc('/data/p53/(range :from -10 :to 2)')).toEqual([[p53, [0, 2]]]);
  });

  it('range: backwards is an error, in index order for labels out of order (EC-7)', () => {
    expect(loc('/data/p53/(range :from 20 :to 10)')).toEqual(['seq.query.range.backwards']);
    // IMGT: 112A sits before 112 in index order, so 112 to 112A is backwards, 111A to 112 is not
    expect(loc('/data/ab/(range :from 112 :to 112A :numbering imgt)')).toEqual(['seq.query.range.backwards']);
    expect(loc('/data/ab/(range :from 111A :to 112 :numbering imgt)')).toEqual([['/data/ab', [4, 7]]]);
    // snapped ends are compared as labels
    expect(loc('/data/ab/(range :from 84 :to 81 :numbering kabat)')).toEqual(['seq.query.range.backwards']);
    const hint = run('/data/p53/(range :from 20 :to 10)').issues[0].hint;
    expect(hint).toBe(QueryIssues.codes['seq.query.range.backwards'].hint);
  });

  it('range: wraps on a cyclic system as two ranges, short of a whole lap (EC-9, EC-10, EC-14)', () => {
    expect(loc('/data/plasmid/(range :from 91 :to 10)')).toEqual([['/data/plasmid', [0, 10, 90, 100]]]);
    expect(loc('/data/plasmid/(range :from 11 :to 10)')).toEqual(['seq.query.range.lap']);
    expect(loc('/data/plasmid/(range :from 12 :to 10)')).toEqual([['/data/plasmid', [0, 10, 11, 100]]]);
  });

  it('range: no insertion code; labels outside int32 refused (EC-4, EC-5)', () => {
    expect(loc('/data/p53/(range :from 1 :to 5 :ins A)')).toEqual(['seq.query.label.insertion']);
    expect(loc('/data/p53/(range :from 1 :to 3000000000)')).toEqual(['seq.query.label.int32']);
    expect(loc('/data/p53/(range :from 1.5)')).toEqual(['seq.query.label.syntax']);
  });

  it('at: one position by label, with insertion codes (EC-5)', () => {
    expect(loc('/data/ab/(at :label 82B :numbering kabat)')).toEqual([['/data/ab', [4, 5]]]);
    expect(loc('/data/ab/(at :label 82 :ins B :numbering kabat)')).toEqual([['/data/ab', [4, 5]]]);
    for (const none of ['""', '?', '.'])
      expect(loc(`/data/ab/(at :label 82 :ins ${none} :numbering kabat)`)).toEqual([['/data/ab', [2, 3]]]);
    expect(loc('/data/ab/(at :label 82A :ins B :numbering kabat)')).toEqual(['seq.query.label.insertion']);
    expect(loc('/data/ab/(at :label 82 :ins 1 :numbering kabat)')).toEqual(['seq.query.label.syntax']);
    expect(loc('/data/p53/(at :label 94)')).toEqual([[p53, [93, 94]]]);
  });

  it('at: an unknown single label is an error naming the labels that exist (EC-6)', () => {
    const r = run('/data/ab/(at :label 84 :numbering kabat)');
    expect(codes(r)).toEqual(['seq.query.label.unknown']);
    expect(r.issues[0]).toMatchObject({
      severity: 'error',
      address: '/data/ab',
      hint: 'labels near it: 82B 82C 83 85 86 87',
    });
    expect(r.result).toBeUndefined();
  });

  it("numberings: a unit's own before its molecule's; an unknown one lists those there are (EC-8)", () => {
    expect(loc('/data/abunit/(at :label 1 :numbering auth)')).toEqual([['/data/abunit', [0, 1]]]);
    expect(loc('/data/abunit/(at :label 82A :numbering kabat)')).toEqual([['/data/abunit', [1, 2]]]);
    // the molecule's 80 is outside the unit
    expect(loc('/data/abunit/(at :label 80 :numbering kabat)')).toEqual(['seq.query.label.unknown']);
    const r = run('/data/abunit/(range :numbering pdb)');
    expect(codes(r)).toEqual(['seq.query.numbering.unknown']);
    expect(r.issues[0].hint).toBe('numberings: index, ordinal, auth, kabat, imgt');
  });

  it('all: a whole system; a form over several systems gives a part each', () => {
    expect(loc('/data/plasmid/(all)')).toEqual([['/data/plasmid', [0, 100]]]);
    expect(loc('/data/empty/(all)')).toEqual([]);
  });

  it('refuses a prefix that is not a system, unknown params and a missing label', () => {
    expect(loc('/data/genes/(range :from 1 :to 2)')).toEqual(['seq.query.form.prefix']);
    expect(loc('/data/p53/(range :from 1 :until 2)')).toEqual(['seq.query.form.param']);
    expect(loc('/data/p53/(at)')).toEqual(['seq.query.form.param']);
    expect(loc('/data/p53/(all :x 1)')).toEqual(['seq.query.form.param']);
    expect(loc('/data/p53/(at :label 94 :numbering [a])')).toEqual(['seq.query.form.param']);
    // a location form's result is not nodes: nothing may follow it but a resolver
    expect(codes(run('/data/p53/(all)/domains'))).toEqual(['repr.tree.address.resolver']);
  });
});

describe('query forms: elements and nodes', () => {
  it('item: by id or by fields; a numeric key never matches a string key', () => {
    expect(keys('/data/p53/domains/(item :id DBD)')).toEqual(['DBD']);
    expect(keys('/data/p53/domains/(item :label DNA-binding)')).toEqual(['DBD']);
    expect(keys('/data/p53/domains/(item :type domain)')).toEqual(['DBD', 'TAD', 'TET']);
    expect(keys('/data/p53/domains/(item :type [domain site] :label [one site])')).toEqual([1, '1']);
    expect(keys('/data/p53/domains/(item :id 1)')).toEqual([1]);
    expect(keys('/data/p53/domains/(item :id "1")')).toEqual(['1']);
    expect(keys('/data/p53/domains/(item :id nope)')).toEqual([]);
    expect(keys('/data/p53/domains/(item)')).toEqual(['seq.query.form.param']);
    expect(keys('/data/p53/(item :id DBD)')).toEqual(['seq.query.form.prefix']);
  });

  it('items: every item, or those overlapping a location on the target', () => {
    expect(keys('/data/p53/domains/(items)')).toEqual([1, '1', 'DBD', 'TAD', 'TET']);
    expect(keys('/data/p53/domains/(items :overlapping (range :from 50 :to 100))')).toEqual(['DBD', 'TAD']);
    // an address as the location; a zero-length item where its point is inside
    expect(keys('/data/p53/domains/(items :overlapping /data/p53/(at :label 301))')).toEqual([1]);
    // a location on another system overlaps nothing here, and says so (no projection yet)
    const other = run('/data/p53/domains/(items :overlapping /data/plasmid/(all))');
    expect(other.result?.type === 'elements' && other.result.elements.keys).toEqual([]);
    expect(other.issues.map((i) => [i.code, i.severity, i.address])).toEqual([
      ['seq.query.location.unreachable', 'warning', '/data/plasmid'],
    ]);
    expect(keys('/data/p53/domains/(items :overlapping /data/p53/domains/(item :id DBD))')).toEqual([
      'seq.query.form.type',
    ]);
  });

  it('node: a dynamic node held under the prefix, addressed by its form', () => {
    const r = run('/data/genes/(node :key T1)');
    expect(r.result?.type === 'nodes' && r.result.nodes.map((m) => [printAddress(m.address), m.node.kind])).toEqual([
      ['/data/genes/(node :key T1)', 'gene'],
    ]);
    const missing = run('/data/genes/(node :key T9)');
    expect(missing.result?.type === 'nodes' && missing.result.nodes).toEqual([]);
    expect(missing.issues.map((i) => [i.code, i.severity])).toEqual([['seq.query.node.missing', 'warning']]);
    expect(codes(run('/data/genes/(node)'))).toEqual(['seq.query.form.param']);
  });

  it('an address of ids gives nodes', () => {
    const r = run('/views/main/*');
    expect(r.result?.type === 'nodes' && r.result.nodes.length).toBe(3);
  });
});

describe('query forms: conversions', () => {
  it('extent: elements to a location, a wrapped item as two ranges (EC-9)', () => {
    expect(loc('/data/p53/domains/(extent :of (item :type domain))')).toEqual([[p53, [0, 61, 93, 292, 324, 355]]]);
    expect(loc('/data/plasmid/features/(extent :of (item :id ori))')).toEqual([['/data/plasmid', [0, 5, 90, 100]]]);
    expect(loc('/(extent :of /data/p53/domains/(item :id TAD))')).toEqual([[p53, [0, 61]]]);
    expect(loc('/data/p53/(extent :of (all))')).toEqual(['seq.query.form.type']);
  });

  it('round trip: the items overlapping an extent hold the elements', () => {
    expect(keys('/data/plasmid/features/(items :overlapping (range :from 98 :to 2))')).toEqual(['end', 'ori']);
  });
});

describe('query forms: combinators over one type', () => {
  it('and, or, not over locations', () => {
    expect(loc('/data/p53/(or :of [(range :from 1 :to 10) (range :from 5 :to 20) (at :label 30)])')).toEqual([
      [p53, [0, 20, 29, 30]],
    ]);
    expect(loc('/data/p53/(and :of [(range :from 1 :to 10) (range :from 5 :to 20)])')).toEqual([[p53, [4, 10]]]);
    expect(loc('/data/p53/(not :of (range :from 2 :to 392))')).toEqual([[p53, [0, 1, 392, 393]]]);
  });

  it('and, or, not over elements and nodes', () => {
    expect(keys('/data/p53/domains/(and :of [(item :type domain) (items :overlapping (range :to 100))])')).toEqual([
      'DBD',
      'TAD',
    ]);
    expect(keys('/data/p53/domains/(not :of (item :type domain))')).toEqual([1, '1']);
    expect(keys('/data/p53/domains/(or :of [(item :id TET) /data/plasmid/features/(item :id ori)])')).toEqual([
      'TET',
      'ori',
    ]);
    const n = run('/views/main/(not :of (or :of [/views/main/protein (view)]))');
    expect(n.result?.type === 'nodes' && n.result.nodes.map((m) => printAddress(m.address))).toEqual([
      '/views/main/gene',
    ]);
  });

  it('a combinator that mixes types is an issue, not a guess', () => {
    expect(codes(run('/data/p53/domains/(or :of [(item :id TET) (extent :of (item :id TET))])'))).toEqual([
      'seq.query.form.type',
    ]);
    expect(codes(run('/data/p53/(or :of (all))'))).toEqual(['seq.query.form.param']);
  });
});

describe('query forms: staged and bounded evaluation', () => {
  it('over data read by region, says what is pending and never calls it empty', () => {
    const r = run('/data/p53/peaks/(items :overlapping (range :from 1 :to 100))');
    expect(r.result?.type === 'elements' && r.result.elements.keys.map((k) => k.key)).toEqual(['a', 'b']);
    expect(r.complete).toBe(false);
    expect(r.pending.parts.map((p) => [p.system, [...p.ranges]])).toEqual([[p53, [50, 100]]]);
    const loaded = run('/data/p53/peaks/(items :overlapping (range :from 1 :to 50))');
    expect([loaded.complete, loaded.pending.parts]).toEqual([true, []]);
    // nested: the pending regions of an operand reach the outer result
    const nested = run('/data/p53/peaks/(extent :of (item :id b))');
    expect(nested.complete).toBe(false);
    expect((nested.pending as Location).parts.length).toBe(1);
  });

  it('a result over the limit is an issue, never a truncation', () => {
    const r = run('/data/p53/domains/(items)', 4);
    expect([r.result, codes(r)]).toEqual([undefined, ['repr.tree.address.limit_exceeded']]);
    expect(run('/data/p53/domains/(items)', 5).result?.type).toBe('elements');
  });

  it('stops on an aborted signal inside a long read', () => {
    const many = Array.from({ length: 3000 }, (_, i) => ({ key: i, start: i, end: i + 1 }));
    const c = new AbortController();
    const model: QueryModel = {
      ...Model,
      annotation: (a) =>
        a === '/data/p53/domains'
          ? { address: a, target: p53, fields: [], items: () => (c.abort(), { items: many }) }
          : Model.annotation(a),
    };
    const r = evaluateQuery(Spec, model, '/data/p53/domains/(items)', { limit: Infinity, signal: c.signal });
    expect(codes(r)).toEqual(['repr.tree.address.cancelled']);
  });

  /** An annotation of `n` items whose read runs `onRead` (aborts a signal, moves a clock). */
  const readsOf = (n: number, onRead: () => void): QueryModel => {
    const rows = Array.from({ length: n }, (_, i) => ({ key: i, start: i, end: i + 1 }));
    return {
      ...Model,
      annotation: (a) =>
        a === '/data/p53/domains'
          ? { address: a, target: p53, fields: [], items: () => (onRead(), { items: rows }) }
          : Model.annotation(a),
    };
  };

  it('nested forms share one visit counter: two reads under CancelEvery each still stop', () => {
    const c = new AbortController();
    const model = readsOf(600, () => c.abort());
    const one = evaluateQuery(Spec, model, '/data/p53/domains/(items)', { limit: Infinity, signal: c.signal });
    expect(codes(one)).toEqual([]);
    const c2 = new AbortController();
    const both = evaluateQuery(
      Spec,
      readsOf(600, () => c2.abort()),
      '/data/p53/domains/(or :of [(items) (items)])',
      {
        limit: Infinity,
        signal: c2.signal,
      },
    );
    expect(codes(both)).toEqual(['repr.tree.address.cancelled']);
  });

  it('stops at a deadline passed during a read', () => {
    let now = 0;
    vi.spyOn(performance, 'now').mockImplementation(() => now);
    const r = evaluateQuery(
      Spec,
      readsOf(2000, () => (now = 10)),
      '/data/p53/domains/(items)',
      {
        limit: Infinity,
        deadline: 5,
      },
    );
    expect(codes(r)).toEqual(['repr.tree.address.deadline_exceeded']);
  });

  it('two staged segments: the result is incomplete (the union of their pending regions is the core’s)', () => {
    const pend: FormResolver = ({ prefix }) => ({
      result: prefix,
      complete: false,
      pending: rangesLocation('/data/x', [0, 1]),
    });
    const t = queryResolvers(Model, { pend });
    const s = selectAddress(Spec, '/data/p53/(pend)/peaks/(items :overlapping (range :to 100))', {
      limit: 9,
      resolvers: t,
    });
    expect(s.complete).toBe(false);
    // today the selection keeps the last segment's pending regions
    expect((s.pending as Location).parts.map((p) => [p.system, [...p.ranges]])).toEqual([[p53, [50, 100]]]);
  });
});

describe('query forms: chains, rects, labels (review round)', () => {
  const chainA = '/data/1tup/(chain :id A)';

  it('chain: a data-made system addressed by its form, a prefix location forms read', () => {
    const r = run('/data/1tup/(chain :id A)');
    expect(r.result?.type === 'nodes' && r.result.nodes.map((m) => printAddress(m.address))).toEqual([chainA]);
    expect(loc('/data/1tup/(chain :id A)/(at :label 94A :numbering auth)')).toEqual([[chainA, [3, 4]]]);
    const b = run('/data/1tup/(chain :id B)');
    expect(b.issues.map((i) => i.code)).toEqual(['seq.query.node.missing']);
    expect(codes(run('/data/1tup/(chain :id A :x 1)'))).toEqual(['seq.query.form.param']);
    // an annotation on a chain: its nested location reads the chain as its prefix
    expect(keys('/data/1tup/sites/(items :overlapping (range :from 94A :to 95 :numbering auth))')).toEqual(['s1']);
  });

  it('row: an alignment’s rows by key, positions of its rows system; a rect over rows × columns', () => {
    expect(loc('/data/aln/(row :key m)')).toEqual([['/data/aln/(rows)', [1, 2]]]);
    expect(loc('/data/aln/(row :key [r h])')).toEqual([['/data/aln/(rows)', [0, 1, 2, 3]]]);
    expect(codes(run('/data/aln/(row :key x)'))).toEqual(['seq.query.row.unknown']);
    const r = run('/data/aln/(rect :a /data/aln/(row :key [h m]) :b (range :from 3 :to 6))');
    expect(
      r.result?.type === 'location' &&
        r.result.location.rects.map((x) => [x.a.system, [...x.a.ranges], x.b.system, [...x.b.ranges]]),
    ).toEqual([['/data/aln/(rows)', [0, 2], '/data/aln', [2, 6]]]);
  });

  it('rects subtract: a location less `(not :of …)` within `and`, a rect less a rect exact', () => {
    const r = run(
      '/data/p53/(and :of [(rect :a (range :from 1 :to 10) :b (range :from 1 :to 10)) (not :of (rect :a (range :from 4 :to 6) :b (range :from 4 :to 6)))])',
    );
    expect(r.issues).toEqual([]);
    const rects = r.result?.type === 'location' ? r.result.location.rects : [];
    let cells = 0;
    for (const x of rects) {
      let w = 0;
      let h = 0;
      for (let k = 0; k < x.a.ranges.length; k += 2) w += x.a.ranges[k + 1] - x.a.ranges[k];
      for (let k = 0; k < x.b.ranges.length; k += 2) h += x.b.ranges[k + 1] - x.b.ranges[k];
      cells += w * h;
    }
    expect(cells).toBe(100 - 9);
    // parts as before: a range less a range
    expect(loc('/data/p53/(and :of [(range :from 1 :to 10) (not :of (range :from 4 :to 6))])')).toEqual([
      [p53, [0, 3, 6, 10]],
    ]);
  });

  it('rect: two locations on two axes, each one system; sides never wrap (EC-10)', () => {
    const r = run('/data/p53/(rect :a (range :from 10 :to 25) :b /data/plasmid/(range :from 40 :to 55))');
    expect(
      r.result?.type === 'location' &&
        r.result.location.rects.map((x) => [x.a.system, [...x.a.ranges], x.b.system, [...x.b.ranges]]),
    ).toEqual([[p53, [9, 25], '/data/plasmid', [39, 55]]]);
    expect(codes(run('/data/plasmid/(rect :a (range :from 91 :to 10) :b (all))'))).toEqual(['seq.query.range.wrap']);
    // the same range wraps outside a rect
    expect(loc('/data/plasmid/(range :from 91 :to 10)')).toEqual([['/data/plasmid', [0, 10, 90, 100]]]);
    expect(codes(run('/data/p53/(rect :a (all) :b /data/p53/domains/(item :id TAD))'))).toEqual([
      'seq.query.form.type',
    ]);
    expect(codes(run('/data/p53/(rect :a (or :of [(all) /data/plasmid/(all)]) :b (all))'))).toEqual([
      'seq.query.form.type',
    ]);
    expect(codes(run('/data/p53/(not :of (rect :a (all) :b (all)))'))).toEqual(['seq.query.form.unsupported']);
    expect(codes(run('/data/p53/domains/(items :overlapping /data/p53/(rect :a (all) :b (all)))'))).toEqual([
      'seq.query.form.unsupported',
    ]);
  });

  it('snaps first, then one index-order rule: a cyclic snapped range wraps; IMGT snapped pairs (EC-6, EC-7)', () => {
    // no label at 5 and 96: 96 starts at 97, 5 ends at 4
    expect(loc('/data/plasmid/(range :from 96 :to 5 :numbering gapped)')).toEqual([['/data/plasmid', [0, 4, 96, 100]]]);
    // IMGT (111, 111A, 112A, 112, 113): 111B starts at 112A, 112B ends at 112
    expect(loc('/data/ab/(range :from 111B :to 112B :numbering imgt)')).toEqual([['/data/ab', [5, 7]]]);
    expect(loc('/data/ab/(range :from 112B :to 112B :numbering imgt)')).toEqual([]);
    expect(loc('/data/ab/(range :from 112B :to 111B :numbering imgt)')).toEqual(['seq.query.range.backwards']);
    // a snapped end over a hole: 84 to 84 holds nothing (82C and 85 around a position with no label)
    expect(loc('/data/ab/(range :from 84 :to 84 :numbering kabat)')).toEqual([]);
  });

  it('.n labels are quoted; 100.10 is not 100.1 (EC-2, EC-3)', () => {
    expect(loc('/data/offset/(at :label "100.1" :numbering dot)')).toEqual([['/data/offset', [1, 2]]]);
    expect(loc('/data/offset/(at :label "100.10" :numbering dot)')).toEqual([['/data/offset', [2, 3]]]);
    const bare = run('/data/offset/(at :label 100.1 :numbering dot)');
    expect([codes(bare), bare.issues[0].hint]).toEqual([['seq.query.label.syntax'], 'a .n label is quoted: "100.1"']);
    expect(loc('/data/offset/(range :from "100.1" :to 101 :numbering dot)')).toEqual([['/data/offset', [1, 4]]]);
  });

  it('a boundary at the system end overlaps a range ending there, so (all) holds it', () => {
    expect(keys('/data/plasmid/features/(items :overlapping (all))')).toEqual(['end', 'gene', 'ori']);
    expect(keys('/data/plasmid/features/(items :overlapping (range :from 50 :to 99))')).toEqual(['ori']);
  });

  it('item: an unknown field is an issue', () => {
    const r = run('/data/p53/domains/(item :colour red)');
    expect([codes(r), r.issues[0].hint]).toEqual([['seq.query.form.param'], 'params: id, label, type']);
  });

  it('refuses an operand naming a ref: refs are not resolved yet', () => {
    expect(codes(run('/data/p53/domains/(items :overlapping /refs/dbd)'))).toEqual(['seq.query.ref.unresolved']);
  });

  it("refuses an item's extent outside int32", () => {
    const model: QueryModel = {
      ...Model,
      annotation: (a) =>
        a === '/data/p53/domains'
          ? { address: a, target: p53, fields: [], items: () => ({ items: [{ key: 'x', start: 0, end: 2 ** 31 }] }) }
          : Model.annotation(a),
    };
    expect(codes(run('/data/p53/domains/(extent :of (item :id x))', 9, model))).toEqual(['seq.query.position.invalid']);
  });
});

describe('query forms: the table', () => {
  it('registers every kind, beside a caller’s own', () => {
    const extra = { mine: () => ({ result: { kind: 'other' as const, type: 'mine', value: 1 } }) };
    const t = queryResolvers(Model, extra);
    expect(Object.keys(t).sort()).toEqual([...QueryFormKinds, 'mine'].sort());
    const s = selectAddress(Spec, '/data/p53/(mine)', { limit: 9, resolvers: t });
    expect(s.result).toEqual({ kind: 'other', type: 'mine', value: 1 });
    // a nested result of a type the queries do not know is an issue
    expect(selectAddress(Spec, '/data/p53/(not :of (mine))', { limit: 9, resolvers: t }).issues[0].code).toBe(
      'seq.query.form.type',
    );
  });

  it('refuses a caller resolver named as a query form: a form is never replaced silently', () => {
    expect(() => queryResolvers(Model, { range: () => ({}) })).toThrow(/clashes/);
  });

  it('reads a Named as a result only when its value has the type’s shape', () => {
    expect(fromNamed({ kind: 'other', type: 'location', value: 5 })).toBeUndefined();
    expect(fromNamed({ kind: 'other', type: 'elements', value: { keys: [{ owner: 1 }] } })).toBeUndefined();
    expect(fromNamed({ kind: 'other', type: 'location', value: rangesLocation('/a', [0, 1]) })?.type).toBe('location');
  });

  it('prints a string that looks like a number quoted, so it reads back a string', () => {
    const form = { kind: 'item', params: { id: '1' } };
    const text = printAddress(['data', 'p53', 'domains', form]);
    expect(text).toBe('/data/p53/domains/(item :id "1")');
    expect(keys(text)).toEqual(['1']);
    expect(printAddress(['x', { kind: 'at', params: { label: '100.10' } }])).toBe('/x/(at :label "100.10")');
  });

  it('reads and prints the forms canonically', () => {
    const a = '/data/p53/domains/(items :overlapping (range :from 94 :numbering uniprot :to 94A))';
    const p = parseAddress(a).pattern!;
    expect(printAddress(p.segments)).toBe(
      '/data/p53/domains/(items :overlapping (range :from 94 :numbering uniprot :to 94A))',
    );
  });

  it('every code used is cataloged, and every cataloged code is used', async () => {
    const { readFileSync, readdirSync } = await import('node:fs');
    const dir = new URL('..', import.meta.url);
    const source = readdirSync(dir)
      .filter((f) => f.endsWith('.ts') && f !== 'issues.ts')
      .map((f) => readFileSync(new URL(f, dir), 'utf8'))
      .join('\n');
    const used = new Set(source.match(/'seq\.query\.[a-z0-9_.]+'/g)?.map((s) => s.slice(1, -1)));
    const cataloged = Object.keys(QueryIssues.codes);
    for (const c of used) expect(cataloged, c).toContain(c);
    for (const c of cataloged) expect(used.has(c), c).toBe(true);
  });

  it('reads each annotation once per form', () => {
    const before = Annotations.domains.reads;
    run('/data/p53/domains/(item :type domain)');
    expect(Annotations.domains.reads - before).toBe(1);
  });
});

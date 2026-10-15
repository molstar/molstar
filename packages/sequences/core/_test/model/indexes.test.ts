/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Build-once indexes and loci (M3): "look up, never search" (seqspec:model#locations-and-loci). The
// scaling tests compare a query's cost at two sizes (the best of a few runs), never a time against a constant; they
// read the clock, so they are the slow tier's.

import { describe, expect, it } from 'vitest';
import {
  type AlignmentId,
  type Annotation,
  Annotations,
  ColumnRole,
  ColumnRoles,
  type Issue,
  type ItemIndex,
  type ItemLoci,
  Loci,
  OrderedSet,
  Pairs,
  Providers,
  type RowIndex,
  RowLocation,
  type SeqIndex,
  SeqLocation,
  type SeqLoci,
  type SeqSet,
  SeqSets,
  type UnitId,
} from '../../src/model/index.js';
import { brs1 } from './fixtures.js';

const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
const strings = (set: SeqSet, ids: ArrayLike<number>) => Array.from(ids, (id) => set.strings.get(id));

/** The best per-call time of `run` (which makes `calls` calls) over a few runs, in ms. */
function perCall(run: () => void, calls: number): number {
  let best = Number.POSITIVE_INFINITY;
  for (let k = 0; k < 5; k++) {
    const t = performance.now();
    run();
    best = Math.min(best, (performance.now() - t) / calls);
  }
  return best;
}

/** A visitor that does nothing: the scaling tests time the query. */
const visit = () => {};

/** A seeded generator (mulberry32). */
function random(seed: number) {
  let a = seed;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** `n` chains over one molecule, each bound to its own row of an alignment by an identity mapping. */
function chains(n: number) {
  const SEQ = 'MKVLAAGIVGLPNVGKSTLFNALTKAGIEAANYPFCTIEPNTGVVPMPDP';
  return SeqSets.build({
    molecules: [{ key: 'M', residues: SEQ }],
    units: Array.from({ length: n }, (_, k) => ({ key: `c${k}`, molecule: 'M' })),
    alignments: [{ key: 'aln', flavor: 'aligned', gapped: Array.from({ length: n }, () => `-${SEQ}`) }],
    mappings: Array.from({ length: n }, (_, k) => ({
      key: `m${k}`,
      from: `aln/${k}`,
      to: `c${k}`,
      method: 'identity' as const,
    })),
  });
}

describe('rows slots', () => {
  const records = (n: number) =>
    SeqSets.build({
      molecules: [],
      alignments: [
        {
          key: 'R',
          flavor: 'aligned',
          gapped: [],
          lazy: { rowCount: n, columnCount: 5, rows: new Map(), labelOf: () => null },
        },
      ],
      annotations: [
        {
          key: 'meta',
          geometry: 'rows',
          target: 'R',
          columns: {
            row: Array.from({ length: n }, (_, i) => n - 1 - i),
            species: Array.from({ length: n }, (_, i) => `s${i % 7}`),
          },
        },
      ],
    });

  it('give 50k records a slot each, through an index', () => {
    const { set, issues } = records(50_000);
    expect(issues).toEqual([]);
    const a = set.annotations[0] as Annotation;
    expect([a.count, a.slots.row.length, a.slot[49_999], a.slots.row[0], a.slots.row[49_999]]).toEqual([
      50_000, 50_000, 49_999, 49_999, 0,
    ]);
    expect(Annotations.rowValue(a, set, 3, 'species')).toBe('s3');
  });

  it('build in time linear in the records: 10 times the records cost under 30 times as much', {
    tags: ['slow'],
  }, () => {
    const cost = (n: number) => perCall(() => records(n), 1);
    expect(cost(50_000) / cost(5_000)).toBeLessThan(30);
  });
});

describe('pairs: touching through a running maximum, endUnits', () => {
  it('visit what a brute force over every pair finds, at every position of every unit', () => {
    const next = random(7);
    const n = 300;
    const pick = <T>(xs: readonly T[]) => xs[Math.floor(next() * xs.length)];
    const start: number[] = [];
    const stop: number[] = [];
    const bStart: (number | null)[] = [];
    for (let i = 0; i < n; i++) {
      start.push(Math.floor(next() * 100));
      stop.push(Math.min(110, start[i] + 1 + Math.floor(next() * (next() < 0.1 ? 60 : 6))));
      bStart.push(next() < 0.15 ? null : Math.floor(next() * 89));
    }
    const { set, issues } = SeqSets.build({
      ...brs1(),
      annotations: [
        {
          key: 'p',
          geometry: 'pairs',
          columns: {
            target: Array.from({ length: n }, () => pick(['1brs/A', '1brs/B', '1brs/entity:1'])),
            start,
            stop,
            b_target: Array.from({ length: n }, () => pick(['1brs/D', '1brs/entity:2', '1brs/A'])),
            b_start: bStart,
          },
        },
      ],
    });
    expect(issues).toEqual([]);
    const p = set.annotations[0] as Annotation;
    const pc = p.pairs as NonNullable<Annotation['pairs']>;
    const reaches = (unit: number, row: number, side: number) => {
      const t = side ? pc.bTarget : p.target;
      return t.unit[row] >= 0 ? t.unit[row] === unit : set.units.molecule[unit] === t.molecule[row];
    };
    let checked = 0;
    for (let u = 0; u < set.units.count; u++) {
      for (let i = 0; i <= 110; i++) {
        const got: number[] = [];
        Pairs.touching(p, u as UnitId, i as SeqIndex, (row, side) => got.push(row * 2 + side), null);
        const want: number[] = [];
        for (let r = 0; r < p.count; r++) {
          for (let side = 0; side < 2; side++) {
            const s = side ? pc.bStart[r] : p.start[r];
            const t = side ? pc.bStop[r] : p.stop[r];
            if (reaches(u, r, side) && (s < 0 || (s <= i && i < Math.max(t, s + 1)))) want.push(r * 2 + side);
          }
        }
        expect(
          got.sort((a, b) => a - b),
          `unit ${u} at ${i}`,
        ).toEqual(want);
        checked += want.length;
      }
    }
    expect(checked).toBeGreaterThan(1000);
    // endUnits: an endpoint's unit, or every unit over its molecule.
    const ends = (row: number, side: number) => {
      const out: string[] = [];
      Pairs.endUnits(p, row, side, (u) => out.push(set.strings.get(set.units.key[u])), null);
      return out;
    };
    for (let r = 0; r < 5; r++) {
      const key = set.strings.get(set.units.key[p.target.unit[r]] ?? -1);
      expect(ends(r, 0)).toEqual(p.target.unit[r] >= 0 ? [key] : ['1brs/A', '1brs/B', '1brs/C']);
    }
    expect(ends(-1, 0)).toEqual([]);
  });

  it('cost about the same per hover near the end of a unit with 32 times the pairs', { tags: ['slow'] }, () => {
    const L = 20_000;
    const contacts = (n: number) => {
      const set = SeqSets.build({
        molecules: [
          { key: 'u', residues: 'A'.repeat(L), ownUnit: true },
          { key: 'w', residues: 'AAAA', ownUnit: true },
        ],
        annotations: [
          {
            key: 'c',
            geometry: 'pairs',
            target: 'u/unit',
            columns: {
              start: Array.from({ length: n }, (_, i) => Math.floor((i * L) / n)),
              b_target: Array.from({ length: n }, () => 'w/unit'),
              b_start: Array.from({ length: n }, () => null),
            },
          },
        ],
      }).set;
      return set.annotations[0] as Annotation;
    };
    const cost = (p: Annotation) => {
      return perCall(() => {
        for (let k = 0; k < 2000; k++) Pairs.touching(p, 0 as UnitId, (L - 1 - (k % 100)) as SeqIndex, visit, null);
      }, 2000);
    };
    const small = contacts(2_000);
    const large = contacts(64_000);
    cost(small);
    expect(cost(large) / cost(small)).toBeLessThan(8);
  });
});

describe('toRowLocations through the unit index', () => {
  it("give each of 1,000 chains its own row's residue", () => {
    const { set, issues } = chains(1000);
    expect(issues).toEqual([]);
    const l = SeqLocation.create(set);
    const out = RowLocation.create(set);
    const seen: string[] = [];
    for (const u of [0, 1, 517, 999]) {
      SeqLocation.set(l, u, 7);
      SeqLocation.toRowLocations(l, 0 as AlignmentId, out, (r) => seen.push(`${r.row}:${r.index}`), null);
    }
    expect(seen).toEqual(['0:7', '1:7', '517:7', '999:7']);
  });

  it('cost about the same per hover with 10 times the chains', { tags: ['slow'] }, () => {
    const cost = (n: number) => {
      const set = chains(n).set;
      const l = SeqLocation.create(set);
      const out = RowLocation.create(set);
      const calls = 20_000;
      return perCall(() => {
        for (let k = 0; k < calls; k++) {
          SeqLocation.set(l, k % n, k % 50);
          SeqLocation.toRowLocations(l, 0 as AlignmentId, out, visit, null);
        }
      }, calls);
    };
    cost(100);
    expect(cost(1000) / cost(100)).toBeLessThan(4);
  });
});

describe('items: itemRows, itemOf and remapped item loci', () => {
  const ranges = (ids: (string | null)[]) =>
    SeqSets.build({
      ...brs1(),
      annotations: [
        { key: 'f', geometry: 'ranges', target: '1brs/A', columns: { start: ids.map((_, i) => i % 100), id: ids } },
      ],
    }).set;

  it("give an item's rows and an item by its key", () => {
    const set = ranges(['a', null, 'a', 'b', 'a']);
    const f = set.annotations[0] as Annotation;
    const byItem = Annotations.itemRows(f);
    expect([Array.from(byItem.start), Array.from(byItem.rows)]).toEqual([
      [0, 3, 4, 5],
      [0, 2, 4, 1, 3],
    ]);
    expect(Annotations.itemRows(f)).toBe(byItem);
    expect([Annotations.itemOf(f, set.strings.find('b')), Annotations.itemOf(f, set.strings.find('zz'))]).toEqual([
      2, -1,
    ]);
  });

  it('remap item loci in time linear in the items', { tags: ['slow'] }, () => {
    const cost = (n: number) => {
      const ids = Array.from({ length: n }, (_, i) => `k${i}`);
      const from = ranges(ids);
      const to = ranges([...ids].reverse());
      const loci: ItemLoci = {
        kind: 'item-loci',
        set: from,
        elements: [{ annotation: 0 as never, items: OrderedSet.ofRange(0 as ItemIndex, n as ItemIndex) }],
      };
      let out: ItemLoci | null = null;
      const t = perCall(() => {
        out = Loci.remap(loci, from, to);
      }, 1);
      expect(OrderedSet.size((out as unknown as ItemLoci).elements[0].items)).toBe(n);
      return t;
    };
    cost(2000);
    expect(cost(16_000) / cost(2000)).toBeLessThan(24);
  });
});

describe('loci and labels', () => {
  it("forEachLocation sets the out location's set", () => {
    const v1 = SeqSets.build(brs1()).set;
    const v2 = SeqSets.withChanges(v1, { upsert: { molecules: [{ key: 'x', residues: 'MK' }] } }).set;
    const loci: SeqLoci = {
      kind: 'seq-loci',
      set: v2,
      elements: [{ unit: 0 as UnitId, indices: OrderedSet.ofRange(3 as SeqIndex, 5 as SeqIndex) }],
    };
    const sets: SeqSet[] = [];
    Loci.forEachLocation(loci, SeqLocation.create(v1), (l) => sets.push(l.set), null);
    expect(sets.every((s) => s === v2) && sets.length).toBe(2);
  });

  it('leave the string table alone over 1,000 labels of rows not held', () => {
    const { set } = SeqSets.build({
      molecules: [],
      alignments: [
        {
          key: 'big',
          flavor: 'aligned',
          gapped: [],
          lazy: {
            rowCount: 10_000,
            columnCount: 3,
            rows: new Map([[0, { label: 'held', gapped: 'MKV' }]]),
            labelOf: (r) => `row${r}`,
          },
        },
      ],
    });
    const a = set.alignments[0];
    if (!a) throw new Error('not built');
    const size = set.strings.size;
    let text = 0;
    for (let r = 1; r <= 1000; r++) {
      if (a.labelOf(r as RowIndex) >= 0) text = -1;
      text += a.labelText(r as RowIndex)?.length ?? 0;
    }
    expect([set.strings.size, a.labelText(5 as RowIndex), a.labelText(0 as RowIndex), text > 0]).toEqual([
      size,
      'row5',
      'held',
      true,
    ]);
  });
});

describe("column roles (the model's half)", () => {
  it("read a geometry's own columns, keep descriptive ones of other geometries as fields, refuse placing ones", () => {
    expect([ColumnRoles.rows.row, ColumnRoles.ranges.alt, ColumnRoles.pairs.b_beg]).toEqual([
      ColumnRole.target,
      ColumnRole.descriptive,
      ColumnRole.position,
    ]);
    const { set, issues } = SeqSets.build({
      ...brs1(),
      alignments: [{ key: 'msa', flavor: 'aligned', gapped: ['MK', 'MV'] }],
      annotations: [
        {
          key: 'meta',
          geometry: 'rows',
          target: 'msa',
          columns: { row: [0, 1], label: ['a', 'b'], category: ['x', 'y'], value: [1, 2] },
        },
        { key: 'r', geometry: 'ranges', target: '1brs/A', columns: { start: [1], row: [0] } },
        { key: 'v', geometry: 'values', target: '1brs/A', columns: { index: [3], value: [0.5], id: ['q'] } },
      ],
    });
    expect(codes(issues)).toEqual(['model.columns.unknown /annotations/1/columns/row']);
    const [meta, r, v] = set.annotations as Annotation[];
    expect([
      strings(
        set,
        meta.fields.map((f) => f.name),
      ),
      meta.category,
      meta.sparse,
    ]).toEqual([['label', 'category', 'value'], null, null]);
    expect(Annotations.rowValue(meta, set, 1, 'category')).toBe('y');
    expect([r.count, r.fields.length]).toEqual([1, 0]);
    expect([
      strings(
        set,
        v.fields.map((f) => f.name),
      ),
      Array.from(v.itemKey),
    ]).toEqual([['id'], [-1]]);
  });
});

describe('coverage through an ambiguous path', () => {
  it('is left out with model.coverage.ambiguous, not as unmapped', () => {
    const barstar = brs1().molecules[1].residues as string;
    const { set, issues } = SeqSets.build({
      ...brs1(),
      molecules: [...brs1().molecules, { key: 'P', residues: barstar }],
      mappings: [
        { key: 'm1', from: 'P', to: '1brs/entity:2', method: 'identity' },
        { key: 'm2', from: 'P', to: '1brs/entity:2', method: 'segments', segments: [[0, 1, 88]] },
      ],
    });
    expect(issues).toEqual([]);
    const r = Providers.coverage(set, { key: 'cov', of: ['1brs/D'], reference: 'P' });
    expect(codes(r.issues)).toEqual(['model.coverage.ambiguous /reference']);
  });
});

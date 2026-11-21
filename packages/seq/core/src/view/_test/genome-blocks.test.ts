/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks over genomes in the resolved view (_spec/README.md "Resolved view", Blocks and Places): a block's direction
// and scale, one system's blocks merged only on one direction, a section of several assemblies' groups, the closing
// run's count, the built-in fields' domains, the places index (`placesOf`) and units as `{ system, position }`.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  blockPositionAt,
  blockUnitOf,
  placeOfUnit,
  placesOf,
  type ResolvedSection,
  resolveView,
  sameBlocks,
  unitOfPlace,
} from '../index.js';

beforeAll(async () => {
  await loadGenome();
});

/** One system shown as written blocks (`[range, direction, scale]` each) in a section. */
function piecesSpec(
  blocks: readonly (readonly [readonly [number, number], 'section' | 'reversed', number])[],
  dir = 'genomic',
): Node {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(250) });
  b.views.view('v', (v) =>
    v.section('a', { of: s, direction: dir } as never, (sec) => {
      blocks.forEach(([range, direction, scale], i) =>
        sec.block(`b${i}`, { of: s, range: [...range], direction, scale } as never),
      );
    }),
  );
  return b.build();
}

const sectionOf = (spec: Node, config = {}): ResolvedSection =>
  resolveView(resolveSpec(spec, { config }), '/views/v')!.sections[0];

describe('a block’s direction and scale (EC-206)', () => {
  it('reads `reversed` against the section, its units its positions times its scale', () => {
    const s = sectionOf(
      piecesSpec([
        [[0, 100], 'section', 1],
        [[100, 200], 'reversed', 1],
        [[200, 300], 'section', 10],
      ]),
    );
    const b = s.blocks!;
    expect(b.map((x) => [x.axis.dir, x.scale ?? 1, x.offset])).toEqual([
      [1, 1, 0],
      [-1, 1, 100],
      [1, 10, 200],
    ]);
    // blocks of one system that read two ways are on no system: units, the section's address
    expect([s.axis.system, s.axis.length]).toEqual(['/views/v/a', 1200]);
  });

  it('flips with the section: a reversed section turns every block, a reversed block reads forward', () => {
    const s = sectionOf(
      piecesSpec(
        [
          [[0, 100], 'section', 1],
          [[100, 200], 'reversed', 1],
        ],
        'reversed',
      ),
    );
    expect(s.blocks!.map((x) => x.axis.dir)).toEqual([-1, 1]);
  });

  it('merges one system’s blocks into one axis only when they read one way at scale 1', () => {
    const one = sectionOf(
      piecesSpec([
        [[0, 100], 'section', 1],
        [[300, 400], 'section', 1],
      ]),
    );
    expect([one.axis.system, one.axis.ranges]).toEqual(['/data/s', [0, 100, 300, 400]]);
    const scaled = sectionOf(
      piecesSpec([
        [[0, 100], 'section', 1],
        [[300, 400], 'section', 2],
      ]),
    );
    expect(scaled.axis.system).toBe('/views/v/a');
  });

  it('keeps a reversed copy apart from a forward one (copies merge on equal direction and scale)', () => {
    const s = sectionOf(
      piecesSpec([
        [[0, 100], 'section', 1],
        [[0, 100], 'reversed', 1],
        [[0, 100], 'section', 1],
      ]),
    );
    expect(s.blocks!.map((x) => [x.label, x.axis.dir])).toEqual([
      ['b0 ×2', 1],
      ['b1', -1],
    ]);
  });
});

const peach = [
  { name: 'Pp01', length: 40 },
  { name: 'Pp02', length: 30 },
];
const grape = [
  { name: 'chr2', length: 20 },
  { name: 'chr1', length: 50 },
];

function twoAssemblies(): Node {
  const b = createSeqSpec();
  b.data.assembly('peach', assembly('peach', peach));
  b.data.assembly('grape', assembly('grape', grape));
  b.views.view('v', (v) => genome(v, 'two', ['/data/peach', '/data/grape']));
  return b.build();
}

describe('a section of several assemblies (EC-210)', () => {
  it('a group each, its blocks a run in order, labelled by its assembly', () => {
    const s = sectionOf(twoAssemblies());
    expect(s.groups).toEqual([
      { assembly: '/data/peach', label: 'peach', first: 0, count: 2 },
      { assembly: '/data/grape', label: 'grape', first: 2, count: 2 },
    ]);
    expect(s.blocks!.map((b) => b.group)).toEqual([0, 0, 1, 1]);
  });

  it('a section of one assembly is one group (its name the ring’s title)', () => {
    const b = createSeqSpec();
    b.data.assembly('g', assembly('GRCh38', grape));
    b.views.view('v', (v) => genome(v, 'g', '/data/g'));
    expect(sectionOf(b.build()).groups?.map((g) => g.label)).toEqual(['GRCh38']);
  });

  it('counts the blocks past `maxAssemblyBlocks` for the closing run', () => {
    const s = sectionOf(twoAssemblies(), { maxAssemblyBlocks: 3 });
    expect([s.blocks!.length, s.more]).toEqual([3, 1]);
  });
});

describe('the built-in fields’ domains (EC-218, EC-191)', () => {
  it('`block` the blocks’ labels in order; `system` an assembly’s names in its order before any `include`', () => {
    const b = createSeqSpec();
    b.data.assembly('g', assembly('g', grape));
    b.views.view('v', (v) => genome(v, 'g', '/data/g', { include: ['chr2'] }));
    const s = sectionOf(b.build());
    expect(s.domains).toEqual({ block: ['chr2'], system: ['chr1', 'chr2'] });
    expect(s.blocks![0].name).toBe('chr2');
  });
});

describe('places (`placesOf`) and `{ system, position }`', () => {
  const mixed = () =>
    sectionOf(
      piecesSpec([
        [[0, 100], 'section', 1],
        [[100, 250], 'reversed', 3],
        [[50, 120], 'section', 0.5],
        [[400, 500], 'reversed', 1],
      ]),
    );

  it('finds every block showing a position (overlapping blocks: each), and none where no block shows it', () => {
    const s = mixed();
    const sys = '/data/s';
    expect(placesOf(s, sys, 60).map((p) => p.block)).toEqual([0, 2]);
    expect(placesOf(s, sys, 110).map((p) => p.block)).toEqual([1, 2]);
    expect(placesOf(s, sys, 300)).toEqual([]);
    expect(placesOf(s, '/data/none', 10)).toEqual([]);
    // a reversed, scaled block: position 100 at its far end, 3 units a position
    const [p] = placesOf(s, sys, 249);
    expect(p).toEqual({ block: 1, unit: 100 + 1 * 3 });
  });

  it('a position’s unit and back round-trip on every block, whatever its direction and scale (property)', () => {
    const s = mixed();
    const rng = new PCG(7);
    for (let i = 0; i < 400; i++) {
      const b = s.blocks![rng.int() % s.blocks!.length];
      const r = b.axis.ranges!;
      const position = r[0] + ((rng.int() % 10_000) / 10_000) * (r[1] - r[0]);
      const unit = blockUnitOf(b, position)!;
      expect(blockPositionAt(b, unit)).toBeCloseTo(position, 9);
      const at = placeOfUnit(s, unit)!;
      expect(at.system).toBe('/data/s');
      expect(at.position).toBeCloseTo(position, 9);
      // the unit of a place: the block nearest the unit it came from (overlapping blocks show it twice)
      expect(unitOfPlace(s, at, unit)).toBeCloseTo(unit, 9);
    }
  });

  it('a unit on a block’s edge reads the next block’s first as a start and the one before’s last as an end', () => {
    const s = mixed();
    expect(placeOfUnit(s, 100, 'start')).toEqual({ system: '/data/s', position: 250 });
    expect(placeOfUnit(s, 100, 'end')).toEqual({ system: '/data/s', position: 100 });
  });

  it('a position no block shows goes to the nearest shown edge on its system', () => {
    const s = mixed();
    // 300 lies between [100, 250) and [400, 500): its nearest edge is 250, the reversed block's start (unit 100)
    expect(unitOfPlace(s, { system: '/data/s', position: 300 })).toBe(100);
  });

  it('a section without blocks is one place on its own axis', () => {
    const b = createSeqSpec();
    const sq = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
    b.views.view('v', (v) => v.section('a', { of: sq, direction: 'reversed' }));
    const s = sectionOf(b.build());
    expect(placesOf(s, '/data/s', 30)).toEqual([{ block: -1, unit: 70 }]);
  });
});

describe('equal block lists', () => {
  it('two sections of one assembly are equal; a reversed block or a scale makes them unequal', () => {
    const b = createSeqSpec();
    b.data.assembly('g', assembly('g', grape));
    b.views.view('v', (v) => {
      genome(v, 'a', '/data/g');
      genome(v, 'b', '/data/g', {}, (s) =>
        s.block('r', { of: '/data/g/(sequence :name chr1)', direction: 'reversed' } as never),
      );
      genome(v, 'c', '/data/g', { layout: 'rings' });
    });
    const [a, r, c] = resolveView(resolveSpec(b.build()), '/views/v')!.sections;
    expect(sameBlocks(a, c)).toBe(true);
    expect(sameBlocks(a, r)).toBe(false);
  });
});

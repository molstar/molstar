/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for NG1's hot paths, on RB0b's stepped `heapScaling` (a step is one hover change, 50,000 and 200,000 of
// them): a host's loci through `Extents.addLoci`, `project` per block into axis runs (sorted and merged), with (a) a
// direct block, (b) a residue onto a gene block through the codon route, (c) a link to a second section's manager, (d)
// the `containing` closure first; (e) the pointer's `hoverAt` across positions; (f) `rangesOf` from a feature's
// fragments and axis runs as block runs.
// Counted, not timed: a repeated highlight builds no string, and a unit's lanes probe O(n log n). The positive control
// is the template's: (a)'s step with a leak, caught by the same measurement (NG-11).

import { projectionOf, type ResolvedView } from '@molstar/seq-core/format';
import { Extents, type Loci, OrderedSet, type SeqSet, type UnitId } from '@molstar/seq-core/model';
import { beforeAll, describe, expect, it } from 'vitest';
import { constructs, gene, open } from '../../core/_test/format/project-fixtures.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import {
  AxisRuns,
  AxisSink,
  axesOf,
  BlockRuns,
  extentOf,
  highlightRuns,
  type Located,
  prepareLoci,
} from '../src/axis.js';
import { HoverResult, type HoverTarget, hoverAt, newLight, nextHighlight, rangesOf } from '../src/events.js';
import type { AxisMeta } from '../src/index.js';
import { axisView } from '../src/labels.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 50_000;
type Opened = Awaited<ReturnType<typeof open>>;
let g: Opened;
let managers: AxisMeta[];

const seq = (set: SeqSet, unit: number, i: number): Loci => ({
  kind: 'seq-loci',
  set,
  elements: [{ unit: unit as UnitId, indices: OrderedSet.ofSingleton(i) as never }],
});
/** One prebuilt loci per position, so a step allocates none. */
const prebuilt = (set: SeqSet, unit: number, n: number) => Array.from({ length: n }, (_, i) => seq(set, unit, i));
const axes = (view: ResolvedView, set: SeqSet) =>
  Array.from({ length: view.sections.count }, (_, s) => axesOf(view, set, s)).flat();

beforeAll(async () => {
  g = await open(gene({ chromosome: true }));
  managers = axes(g.view, g.set);
}, 60_000);

/** The unit id behind a coordinate system key. */
const unitOf = (o: Opened, key: string) => o.system(key).id;

describe('the highlight path allocates nothing per change', () => {
  /** A highlight change's step, and what it lit so far. */
  const stepOf = (o: Opened, loci: Loci[], closure: 'element' | 'containing', ms: AxisMeta[]) => {
    const [p, ext, sink, out] = [projectionOf(o.view, o.set), Extents.create(8), new AxisSink(), AxisRuns.create()];
    const lit = { n: 0 };
    const step = (i: number) => {
      if (!prepareLoci(p, ext, loci[i % loci.length], closure)) return;
      for (let m = 0; m < ms.length; m++) lit.n += highlightRuns(p, ext, ms[m], sink, out).count;
    };
    return { step, lit };
  };
  const run = (o: Opened, loci: Loci[], closure: 'element' | 'containing', ms: AxisMeta[]) => {
    const { step, lit } = stepOf(o, loci, closure, ms);
    const result = heapScaling(N, step);
    expect(lit.n).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  };
  const gene120 = () => prebuilt(g.set, unitOf(g, 'gene/unit'), 120);

  it('(a) a nucleotide on its direct block', () => run(g, gene120(), 'element', [managers[0]]));

  it('(a) with a leak is caught (the positive control)', () => {
    const { step } = stepOf(g, gene120(), 'element', [managers[0]]);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it('(b) a residue on the gene block through the codon route', () => {
    run(g, prebuilt(g.set, unitOf(g, 'protein'), 20), 'element', [managers[0]]);
  });

  it('(c) a residue to every manager: the protein, the gene and the chromosome across the link', () => {
    run(g, prebuilt(g.set, unitOf(g, 'protein'), 20), 'element', managers);
  });

  it('(d) a nucleotide with the containing closure first, then projected', () => {
    run(g, gene120(), 'containing', managers);
  });

  it('a source used twice, at the end of its other section', async () => {
    const c = await open(constructs());
    run(c, prebuilt(c.set, unitOf(c, 'ub'), 76), 'element', axes(c.view, c.set));
  });

  it("(e) the pointer's hoverAt across positions, deduped before it locates", () => {
    const meta = managers[0];
    const view = axisView({ start: 0, end: 0, left: 0, right: 0 }, 1, -1, meta.length, 10 * meta.length + 20, 10, 10);
    const [state, out, layer, manager]: [HoverTarget, Located, object, object] = [
      { manager: null, position: -1, layer: null },
      { block: -1, index: -1 },
      {},
      {},
    ];
    let moved = 0;
    const result = heapScaling(N, (i) => {
      // Three pointer moves in a position: only the first is a change.
      if (hoverAt(state, manager, meta, view, 10 + 10 * (i % meta.length) + (i % 3), layer, out) === HoverResult.moved)
        moved++;
    });
    expect(moved).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it("(f) a feature's fragments as axis runs, and axis runs as block runs", () => {
    const meta = managers[0];
    const detail = {
      feature: {
        svs: { layer: 0 },
        locations: [
          {
            fragments: [
              { start: 3, end: 10 },
              { start: 30, end: 60 },
            ],
          },
        ],
      },
    };
    const [axis, blocks] = [AxisRuns.create(), BlockRuns.create()];
    let seen = 0;
    const result = heapScaling(N, (i) => {
      detail.feature.locations[0].fragments[0].start = 3 + (i % 5);
      seen += extentOf(meta, rangesOf(detail, 3, axis), blocks).count;
    });
    expect(seen).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it('a repeated highlight builds no string and dispatches nothing', () => {
    const [p, ext, sink, out, light] = [
      projectionOf(g.view, g.set),
      Extents.create(8),
      new AxisSink(),
      AxisRuns.create(),
      newLight(),
    ];
    const loci = prebuilt(g.set, unitOf(g, 'protein'), 20)[4];
    const sent: unknown[] = [];
    for (let i = 0; i < 100; i++) {
      prepareLoci(p, ext, loci, 'element');
      const next = nextHighlight(light, highlightRuns(p, ext, managers[0], sink, out));
      if (next !== undefined) sent.push(next);
    }
    expect(sent).toHaveLength(1);
    expect(typeof sent[0]).toBe('string');
  });
});

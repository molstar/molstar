/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the highlight path's allocation guards (NG1; hover.heap.ts) and their quick twin (hover.quick.ts): each
// builds a step, one hover change, over the gene and constructs fixtures, and what the step lit so far. (a) a
// nucleotide on its direct block, (b) a residue on the gene block through the codon route, (c) a residue to every
// manager across the link, (d) a nucleotide with the `containing` closure first, a source used twice at the end of
// its other section, (e) the pointer's `hoverAt` across positions, (f) `rangesOf` from a feature's fragments and
// axis runs as block runs. The heap guard and the quick twin share them: neither copies a setup.

import { projectionOf, type ResolvedView } from '@molstar/seq-core/format';
import { Extents, type Loci, OrderedSet, type SeqSet, type UnitId } from '@molstar/seq-core/model';
import { constructs, gene, open } from '../../core/_test/format/project-fixtures.js';
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

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;

export type Opened = Awaited<ReturnType<typeof open>>;

/** A case: its step, and what the steps found so far (zero means the step did nothing). */
export interface Case {
  step: (i: number) => void;
  seen: () => number;
}

const seq = (set: SeqSet, unit: number, i: number): Loci => ({
  kind: 'seq-loci',
  set,
  elements: [{ unit: unit as UnitId, indices: OrderedSet.ofSingleton(i) as never }],
});
/** One prebuilt loci per position, so a step allocates none. */
export const prebuilt = (set: SeqSet, unit: number, n: number) =>
  Array.from({ length: n }, (_, i) => seq(set, unit, i));
const axes = (view: ResolvedView, set: SeqSet) =>
  Array.from({ length: view.sections.count }, (_, s) => axesOf(view, set, s)).flat();
/** The unit id behind a coordinate system key. */
export const unitOf = (o: Opened, key: string) => o.system(key).id;

/** The gene fixture (with its chromosome) opened, and the managers of its sections. */
export async function genes() {
  const g = await open(gene({ chromosome: true }));
  return { g, managers: axes(g.view, g.set) };
}

/** A highlight change's step, and what it lit so far. */
function highlighting(o: Opened, loci: Loci[], closure: 'element' | 'containing', ms: AxisMeta[]): Case {
  const [p, ext, sink, out] = [projectionOf(o.view, o.set), Extents.create(8), new AxisSink(), AxisRuns.create()];
  let lit = 0;
  const step = (i: number) => {
    if (!prepareLoci(p, ext, loci[i % loci.length], closure)) return;
    for (let m = 0; m < ms.length; m++) lit += highlightRuns(p, ext, ms[m], sink, out).count;
  };
  return { step, seen: () => lit };
}

/** (a) a nucleotide on its direct block. */
export async function direct(): Promise<Case> {
  const { g, managers } = await genes();
  return highlighting(g, prebuilt(g.set, unitOf(g, 'gene/unit'), 120), 'element', [managers[0]]);
}

/** (b) a residue on the gene block through the codon route. */
export async function codon(): Promise<Case> {
  const { g, managers } = await genes();
  return highlighting(g, prebuilt(g.set, unitOf(g, 'protein'), 20), 'element', [managers[0]]);
}

/** (c) a residue to every manager: the protein, the gene and the chromosome across the link. */
export async function linked(): Promise<Case> {
  const { g, managers } = await genes();
  return highlighting(g, prebuilt(g.set, unitOf(g, 'protein'), 20), 'element', managers);
}

/** (d) a nucleotide with the containing closure first, then projected. */
export async function containing(): Promise<Case> {
  const { g, managers } = await genes();
  return highlighting(g, prebuilt(g.set, unitOf(g, 'gene/unit'), 120), 'containing', managers);
}

/** A source used twice, at the end of its other section. */
export async function twice(): Promise<Case> {
  const c = await open(constructs());
  return highlighting(c, prebuilt(c.set, unitOf(c, 'ub'), 76), 'element', axes(c.view, c.set));
}

/** (e) the pointer's `hoverAt` across positions, three moves in a position: only the first is a change. */
export async function pointer(): Promise<Case> {
  const { managers } = await genes();
  const meta = managers[0];
  const view = axisView({ start: 0, end: 0, left: 0, right: 0 }, 1, -1, meta.length, 10 * meta.length + 20, 10, 10);
  const [state, out, layer, manager]: [HoverTarget, Located, object, object] = [
    { manager: null, position: -1, layer: null },
    { block: -1, index: -1 },
    {},
    {},
  ];
  let moved = 0;
  const step = (i: number) => {
    if (hoverAt(state, manager, meta, view, 10 + 10 * (i % meta.length) + (i % 3), layer, out) === HoverResult.moved)
      moved++;
  };
  return { step, seen: () => moved };
}

/** (f) a feature's fragments as axis runs, and axis runs as block runs. */
export async function fragments(): Promise<Case> {
  const { managers } = await genes();
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
  const step = (i: number) => {
    detail.feature.locations[0].fragments[0].start = 3 + (i % 5);
    seen += extentOf(meta, rangesOf(detail, 3, axis), blocks).count;
  };
  return { step, seen: () => seen };
}

/** A repeated highlight builds no string and dispatches nothing: the strings sent for a hundred identical highlights. */
export async function repeated(): Promise<unknown[]> {
  const { g, managers } = await genes();
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
  return sent;
}

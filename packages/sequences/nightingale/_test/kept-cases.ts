/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of the kept highlights' allocation guard (LS1; kept.heap.ts) and its quick twin (kept.quick.ts): a section
// of the gene fixture with forty highlight runs kept (`AxisMeta.highlights`), and a step that is a host's highlight
// change projected, then kept under them (`keepLit`). The heap guard and the quick twin share it.

import { projectionOf } from '@molstar/seq-core/format';
import { Extents, type Loci, OrderedSet, type UnitId } from '@molstar/seq-core/model';
import { gene, open } from '../../core/_test/format/project-fixtures.js';
import { AxisRuns, AxisSink, axesOf, highlightRuns, keepLit, prepareLoci } from '../src/axis.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;

/** The step, and the runs kept so far (the guard's non-vacuity: forty or more a step). */
export async function kept(): Promise<{ step: (i: number) => void; lit: () => number }> {
  const g = await open(gene({ chromosome: true }));
  const runs = AxisRuns.create();
  for (let k = 0; k < 40; k++) AxisRuns.add(runs, 7 * k, 7 * k + 3);
  const meta = { ...axesOf(g.view, g.set, 0)[0], highlights: AxisRuns.finish(runs) };
  const unit = g.system('gene/unit').id as UnitId;
  const loci: Loci[] = Array.from({ length: 120 }, (_, i) => ({
    kind: 'seq-loci',
    set: g.set,
    elements: [{ unit, indices: OrderedSet.ofSingleton(i) as never }],
  }));
  const [p, ext, sink, out] = [projectionOf(g.view, g.set), Extents.create(8), new AxisSink(), AxisRuns.create()];
  let lit = 0;
  const step = (i: number) => {
    if (prepareLoci(p, ext, loci[i % loci.length], 'element')) highlightRuns(p, ext, meta, sink, out);
    else out.count = 0;
    lit += keepLit(out, meta.highlights).count;
  };
  return { step, lit: () => lit };
}

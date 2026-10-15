/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the retention check (IN1; retention.heap.ts) and its quick twin (retention.quick.ts). The check itself,
// weak references to a disposed viewer's set and view collected after a task apart, is a once-only, async measurement
// with no step: the heap guard alone has it (`used`, `collected`). Its repeatable synchronous core is the projection
// the viewer is used through, `project`: an item's extents, closed over the domain, projected into a section's hover
// marks, then closed again; `projecting` is that as a step on a viewer, and the quick twin measures it. The heap guard
// and the quick twin share `project`.

import type { Loci } from '@molstar/seq-core/model';
import { expect } from 'vitest';
import { case14 } from '../../core/_test/format/cases.js';
import { SectionMarks } from '../src/marks.js';
import { shown } from './projection-fixtures.js';

type Shown = Awaited<ReturnType<typeof shown>>;
const gc = (globalThis as { gc?: () => void }).gc as () => void;

/** A viewer on case 14, a loci of its `domain` item, and a section's marks. */
async function open() {
  const t = await shown(case14());
  return { t, domain: t.ctx.refs.toLoci('domain') as Loci, marks: new SectionMarks() };
}

/** The projection: the domain's extents closed over its container, projected into the hover marks, closed again. */
function project(t: Shown, domain: Loci, marks: SectionMarks) {
  t.p.project(t.extents(domain, 'containing'), 0, marks.hover);
  t.p.closeLoci(domain, 'containing');
}

/** The projection as a step on an open viewer; `projected` reads what the marks hold (not vacuous when above zero). */
export async function projecting(): Promise<{
  step: (i: number) => void;
  projected: () => number;
  dispose: () => void;
}> {
  const { t, domain, marks } = await open();
  return { step: () => project(t, domain, marks), projected: () => marks.hover.count, dispose: () => t.ctx.dispose() };
}

/** A viewer used through its projectors, disposed; weak references to its set and view. `kept`: the host keeps the set. */
export async function used(kept?: unknown[]) {
  const { t, domain, marks } = await open();
  expect(t.ctx.focus.fromLoci(domain, { transitionMs: 0 })).toBe(true);
  project(t, domain, marks);
  const refs = { set: new WeakRef(t.set), view: new WeakRef(t.view) };
  t.ctx.dispose();
  kept?.push(t.set);
  return refs;
}

/** How many of the set and the view are still reachable after collecting, a task apart, up to four times. */
export async function collected(refs: { set: WeakRef<object>; view: WeakRef<object> }) {
  for (let k = 0; k < 4 && refs.set.deref(); k++) {
    await new Promise((r) => setTimeout(r, 0));
    gc();
  }
  return (refs.set.deref() ? 1 : 0) + (refs.view.deref() ? 1 : 0);
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the projector and marks tests (IN1): a context on the example fixtures' files with a view shown, and a
// projector over it reading what a target lights per block (by label), as runs or as a count of positions.

import type { ResolvedView, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { Extents, type Loci, type SeqSet } from '@molstar/seq-core/model';
import { fixtureLoaders, MolstarFixture } from '../../core/_test/format/load-fixtures.js';
import { type Coordinates, createSeqViewer, DefaultSeqViewerSpec, manualScheduler } from '../src/index.js';
import { MarkSet } from '../src/marks.js';
import { Projector } from '../src/projection.js';
import { kernelFree, stubs } from './graph-fixtures.js';
/** A context with the structure extension's fixture showing `doc`'s view `view` (else its first), and a projector. */
export async function shown(doc: SeqViewSpecDoc, view?: string) {
  const scheduler = manualScheduler();
  const registry = [kernelFree(), stubs().entry, { svs: MolstarFixture }];
  const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), registry }, { scheduler, loaders: fixtureLoaders() });
  const v = { ctx };
  const r = await v.ctx.setDocument(doc);
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  if (view) {
    v.ctx.views.setCurrent(view);
    await v.ctx.graph.update();
  }
  const a = v.ctx.graph.assembled.value;
  const set = a?.set as SeqSet;
  const rv = a?.view?.resolved as ResolvedView;
  return { ...v, scheduler, ...over(rv, set) };
}

/** A projector over a resolved view, with readers. */
export function over(view: ResolvedView, set: SeqSet) {
  const p = new Projector();
  p.update(view, set);
  const ext = Extents.create(4);
  const marks = new MarkSet();
  const label = (b: number) => set.strings.get(view.blocks.label[b]);
  const unit = (key: string) =>
    set.units.key.findIndex((k, u) => set.units.alive[u] === 1 && set.strings.get(k) === key);
  /** The target as extents (closed when asked). */
  const extents = (target: Loci | Coordinates | Extents, closure: 'element' | 'containing' = 'element') => {
    if (target === ext) return ext;
    ext.count = 0;
    if ('kind' in target && typeof target.kind === 'string') p.fromLoci(target, ext);
    else if ('runs' in target) p.fromLoci(target, ext);
    else for (let i = 0; i < target.count; i++) Extents.add(ext, at(target, i), target.start[i], target.stop[i]);
    if (closure === 'containing') p.close(ext);
    return ext;
  };
  /** Per block label (every section): its lit runs, merged and clipped as a MarkSet holds them. */
  const lit = (target: Loci | Coordinates | Extents, closure: 'element' | 'containing' = 'element') => {
    const e = extents(target, closure);
    const out: Record<string, number[][]> = {};
    for (let s = 0; s < view.sections.count; s++) {
      p.project(e, s, marks);
      for (let b = 0; b < view.sections.blockCount[s]; b++) {
        const runs: number[][] = [];
        for (let r = marks.from(b); r < marks.to(b); r++) runs.push([marks.start[r], marks.stop[r]]);
        out[label(view.sections.blockStart[s] + b)] = runs;
      }
    }
    return out;
  };
  /** Per block label: how many positions the target lights. */
  const marked = (target: Loci | Coordinates | Extents, closure: 'element' | 'containing' = 'element') =>
    Object.fromEntries(
      Object.entries(lit(target, closure)).map(([k, runs]) => [k, runs.reduce((n, [s, t]) => n + t - s, 0)]),
    );
  /** The positions a target lights on one block, one by one. */
  const positions = (block: string, target: Loci | Coordinates | Extents, closure?: 'element' | 'containing') =>
    (lit(target, closure)[block] ?? []).flatMap(([s, t]) => Array.from({ length: t - s }, (_, k) => s + k));
  return { set, view, p, ext, marks, label, unit, extents, lit, marked, positions };
}

const sys = { kind: 0, id: 0, row: -1 };
function at(e: Extents, i: number) {
  sys.kind = e.kind[i];
  sys.id = e.id[i];
  sys.row = e.row[i];
  return sys;
}

/** Extents of one entry: `[kind, id, row, start, stop]`. */
export function entries(...list: [number, number, number, number, number][]) {
  const e = Extents.create(1);
  for (const [kind, id, row, s, t] of list) Extents.add(e, { kind, id, row }, s, t);
  return e;
}

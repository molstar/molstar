/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's colors and an annotation's rows on units (_spec/README.md "Structures", Colors; Rows): letters, values and
// ranges layers compiled as drawn, light and dark; unobserved positions cut (EC-118), uncolored ones in no run, every
// copy colored (EC-65), data the store does not hold pending and never an empty run (EC-169).

import { DataStore } from '@molstar/seqx-core/data';
import { partOf, positionCount } from '@molstar/seqx-core/query';
import { openModel } from '@molstar/seqx-core/render';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { colorRuns } from '../colors.js';
import { unitRows } from '../rows.js';
import { chainId, Complex, EntStart, Observed, p53Spec, Up, UpLength } from './p53.js';

const unit = (k: number) => `${Complex}/${chainId(k)}`;

async function opened(field = false, settle = true) {
  const store = new DataStore();
  const spec = p53Spec({ view: { field } });
  const { model, view } = await openModel(spec, store);
  const layer = (track: number) => view!.sections[0].tracks[track].layers[0];
  if (settle)
    for (const t of [1, 2]) {
      colorRuns(model, store, layer(t), undefined, { want: (w) => store.want([w]) });
      await store.whenSettled();
    }
  return { model, store, view: view!, layer };
}

describe('colorRuns', () => {
  it('letters: a run per color, on every copy, unobserved positions cut, light and dark (EC-65, EC-118)', async () => {
    const { model, store, layer } = await opened();
    const r = colorRuns(model, store, layer(0));
    expect([r.complete, r.issues]).toEqual([true, []]);
    expect(r.runs.length).toBeGreaterThan(3);
    let lit = 0;
    for (const run of r.runs) {
      expect(typeof run.color === 'number' && typeof run.dark === 'number').toBe(true);
      expect(run.target.parts.map((p) => p.system)).toEqual([0, 1, 2, 3].map(unit));
      lit += positionCount(run.target);
    }
    // the observed part of the entity on each of the four units, every residue colored once
    expect(lit).toBe(4 * (Observed[1] - Observed[0]));
    expect(r.runs.some((x) => x.color !== x.dark)).toBe(true);
  });

  it('values: an open scale over the data’s extent; a missing value in no run', async () => {
    const { model, store, layer } = await opened();
    const r = colorRuns(model, store, layer(2));
    expect(r.complete).toBe(true);
    // position 100 of P04637 is NaN: no run lights its residue on a unit; its neighbors are lit
    const covered = (p: number) =>
      r.runs.some((x) => {
        const g = partOf(x.target, unit(0))?.ranges ?? [];
        for (let i = 0; i + 1 < g.length; i += 2) if (g[i] <= p - EntStart && p - EntStart < g[i + 1]) return true;
        return false;
      });
    expect([covered(99), covered(100), covered(101)]).toEqual([true, false, true]);
    expect(new Set(r.runs.map((x) => x.color)).size).toBe(r.runs.length);
  });

  it('ranges: a uniform color one run; a field’s colors one run each value', async () => {
    const uniform = await opened(false);
    const r = colorRuns(uniform.model, uniform.store, uniform.layer(1));
    expect(r.runs.map((x) => x.color)).toEqual([0xff0000ff]);
    // DBD 99..120 on P04637: the entity's 6..27 on each unit
    expect(Array.from(partOf(r.runs[0].target, unit(3))!.ranges)).toEqual([99 - EntStart, 120 - EntStart]);
    const field = await opened(true);
    const f = colorRuns(field.model, field.store, field.layer(1));
    // only `binding` reaches the structure (DBD; TET is past the entity)
    expect(f.runs.length).toBe(1);
  });

  it('data the store does not hold is pending, never an empty run (EC-169)', async () => {
    const { model, store, layer } = await opened(false, false);
    const r = colorRuns(model, store, layer(2));
    expect(r.complete).toBe(false);
    expect(Array.from(partOf(r.pending, Up)!.ranges)).toEqual([0, UpLength]);
    expect(r.runs).toEqual([]);
  });

  it('a layer with no colors to give is said', async () => {
    const { model, store, view } = await opened();
    const fake = { ...view.sections[0].tracks[0].layers[0], type: 'ruler', fallbacks: [] };
    expect(colorRuns(model, store, fake).issues.map((i) => i.code)).toEqual(['seq.structure.layer']);
    expect(resolveViews(model).length).toBe(1);
  });
});

describe('unitRows', () => {
  it('each item reaching a unit, with its fields; items past the structure left out', async () => {
    const { model, store } = await opened();
    const r = unitRows(model, store, '/data/domains');
    expect(r.complete).toBe(true);
    expect(r.rows.map((x) => [x.key, x.fields.category])).toEqual([['DBD', 'binding']]);
    expect(r.rows[0].target.parts.map((p) => p.system)).toEqual([0, 1, 2, 3].map(unit));
    expect(unitRows(model, store, '/data/scores').issues.map((i) => i.code)).toEqual(['seq.structure.rows']);
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b's heap guards (seqspec:sources#the-viewers-side; the memory rule): a pan's pass and composer frame over a
// template's 1,000 instance tracks allocate nothing once laid out (the instances are tracks like any other); a frame
// that allocates an object a step fails the same guard (the positive control), warmed and settled by count
// (`frameScaling`); and an instance that leaves (its node evicted) is not retained: its resolved layer and prepared
// state are freed once a later update has run (a WeakRef check, so it has no quick twin). The pan is
// templates-cases.ts's, shared with the quick twin (templates.quick.ts).

import { describe, expect, it, vi } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { INSTANCES, templated } from './templates-cases.js';

const N = 1_000;
const gc = (globalThis as { gc?: () => void }).gc as () => void;

describe("a template's 1,000 instance tracks", () => {
  it("a pan's pass and frame allocate nothing", async () => {
    const g = await templated();
    expect(g.t.ctx.graph.assembled.value?.view?.resolved.tracks.count).toBeGreaterThan(INSTANCES);
    const r = frameScaling(g.frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await templated();
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
  }, 300_000);
});

describe('an instance that leaves', () => {
  it('is not retained once a later update has run', async () => {
    vi.resetModules();
    const { openTemplates, range } = await import('./template-fixtures.js');
    const t = await openTemplates({ genes: [range('a', 1)] });
    const find = (node: string) => {
      const v = t.ctx.graph.assembled.value?.view;
      for (let k = 0; v && k < v.resolved.layers.count; k++)
        if (v.resolved.layers.node[k] === node)
          return { state: v.prepared[k].state as object, themes: v.prepared[k].themes };
      return null;
    };
    const refs = await (async () => {
      t.manager.nodes.add('genes', [range('b', 9)]);
      await t.frame();
      const b = find('genes/b');
      if (!b) throw new Error('no instance');
      return { state: new WeakRef(b.state), themes: new WeakRef(b.themes) };
    })();
    t.manager.nodes.release(['genes/b']);
    await t.frame();
    t.manager.nodes.add('genes', [range('c', 12)]);
    await t.frame();
    expect(find('genes/b')).toBeNull();
    for (let k = 0; k < 4 && refs.state.deref(); k++) {
      await new Promise((r) => setTimeout(r, 0));
      gc();
    }
    expect([refs.state.deref(), refs.themes.deref()]).toEqual([undefined, undefined]);
    t.ctx.dispose();
  }, 60_000);
});

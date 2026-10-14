/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The DOM parts' subscriptions: a state replays its value on subscribe, which is not a change, so a part (or `derived`,
// `watch`) built over several states runs its work once at mount and once per change, never once per state's replay.

import { afterEach, describe, expect, it, vi } from 'vitest';
import { derived, later, watch } from '../../src/dom/core.js';
import * as D from '../../src/dom/index.js';
import { SeqViewerConfig, type SeqViewerContext, shownSections } from '../../src/index.js';
import { sameData, sameShallow } from '../../src/same.js';
import { SearchBox } from '../../src/search-box.js';
import { live, state } from '../../src/state.js';
import { at, flush, openedContext, pointer, stubSurface, uiDoc, WIDTH } from './fixtures.js';

const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
  vi.restoreAllMocks();
});

describe('later, watch and derived', () => {
  it('`later` skips the replay and calls back on each change', () => {
    const a = state(1);
    const got: number[] = [];
    const sub = later(a, (v) => got.push(v));
    expect(got).toEqual([]);
    a.set(2);
    a.set(3);
    expect(got).toEqual([2, 3]);
    sub.unsubscribe();
    a.set(4);
    expect(got).toEqual([2, 3]);
  });

  it('`watch` runs once at mount over any number of states, then once per change', () => {
    const [a, b] = [state(1), state('x')];
    const l = live({ n: 0 });
    let runs = 0;
    const subs = watch([a, b, l.version], () => runs++);
    expect(runs).toBe(1);
    a.set(2);
    b.set('y');
    l.bump();
    expect(runs).toBe(4);
    for (const s of subs) s.unsubscribe();
    a.set(3);
    expect(runs).toBe(4);
  });

  it('`derived` computes once at mount over several states and once per change; apply only when the value changed', () => {
    const [a, b, c] = [state(3), state(2), state(1)];
    let [computes, applies] = [0, 0];
    const subs = derived(
      [a, b, c],
      () => {
        computes++;
        return a.value + b.value + c.value > 5;
      },
      () => applies++,
    );
    expect([computes, applies]).toEqual([1, 1]);
    // Another sum, the same answer: computed, not applied.
    a.set(4);
    expect([computes, applies]).toEqual([2, 1]);
    b.set(0);
    c.set(0);
    expect([computes, applies]).toEqual([4, 2]);
    for (const s of subs) s.unsubscribe();
  });

  it('`derived` sees a changed handler and an undefined that becomes null, which a JSON key does not', () => {
    const a = state(0);
    const [f, g] = [() => 1, () => 2];
    let value: { run: () => number; list: (number | undefined)[] } = { run: f, list: [1, undefined] };
    const seen: unknown[] = [];
    derived(
      [a],
      () => value,
      (v) => seen.push(v),
      sameData,
    );
    expect(seen.length).toBe(1);
    // The same content again: a new object of the same members is no change.
    value = { run: f, list: [1, undefined] };
    a.set(1);
    expect(seen.length).toBe(1);
    // Only the handler differs: JSON.stringify omits functions.
    value = { run: g, list: [1, undefined] };
    a.set(2);
    expect(seen.length).toBe(2);
    // `undefined` to `null` in an array: JSON writes both as null.
    value = { run: g, list: [1, null as never] };
    a.set(3);
    expect(seen.length).toBe(3);
    // A list of records is compared by `sameData`.
    let items = [{ ref: 'a', label: 'A' }];
    const lists: unknown[] = [];
    derived(
      [a],
      () => items.map((i) => ({ ...i })),
      (v) => lists.push(v),
      sameData,
    );
    items = [{ ref: 'a', label: 'A' }];
    a.set(4);
    expect(lists.length).toBe(1);
    items = [{ ref: 'a', label: 'B' }];
    a.set(5);
    expect(lists.length).toBe(2);
  });
});

describe('sameShallow and sameData', () => {
  const f = () => {};
  it('sameShallow: one level, by identity', () => {
    expect(sameShallow(undefined, undefined)).toBe(true);
    expect(sameShallow(undefined, null)).toBe(false);
    expect(sameShallow([1, 'a'], [1, 'a'])).toBe(true);
    expect(sameShallow([1, undefined], [1, null])).toBe(false);
    expect(sameShallow({ a: 1, b: 2 }, { b: 2, a: 1 })).toBe(true);
    expect(sameShallow({ a: 1, b: undefined }, { a: 1 })).toBe(false);
    expect(sameShallow({ run: f }, { run: f })).toBe(true);
    expect(sameShallow({ run: f }, { run: () => {} })).toBe(false);
    expect(sameShallow([{ a: 1 }], [{ a: 1 }])).toBe(false);
    expect(sameShallow(new Set([1]), new Set([1]))).toBe(false);
  });
  it('sameData: every level of plain arrays and records', () => {
    expect(sameData([{ a: [1, { b: 2 }] }], [{ a: [1, { b: 2 }] }])).toBe(true);
    expect(sameData([{ a: [1, { b: 2 }] }], [{ a: [1, { b: 3 }] }])).toBe(false);
    expect(sameData({ run: f }, { run: () => {} })).toBe(false);
    expect(sameData({ run: f }, { run: f })).toBe(true);
    expect(sameData([1], [1, 2])).toBe(false);
  });
});

async function mountedPart<P extends object>(Part: D.DomPart<P>, props: (ctx: SeqViewerContext) => P) {
  const { ctx } = await openedContext();
  const host = document.createElement('div');
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  const ui = D.createDomUI(ctx, { surface: stubSurface });
  cleanups.push(() => {
    host.remove();
    ctx.dispose();
  });
  return {
    ctx,
    ui,
    host,
    mount: () => {
      const p = Part(host, ctx, { ui, ...props(ctx) } as never);
      cleanups.unshift(() => p.dispose());
      return p;
    },
  };
}

describe('a part over several states', () => {
  it('the search box writes once at mount (its view and the shown document) and once per change', async () => {
    const t = await mountedPart(D.Search as D.DomPart<object>, () => ({ search: async () => ({ results: [] }) }));
    const status = vi.spyOn(SearchBox.prototype, 'status');
    t.mount();
    expect(status).toHaveBeenCalledTimes(1);
    status.mockClear();
    // The document the part shows changes: one write.
    await t.ctx.setDocument(uiDoc());
    await flush();
    expect(status).toHaveBeenCalledTimes(1);
  });

  it('the track menu is built once at mount (shown, expanded, assembled) and once per change', async () => {
    const t = await mountedPart(D.TrackActions as D.DomPart<object>, () => ({
      section: 'chains.main',
      track: 't.dom',
    }));
    const can = vi.spyOn(t.ctx.actions, 'can');
    t.mount();
    // One build of the menu asks `can` for each of its actions.
    const build = can.mock.calls.length;
    expect(build).toBeGreaterThan(0);
    const expanded = t.ctx.ui.expanded as unknown as { set(v: ReadonlySet<string>, o: string): void };
    expanded.set(new Set(['t.dom']), 'api');
    expect(can.mock.calls.length).toBe(2 * build);
    expanded.set(new Set(), 'api');
    expect(can.mock.calls.length).toBe(3 * build);
  });

  it("the track menu is written again when a host's handler changes, though no label does", async () => {
    const { ctx } = await openedContext();
    const ran: string[] = [];
    const entry = (tag: string) => ({ id: 'log', label: 'Log it', run: () => void ran.push(tag) });
    ctx.config.set(SeqViewerConfig.UI.TrackActions, ['hide', entry('one')]);
    const host = document.createElement('div');
    document.body.appendChild(host);
    const m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
    cleanups.unshift(() => {
      m.dispose();
      host.remove();
      ctx.dispose();
    });
    const h = at(ctx, 'chains.main', 't.dom');
    pointer(ctx)('move', 'chains.main', h.header, h.y);
    const menuItems = () => {
      (host.querySelector('.seqv-track-actions button') as HTMLButtonElement).click();
      return [...host.querySelectorAll<HTMLButtonElement>('[role=menuitem]')];
    };
    menuItems()[1].click();
    expect(ran).toEqual(['one']);
    // Same ids and labels, another handler: a JSON key would keep the first.
    ctx.config.set(SeqViewerConfig.UI.TrackActions, ['hide', entry('two')]);
    await flush();
    menuItems()[1].click();
    expect(ran).toEqual(['one', 'two']);
  });
});

describe("a section's chrome", () => {
  it('is computed once per change for the header, the switches and the chips together (DOM)', async () => {
    const { ctx } = await openedContext();
    const sections = shownSections(ctx).length;
    expect(sections).toBeGreaterThan(0);
    // `sectionInfo` asks the registry for the presets once per section it computes.
    const list = vi.spyOn(ctx.registry, 'list');
    const computes = () => list.mock.calls.filter((c) => c[0] === 'presets').length;
    const host = document.createElement('div');
    host.style.width = `${WIDTH}px`;
    document.body.appendChild(host);
    const m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
    cleanups.unshift(() => {
      m.dispose();
      host.remove();
      ctx.dispose();
    });
    await flush();
    expect(host.querySelector('.seqv-section-header')).not.toBeNull();
    expect(host.querySelector('.seqv-switches')).not.toBeNull();
    expect(computes()).toBe(sections);
    expect(ctx.actions.run('section.style', { section: 'chains.main', value: 'railed' })).toBe(true);
    await flush();
    expect(computes()).toBe(2 * sections);
  });
});

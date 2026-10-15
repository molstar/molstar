/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The track menu (seqspec:interaction#the-track-menu): the browser's context-menu request (`ctx.input.menu`, the
// engine's `contextmenu`) is `interaction.menu`, hit at its point with the track's ref and the place; a subscriber's
// `preventDefault` shows its menu instead of the browser's. From the keyboard (the menu key, Shift+F10, a keyboard
// context-menu request) it is asked at the cursor's cell, origin `keyboard`, once per ask: the browser's own request
// that follows a key's on the same section, with no press, later key, blur or leave between, is not asked again and has
// its default prevented as the key's was.

import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { Hit, HitRegion, Inputs, type SeqInputEvent, type SeqNativeEvent } from '../src/index.js';
import { type Driven, driven } from './input-fixtures.js';

const S = 'family.main';

/** A DOM event's stand-in: records `preventDefault`. */
function native(type = 'contextmenu') {
  const e = {
    type,
    defaultPrevented: false,
    preventDefault() {
      e.defaultPrevented = true;
    },
  };
  return e satisfies SeqNativeEvent;
}

/** case 3 in 1000 × 300: the navigator on top, the pinned alignment rows, one scrolled track, then empty rows. */
async function opened(prevent: (e: SeqInputEvent) => boolean = () => false) {
  const t = await driven(case3(), { size: { [S]: [1000, 300] } });
  const asked: SeqInputEvent[] = [];
  const base: boolean[] = [];
  t.ctx.input.events.menu.subscribe((e) => base.push(e.keyboard));
  t.ctx.interaction.menu.subscribe((e) => {
    asked.push(e);
    if (prevent(e)) e.native.preventDefault();
  });
  const m = Inputs.pointer();
  /** The engine's `contextmenu` at a point; `keyboard`: no pointer made it (`pointerId` -1). */
  const menu = (p: { x: number; y: number }, opts: { keyboard?: boolean; button?: number; ctrl?: boolean } = {}) => {
    const n = native();
    Object.assign(m, {
      section: t.ctx.layout.sectionIndex(S),
      x: p.x,
      y: p.y,
      button: opts.button ?? 2,
      buttons: opts.button === undefined ? 2 : 0,
      modifiers: opts.ctrl ? 2 : 0,
      pointerId: opts.keyboard ? -1 : 1,
      native: n,
    });
    t.ctx.input.menu(m);
    m.native = null;
    return n;
  };
  return { t, asked, base, menu };
}
/** The middle of the first scrolled track's data, in surface px. */
const onTrack = (t: Driven) => {
  const l = t.layout(S);
  const e = l.pinnedCount;
  return { x: l.dataLeft + 40, y: l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 };
};
/** The first point down the column at `x` whose hit is in `region` with (or without) a track. */
const find = (t: Driven, x: number, region: number, track: boolean) => {
  const h = Hit.create();
  for (let y = 1; y < t.layout(S).height; y++) {
    t.ctx.layout.hitTest(S, x, y, h);
    if (h.region === region && h.track >= 0 === track) return { x, y };
  }
  throw new Error(`no ${track ? 'track' : 'empty'} point in region ${region}`);
};

describe("the browser's context-menu request (interaction.menu)", () => {
  it("is hit at its point: the track's ref, the part and the place; the default stays unless a subscriber prevents it", async () => {
    const { t, asked, base, menu } = await opened((e) => e.part.track !== null);
    const p = onTrack(t);
    const n = menu(p, { ctrl: true });
    expect(base).toEqual([false]);
    expect(asked).toHaveLength(1);
    const e = asked[0];
    expect(e).toMatchObject({ stream: 'menu', origin: 'pointer', owner: '', buttons: 2, button: 2 });
    expect(e.part).toMatchObject({ kind: 'data', section: S, key: null });
    expect(e.part.track).toBeTruthy();
    expect(e.place).toMatchObject({ x: p.x, y: p.y });
    expect(e.modifiers).toEqual({ shift: false, control: true, alt: false, meta: false });
    expect(e.at?.block).toBeGreaterThanOrEqual(0);
    expect(n.defaultPrevented).toBe(true);
    // After the call the payload keeps no DOM event.
    expect(e.native).not.toBe(n);
    // Where no track is (the data past the last row, the navigator): asked, not prevented by this subscriber.
    const empty = menu(find(t, p.x, HitRegion.data, false));
    const nav = menu(find(t, p.x, HitRegion.navigator, false));
    expect(asked.slice(1).map((a) => [a.part.kind, a.part.track])).toEqual([
      ['data', null],
      ['navigator', null],
    ]);
    expect([empty.defaultPrevented, nav.defaultPrevented]).toEqual([false, false]);
    t.ctx.dispose();
  });

  it("on a track's header: the same track", async () => {
    const { t, asked, menu } = await opened();
    const p = onTrack(t);
    menu(p);
    menu({ x: 20, y: p.y });
    expect(asked.map((a) => a.part.kind)).toEqual(['data', 'header']);
    expect(asked[1].part.track).toBe(asked[0].part.track);
    t.ctx.dispose();
  });

  it('with no subscriber nothing is built and nothing prevented: the browser shows its menu', async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 300] } });
    expect(t.ctx.interaction.observed('menu')).toBe(false);
    const m = Inputs.pointer();
    const n = native();
    Object.assign(m, { section: 0, ...onTrack(t), button: 2, pointerId: 1, native: n });
    t.ctx.input.menu(m);
    expect(n.defaultPrevented).toBe(false);
    t.ctx.dispose();
  });
});

describe('the track menu from the keyboard', () => {
  async function focused(prevent = true) {
    const o = await opened(() => prevent);
    o.t.focus(S);
    o.t.key(S, 'ArrowDown');
    return o;
  }

  it("Shift+F10 and ContextMenu ask at the cursor's cell, origin keyboard; handled where a subscriber prevented it", async () => {
    const { t, asked, base } = await focused();
    const c = t.ctx.ui.keyCursor.value;
    expect(c?.track).toBeTruthy();
    expect(t.key(S, 'F10', { shift: true })).toBe(true);
    expect(t.key(S, 'ContextMenu')).toBe(true);
    expect(base).toEqual([true, true]);
    expect(asked.map((e) => [e.origin, e.owner, e.part.section, e.part.track])).toEqual([
      ['keyboard', 'keyboard', S, c?.track],
      ['keyboard', 'keyboard', S, c?.track],
    ]);
    // At the cursor's cell: a point on the surface, in its bounds.
    const { x, y } = asked[0].place;
    expect([x >= 0 && x <= 1000, y >= 0 && y <= 300]).toEqual([true, true]);
    expect(t.ctx.ui.keyCursor.value).toEqual(c);
    t.ctx.dispose();
  });

  it("not prevented, the key is not handled (the browser's menu shows); F10 alone is the page's", async () => {
    const { t, asked } = await focused(false);
    expect(t.key(S, 'ContextMenu')).toBe(false);
    expect(asked).toHaveLength(1);
    expect(t.key(S, 'F10')).toBe(false);
    expect(asked).toHaveLength(1);
    t.ctx.dispose();
  });

  it("the browser's request after a key's is the same ask: not asked again, its default as the key's was", async () => {
    for (const prevent of [true, false]) {
      const { t, asked, menu } = await focused(prevent);
      t.key(S, 'ContextMenu');
      const n = menu({ x: 0, y: 0 }, { keyboard: true });
      expect(asked).toHaveLength(1);
      expect(n.defaultPrevented).toBe(prevent);
      // Only once: the next is a new ask.
      menu(onTrack(t));
      expect(asked).toHaveLength(2);
      t.ctx.dispose();
    }
  });

  it('a later key, a blur or a leave ends it too: a browser that sent no request swallows no later one', async () => {
    for (const end of ['key', 'blur', 'leave'] as const) {
      const { t, asked, menu } = await focused();
      // Prevented: this browser fires no context-menu request after it.
      expect(t.key(S, 'ContextMenu')).toBe(true);
      if (end === 'key') t.key(S, 'ArrowRight');
      else if (end === 'blur') t.focus(S, false);
      else t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
      if (end === 'blur') t.focus(S);
      const n = menu({ x: 0, y: 0 }, { keyboard: true });
      expect([end, asked.length, asked[1]?.origin]).toEqual([end, 2, 'keyboard']);
      expect(n.defaultPrevented).toBe(true);
      t.ctx.dispose();
    }
  });

  it("a press between ends it: the browser's next request is asked", async () => {
    const { t, asked, menu } = await focused();
    t.key(S, 'ContextMenu');
    const p = onTrack(t);
    t.pointer('down', S, p.x, p.y, { button: 2 });
    menu(p);
    t.pointer('up', S, p.x, p.y, { button: 2 });
    expect(asked.map((e) => e.origin)).toEqual(['keyboard', 'pointer']);
    t.ctx.dispose();
  });

  it("a keyboard request with no key before it (assistive technology) is asked at the cursor's cell", async () => {
    const { t, asked, menu } = await focused();
    const n = menu({ x: 0, y: 0 }, { keyboard: true });
    expect(asked).toHaveLength(1);
    expect(asked[0]).toMatchObject({ origin: 'keyboard', owner: 'keyboard' });
    expect(asked[0].part.track).toBe(t.ctx.ui.keyCursor.value?.track);
    expect(n.defaultPrevented).toBe(true);
    t.ctx.dispose();
  });
});

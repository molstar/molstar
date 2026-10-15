/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The public input streams (seqspec:interaction#streams): `ctx.interaction`'s payloads per part (refs, never indices),
// scoped by `on`, fed after the behaviors (the hover already set, the owner as claimed), observe-only, in Mol*'s button
// and modifier terms, the DOM event only during the call, and nothing built while a stream has no subscriber.

import { describe, expect, it, vi } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import type { Interaction } from '../src/interaction.js';
import {
  Hit,
  HitRegion,
  type SeqClickEvent,
  type SeqDragEvent,
  type SeqInputEvent,
  type SeqKeyEvent,
  type SeqNativeEvent,
  type SeqWheelEvent,
} from '../src/index.js';
import { type Driven, driven } from './input-fixtures.js';

const S = 'family.main';
const open = () => driven(case3(), { size: { [S]: [1000, 300] } });
const out = (t: Driven) => (t.ctx as unknown as { out: Interaction }).out;
/** The middle of the first scrolled track's data, in surface px. */
const onTrack = (t: Driven) => {
  const l = t.layout(S);
  const e = l.pinnedCount;
  return { x: l.dataLeft + 40, y: l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 };
};
const regionAt = (t: Driven, x: number, region: number) => {
  const h = Hit.create();
  for (let y = 1; y < t.layout(S).height; y++) {
    t.ctx.layout.hitTest(S, x, y, h);
    if (h.region === region) return { x, y };
  }
  throw new Error(`no point in region ${region}`);
};

describe('ctx.interaction', () => {
  it("a hover's payload per part: the region, the section's and track's refs, what is under the point", async () => {
    const t = await open();
    const seen: SeqInputEvent[] = [];
    t.ctx.interaction.hover.subscribe((e) => seen.push(e));
    const p = onTrack(t);
    t.pointer('move', S, p.x, p.y);
    t.pointer('move', S, 20, p.y);
    const nav = regionAt(t, p.x, HitRegion.navigator);
    t.pointer('move', S, nav.x, nav.y);
    expect(seen.map((e) => e.part.kind)).toEqual(['data', 'header', 'navigator']);
    const [data, header] = seen;
    expect(data).toMatchObject({ stream: 'hover', origin: 'pointer', place: { x: p.x, y: p.y } });
    expect(data.part.section).toBe(S);
    expect(data.part.track).toBeTruthy();
    expect(header.part.track).toBe(data.part.track);
    expect(data.at?.layer).toBeTruthy();
    expect(data.at?.block).toBeGreaterThanOrEqual(0);
    expect(data.at?.index).toBeGreaterThanOrEqual(0);
    expect(seen[2].part.track).toBeNull();
    // A payload is a consumer's to keep: each a fresh object.
    expect(new Set(seen).size).toBe(3);
    t.ctx.dispose();
  });

  it('`on` keeps a part, a section or a track; unscoped, the stream itself', async () => {
    const t = await open();
    const header: string[] = [];
    const elsewhere: string[] = [];
    const track: string[] = [];
    const p = onTrack(t);
    t.ctx.interaction.on('hover', { section: S, part: 'header' }).subscribe((e) => header.push(e.part.kind));
    t.ctx.interaction.on('hover', { section: 'nope' }).subscribe((e) => elsewhere.push(e.part.kind));
    t.pointer('move', S, p.x, p.y);
    const ref = t.ctx.pointer.current.inside ? (t.ctx.hover.live.current.track as number) : -1;
    expect(ref).toBeGreaterThanOrEqual(0);
    const name = t.ctx.graph.assembled.value?.set.strings.get(
      t.ctx.graph.assembled.value.view?.resolved.tracks.ref[ref] as number,
    ) as string;
    t.ctx.interaction.on('hover', { track: name, part: ['data', 'header'] }).subscribe((e) => track.push(e.part.kind));
    t.pointer('move', S, 20, p.y);
    t.pointer('move', S, p.x, p.y);
    t.pointer('move', S, p.x, regionAt(t, p.x, HitRegion.navigator).y);
    expect(header).toEqual(['header']);
    expect(elsewhere).toEqual([]);
    expect(track).toEqual(['header', 'data']);
    t.ctx.dispose();
  });

  it('after the behaviors: the hover already set, the owner as claimed; observe-only (no claim on a payload)', async () => {
    const t = await open();
    const hovered: string[] = [];
    const presses: SeqInputEvent[] = [];
    const drags: SeqDragEvent[] = [];
    const claimed: string[] = [];
    t.ctx.interaction.hover.subscribe(() => hovered.push(t.ctx.hover.live.current.owner));
    t.ctx.interaction.press.subscribe((e) => presses.push(e));
    t.ctx.interaction.drag.subscribe((e) => drags.push(e));
    t.ctx.input.events.pointer.subscribe((e) => e.input.type === 0 && claimed.push(e.owner), { priority: -100 });
    const nav = regionAt(t, 500, HitRegion.navigator);
    t.drag(S, nav.x, nav.y, nav.x + 40, nav.y);
    expect(hovered[0]).toBe('pointer');
    expect(presses.map((e) => e.owner)).toEqual(claimed);
    expect(claimed[0]).not.toBe('');
    expect(drags.map((e) => e.phase)).toEqual(['start', 'move', 'move', 'move', 'end']);
    expect(drags.every((e) => e.owner === claimed[0] && e.start.kind === 'navigator')).toBe(true);
    expect(drags[4].totalX).toBe(40);
    expect('claim' in presses[0]).toBe(false);
    t.ctx.dispose();
  });

  it("Mol*'s buttons and modifiers; a click's count; the DOM event only during the call", async () => {
    const t = await open();
    const clicks: SeqClickEvent[] = [];
    const presses: SeqInputEvent[] = [];
    const during: (SeqNativeEvent | null)[] = [];
    t.ctx.interaction.click.subscribe((e) => {
      clicks.push(e);
      during.push(e.native);
    });
    t.ctx.interaction.press.subscribe((e) => presses.push(e));
    const p = onTrack(t);
    t.click(S, p.x, p.y, { shift: true });
    t.click(S, p.x, p.y, { shift: true });
    t.pointer('down', S, p.x, p.y, { button: 1 });
    t.pointer('up', S, p.x, p.y, { button: 1 });
    t.pointer('down', S, p.x, p.y, { button: 2, ctrl: true, alt: true });
    t.pointer('up', S, p.x, p.y, { button: 2 });
    expect(presses.map((e) => [e.button, e.buttons])).toEqual([
      [1, 1],
      [1, 1],
      [4, 2],
      [2, 4],
    ]);
    expect(presses[3].modifiers).toEqual({ shift: false, control: true, alt: true, meta: false });
    expect(clicks[0].modifiers).toEqual({ shift: true, control: false, alt: false, meta: false });
    expect(clicks.slice(0, 2).map((e) => e.count)).toEqual([1, 2]);
    // Headless, a stand-in during the call; after it, another (a kept payload holds no DOM event).
    expect(during[0]).not.toBe(clicks[0].native);
    t.ctx.dispose();
  });

  it("a record's DOM event is the payload's during the call and let go after it", async () => {
    const t = await open();
    let seen: SeqNativeEvent | null = null;
    let kept: SeqInputEvent | null = null;
    t.ctx.interaction.hover.subscribe((e) => {
      seen = e.native;
      kept = e;
    });
    const n = { type: 'pointermove', defaultPrevented: false, preventDefault() {}, timeStamp: 1234.5 };
    const p = onTrack(t);
    const m = { type: 1, section: 0, x: p.x, y: p.y, buttons: 0, button: -1, modifiers: 0, time: 0.5 };
    t.ctx.input.pointer({ ...m, pointerId: 0, overX: 0, overY: 0, native: n });
    expect(seen).toBe(n);
    expect((kept as SeqInputEvent | null)?.native).not.toBe(n);
    expect((kept as SeqInputEvent | null)?.time).toBe(1234.5);
    t.ctx.dispose();
  });

  it('wheel, key and leave payloads', async () => {
    const t = await open();
    const wheels: SeqWheelEvent[] = [];
    const keys: SeqKeyEvent[] = [];
    const leaves: SeqInputEvent[] = [];
    t.ctx.interaction.wheel.subscribe((e) => wheels.push(e));
    t.ctx.interaction.key.subscribe((e) => keys.push(e));
    t.ctx.interaction.leave.subscribe((e) => leaves.push(e));
    const p = onTrack(t);
    t.wheel(S, p.x, p.y, 0, 30, { ctrl: true });
    t.focus(S);
    t.key(S, 'ArrowRight');
    t.key(S, 'q');
    t.ctx.input.leave(0);
    expect(wheels).toHaveLength(1);
    expect(wheels[0]).toMatchObject({ deltaX: 0, deltaY: 30, handled: true, part: { kind: 'data', section: S } });
    expect(wheels[0].modifiers.control).toBe(true);
    expect(keys.map((k) => [k.key, k.handled, k.origin, k.part.section])).toEqual([
      ['ArrowRight', true, 'keyboard', S],
      ['q', false, 'keyboard', S],
    ]);
    expect(keys[0].at).toBeNull();
    expect(leaves.map((e) => [e.stream, e.part.section, e.part.kind])).toEqual([['leave', S, 'none']]);
    t.ctx.dispose();
  });

  it("activate: a track header's double click (the `headers` behavior) and the keyboard's Enter", async () => {
    const t = await open();
    const seen: SeqInputEvent[] = [];
    t.ctx.interaction.activate.subscribe((e) => seen.push(e));
    const p = onTrack(t);
    t.click(S, 20, p.y);
    t.click(S, 20, p.y);
    expect(seen).toHaveLength(1);
    expect(seen[0]).toMatchObject({ origin: 'pointer', owner: 'headers', part: { kind: 'header', section: S } });
    expect(seen[0].part.track).toBeTruthy();
    t.focus(S);
    t.key(S, 'ArrowDown');
    t.key(S, 'Home');
    t.key(S, 'Enter');
    expect(seen.map((e) => e.origin)).toEqual(['pointer', 'keyboard']);
    expect(seen[1].part.track).toBe(t.ctx.ui.keyCursor.value?.track);
    t.ctx.dispose();
  });

  it('nothing is built while a stream has no subscriber, and again once it unsubscribed', async () => {
    const t = await open();
    const o = out(t);
    const spies = (['pointer', 'click', 'drag', 'wheel', 'key', 'leave'] as const).map((m) => vi.spyOn(o, m));
    const p = onTrack(t);
    const run = () => {
      t.pointer('move', S, p.x, p.y);
      t.pointer('move', S, 20, p.y);
      t.click(S, p.x, p.y);
      t.drag(S, p.x, p.y, p.x + 30, p.y);
      t.wheel(S, p.x, p.y, 0, 10);
      t.key(S, 'ArrowRight');
      t.ctx.input.leave(0);
    };
    run();
    expect(spies.map((s) => s.mock.calls.length)).toEqual([0, 0, 0, 0, 0, 0]);
    const sub = t.ctx.interaction.on('hover', { part: 'header' }).subscribe(() => {});
    expect(t.ctx.interaction.observed('hover')).toBe(true);
    run();
    expect(spies[0].mock.calls.length).toBeGreaterThan(0);
    sub.unsubscribe();
    expect(t.ctx.interaction.observed('hover')).toBe(false);
    for (const s of spies) s.mockClear();
    run();
    expect(spies.map((s) => s.mock.calls.length)).toEqual([0, 0, 0, 0, 0, 0]);
    t.ctx.dispose();
  });

  it("a UI's chrome reports through the same streams, origin chrome, only while observed", async () => {
    const t = await open();
    const seen: SeqClickEvent[] = [];
    const n = { type: 'click', defaultPrevented: false, preventDefault() {} };
    const part = { kind: 'legend' as const, section: null, track: null, group: null, key: 'colors:x' };
    const place = { x: 3, y: 4, clientX: 103, clientY: 54 };
    t.ctx.interaction.chrome('click', part, n, place, { count: 2 });
    t.ctx.interaction.click.subscribe((e) => seen.push(e));
    t.ctx.interaction.chrome('click', part, n, place, { count: 2, modifiers: 1 });
    expect(seen).toHaveLength(1);
    expect(seen[0]).toMatchObject({ origin: 'chrome', at: null, count: 2, place, part });
    expect(seen[0].modifiers.shift).toBe(true);
    t.ctx.dispose();
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the viewer publishes for a host's wiring (_spec/README.md "HoverModel", A Location; "AddressModel", Actions;
// "ViewModel", The last jump; "Composition", Retire): a hover of a Location with no hit, drawn projected with no card
// and owned (EC-141); a reveal by Location; the last jump's kind, Location and duration, a step back keeping the
// Location of the jump that left the view there; the retire's signal.

import { EmptyLocation, locationsEqual } from '@molstar/seqx-core/query';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';
import { EntStart, onChains, onUp, p53Structure, Up, UpSection } from './structure-fixture.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make() {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewer = v;
  await v.load(p53Structure());
  v.resize(800);
  await v.settle();
  frames.flush();
  return v;
}

describe('what the viewer publishes for the wiring', () => {
  it('setLocation: a hover with no hit, projected onto P04637, no card; owned (EC-141)', async () => {
    const v = await make();
    const link = {};
    const loc = onChains(10, 11, ['A']);
    const version = v.hover.version.value;
    expect(v.hover.setLocation(loc, link)).toBe(true);
    expect(v.hover.setLocation(loc, link)).toBe(false);
    expect(v.hover.version.value).toBe(version + 1);
    expect(v.hover.location()).toBe(loc);
    const up = v.hover.projectedAll()?.parts.find((p) => p.system === Up);
    expect(up?.ranges).toEqual(Int32Array.of(EntStart + 10, EntStart + 11));
    // no card, no words, no address: nothing was hit
    const card = new Describer({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe(v.hover.current);
    expect(card).toBeUndefined();
    expect(v.addresses.hoverAddressNow()).toBeUndefined();
    // the pointer takes it; the link's late clear does nothing
    const pointer = {};
    v.hover.set({ section: UpSection, system: Up, position: 5 }, pointer);
    expect(v.hover.clear(link)).toBe(false);
    expect(v.hover.owner).toBe(pointer);
    // an empty Location clears the owner's own; an undefined owner is refused
    v.hover.setLocation(loc, link);
    expect(v.hover.setLocation(loc, undefined)).toBe(false);
    expect(v.hover.setLocation(EmptyLocation, link)).toBe(true);
    expect(v.hover.active).toBe(false);
  });

  it('revealLocation: shown, nothing moves; else a pan, the zoom kept', async () => {
    const v = await make();
    v.addresses.focusLocation(onUp(100, 140));
    const before = v.view.get(UpSection)!;
    const jumps = v.view.jumps.value;
    expect(v.addresses.revealLocation(onUp(110, 120))).toBe('moved');
    expect(v.view.jumps.value).toBe(jumps);
    expect(v.addresses.revealLocation(onUp(200, 210))).toBe('moved');
    const after = v.view.get(UpSection)!;
    expect(after.end - after.start).toBeCloseTo(before.end - before.start);
    expect(after.start).toBeLessThanOrEqual(200);
    expect(after.end).toBeGreaterThanOrEqual(210);
    expect(v.view.lastJump).toMatchObject({ n: jumps + 1, kind: 'reveal' });
  });

  it('lastJump: kind, Location and duration, told before jumps; back names the earlier jump’s Location', async () => {
    const v = await make();
    const seen: unknown[] = [];
    const off = v.view.jumps.listen((n) => seen.push(v.view.lastJump?.n === n));
    const a = onUp(100, 140);
    v.addresses.focusLocation(a);
    expect(v.view.lastJump).toEqual({ n: v.view.jumps.value, kind: 'focus', location: a, durationMs: 240 });
    v.addresses.focusLocation(onUp(200, 220));
    v.view.back();
    expect(v.view.lastJump?.kind).toBe('back');
    expect(v.view.lastJump?.location).toBe(a);
    v.view.forward();
    expect(locationsEqual(v.view.lastJump!.location!, onUp(200, 220))).toBe(true);
    // a fit names none
    v.runCommand('view.fit', { section: UpSection });
    expect(v.view.lastJump).toMatchObject({ kind: 'fit' });
    expect(v.view.lastJump?.location).toBeUndefined();
    // a hand's move after a jump: the step pushed then names none
    v.addresses.focusLocation(a);
    v.view.pan(UpSection, 3);
    v.addresses.focusLocation(onUp(300, 320));
    v.view.back();
    expect(v.view.lastJump?.kind).toBe('back');
    expect(v.view.lastJump?.location).toBeUndefined();
    off();
    expect(seen.every((s) => s === true)).toBe(true);
  });

  it('lastJump names a step zoom `zoom` and a navigator’s move `move`, with no Location', async () => {
    const v = await make();
    v.addresses.focusLocation(onUp(100, 140));
    v.runCommand('view.zoom-out', { section: UpSection });
    expect(v.view.lastJump).toMatchObject({ kind: 'zoom', location: undefined });
    v.behavior<{ center(section: string, at: number): boolean }>('navigator')?.center(UpSection, 300);
    expect(v.view.lastJump).toMatchObject({ kind: 'move', location: undefined });
  });

  it('retiredSignal: aborted once at retire, the viewer already marked retired', async () => {
    const v = await make();
    const seen: boolean[] = [];
    expect(v.retiredSignal.aborted).toBe(false);
    v.retiredSignal.addEventListener('abort', () => seen.push(v.retired));
    v.retire();
    v.retire();
    expect(v.retiredSignal.aborted).toBe(true);
    expect(seen).toEqual([true]);
  });
});

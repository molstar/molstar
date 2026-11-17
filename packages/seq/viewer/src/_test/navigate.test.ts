/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '@molstar/ui/dom';
import { describe, expect, it, vi } from 'vitest';
import { navigateBehavior } from '../behaviors/navigate.js';
import { bindingTable, GenomeBrowserBindings } from '../bindings.js';
import { changeKindOf } from '../commands.js';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import { GeneSection, ProteinSection, proteinOnGene } from './fixture.js';

/** Frames through a patched `DomHost`: run by the test, never by a clock. */
function frames() {
  const queue: ((t: number) => void)[] = [];
  vi.spyOn(DomHost, 'requestAnimationFrame').mockImplementation((fn) => void queue.push(fn));
  return () => {
    for (const fn of queue.splice(0)) fn(0);
  };
}

function setup(config = viewerConfig(), bindings = bindingTable()) {
  const spec = new SpecModel(proteinOnGene({ geneLength: 24 }));
  const view = new ViewModel(spec, config);
  view.setWidth(ProteinSection, 40);
  view.setWidth(GeneSection, 40);
  const nav = navigateBehavior(view, bindings, config);
  const range = (s: string) => {
    const v = view.get(s)!;
    return [v.start, v.end];
  };
  return { view, nav, range };
}

describe('the navigate behavior', () => {
  it('wheel deltas within a frame are applied once, at the frame', () => {
    const frame = frames();
    const { view, nav, range } = setup();
    view.setRange(ProteinSection, 1, 3);
    // 2 positions over 40 px: 20 px a position
    nav.wheel(ProteinSection, 0, 10, 0, {}, 'cell');
    nav.wheel(ProteinSection, 0, 10, 0, {}, 'cell');
    expect(range(ProteinSection)).toEqual([1, 3]);
    frame();
    expect(range(ProteinSection)).toEqual([2, 4]);
  });

  it('Cmd/Ctrl+wheel and pinch zoom about the pointer; without coalescing a wheel applies at once', () => {
    const { view, nav, range } = setup(viewerConfig({ navigate: { coalesceWheel: false, wheelZoomStep: 2 } }));
    view.setRange(ProteinSection, 0, 4);
    // the pointer at 10 px is offset 1: it stays under the pointer
    nav.wheel(ProteinSection, 0, -100, 10, { mod: true }, 'cell');
    expect(range(ProteinSection)).toEqual([0.5, 2.5]);
    expect(nav.pinch(ProteinSection, 0.5, 20, {}, 'cell')).toBe(true);
    expect(range(ProteinSection)).toEqual([0, 4]);
  });

  it('a drag pans where it is bound to: on the axis by default, on data in the genome browser preset', () => {
    const { view, nav, range } = setup();
    view.setRange(ProteinSection, 1, 3);
    expect(nav.drag(ProteinSection, 20, {}, 'cell')).toBe(false);
    expect(nav.drag(ProteinSection, 20, {}, 'axis')).toBe(true);
    expect(range(ProteinSection)).toEqual([0, 2]);
    // a follower does not lead (EC-114)
    expect(nav.drag(GeneSection, 20, {}, 'axis')).toBe(false);

    const gb = setup(viewerConfig(), bindingTable(GenomeBrowserBindings));
    gb.view.setRange(ProteinSection, 1, 3);
    expect(gb.nav.drag(ProteinSection, -20, {}, 'cell')).toBe(true);
    expect(gb.range(ProteinSection)).toEqual([2, 4]);
  });

  it('the view’s moves are commands that navigate, with keys for the keymap', async () => {
    const { view, nav, range } = setup();
    const signal = new AbortController().signal;
    const byId = new Map(nav.commands.map((c) => [c.id, c]));
    for (const c of nav.commands) expect(changeKindOf(c.id)).toBe('navigate');
    view.setRange(ProteinSection, 1, 3);
    await byId.get('view.pan-right')!.run({}, { section: ProteinSection }, signal);
    expect(range(ProteinSection)[0]).toBeCloseTo(1.2, 10);
    expect(range(ProteinSection)[1]).toBeCloseTo(3.2, 10);
    await byId.get('view.fit')!.run({}, { section: ProteinSection }, signal);
    expect(range(ProteinSection)).toEqual([0, 4]);
    expect(() => byId.get('view.zoom-in')!.run({}, { section: GeneSection }, signal)).toThrow('cannot move');
    expect(() => byId.get('view.zoom-in')!.run({}, { section: undefined }, signal)).toThrow('No section');
  });
});

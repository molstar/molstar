// @vitest-environment jsdom

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { longPress, PressMs, PressSlop } from '../canvas-host.js';

// Touch (accessibility.md): a long press on the canvas highlights what is under the finger, as hover does for a mouse
// (the hover label shows it). A quick tap, a drag, another pointer kind and a lifted finger do not.

function setup() {
  const host = document.createElement('div');
  const context = { canvas: document.createElement('canvas') };
  host.append(context.canvas);
  document.body.append(host);
  const highlightOnly = vi.fn();
  // read asynchronously on every backend
  const identify = vi.fn(() => ({ wait: () => Promise.resolve({ id: 7 }) }));
  const loci = { kind: 'element-loci' };
  const plugin = {
    canvas3dContext: context,
    canvas3d: { asyncIdentify: identify, getLoci: vi.fn(() => ({ loci })) },
    managers: { interactivity: { lociHighlights: { highlightOnly } } },
  } as unknown as PluginContext;
  longPress(plugin);
  const pointer = (type: string, init: { x?: number; y?: number; pointerType?: string } = {}) => {
    const e = new Event(type, { bubbles: true }) as Event & Record<string, unknown>;
    e.pointerType = init.pointerType ?? 'touch';
    e.clientX = init.x ?? 10;
    e.clientY = init.y ?? 20;
    context.canvas.dispatchEvent(e);
  };
  // a backend fallback replaces the canvas element in place
  const replaceCanvas = () => {
    const next = document.createElement('canvas');
    context.canvas.replaceWith(next);
    context.canvas = next;
  };
  return { pointer, replaceCanvas, highlightOnly, identify, loci };
}

describe('long press on the canvas', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => {
    vi.useRealTimers();
    document.body.innerHTML = '';
  });

  it('highlights what is under a finger held still', async () => {
    const t = setup();
    t.pointer('pointerdown');
    await vi.advanceTimersByTimeAsync(PressMs - 1);
    expect(t.highlightOnly).not.toHaveBeenCalled();
    await vi.advanceTimersByTimeAsync(1);
    expect(t.identify).toHaveBeenCalled();
    expect(t.highlightOnly).toHaveBeenCalledWith({ loci: t.loci });
  });

  it('does not for a tap, a lifted finger, a drag past the slop, or a mouse', async () => {
    const t = setup();
    t.pointer('pointerdown');
    t.pointer('pointerup');
    await vi.advanceTimersByTimeAsync(PressMs * 2);
    t.pointer('pointerdown');
    t.pointer('pointermove', { x: 10 + PressSlop + 1 });
    await vi.advanceTimersByTimeAsync(PressMs * 2);
    t.pointer('pointerdown', { pointerType: 'mouse' });
    await vi.advanceTimersByTimeAsync(PressMs * 2);
    expect(t.highlightOnly).not.toHaveBeenCalled();
  });

  it('allows a small tremor of the finger', async () => {
    const t = setup();
    t.pointer('pointerdown');
    t.pointer('pointermove', { x: 10 + PressSlop - 1 });
    await vi.advanceTimersByTimeAsync(PressMs);
    expect(t.highlightOnly).toHaveBeenCalledTimes(1);
  });

  it('follows the canvas a backend fallback puts in its place', async () => {
    const t = setup();
    t.replaceCanvas();
    t.pointer('pointerdown');
    await vi.advanceTimersByTimeAsync(PressMs);
    expect(t.highlightOnly).toHaveBeenCalledWith({ loci: t.loci });
  });
});

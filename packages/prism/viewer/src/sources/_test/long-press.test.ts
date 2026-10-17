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
  const canvas = document.createElement('canvas');
  document.body.append(canvas);
  const highlightOnly = vi.fn();
  const identify = vi.fn(() => ({ id: 7 }));
  const loci = { kind: 'element-loci' };
  const plugin = {
    canvas3dContext: { canvas },
    canvas3d: { identify, getLoci: vi.fn(() => ({ loci })) },
    managers: { interactivity: { lociHighlights: { highlightOnly } } },
  } as unknown as PluginContext;
  longPress(plugin);
  const pointer = (type: string, init: { x?: number; y?: number; pointerType?: string } = {}) => {
    const e = new Event(type, { bubbles: true }) as Event & Record<string, unknown>;
    e.pointerType = init.pointerType ?? 'touch';
    e.clientX = init.x ?? 10;
    e.clientY = init.y ?? 20;
    canvas.dispatchEvent(e);
  };
  return { pointer, highlightOnly, identify, loci };
}

describe('long press on the canvas', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => {
    vi.useRealTimers();
    document.body.innerHTML = '';
  });

  it('highlights what is under a finger held still', () => {
    const t = setup();
    t.pointer('pointerdown');
    vi.advanceTimersByTime(PressMs - 1);
    expect(t.highlightOnly).not.toHaveBeenCalled();
    vi.advanceTimersByTime(1);
    expect(t.identify).toHaveBeenCalled();
    expect(t.highlightOnly).toHaveBeenCalledWith({ loci: t.loci });
  });

  it('does not for a tap, a lifted finger, a drag past the slop, or a mouse', () => {
    const t = setup();
    t.pointer('pointerdown');
    t.pointer('pointerup');
    vi.advanceTimersByTime(PressMs * 2);
    t.pointer('pointerdown');
    t.pointer('pointermove', { x: 10 + PressSlop + 1 });
    vi.advanceTimersByTime(PressMs * 2);
    t.pointer('pointerdown', { pointerType: 'mouse' });
    vi.advanceTimersByTime(PressMs * 2);
    expect(t.highlightOnly).not.toHaveBeenCalled();
  });

  it('allows a small tremor of the finger', () => {
    const t = setup();
    t.pointer('pointerdown');
    t.pointer('pointermove', { x: 10 + PressSlop - 1 });
    vi.advanceTimersByTime(PressMs);
    expect(t.highlightOnly).toHaveBeenCalledTimes(1);
  });
});

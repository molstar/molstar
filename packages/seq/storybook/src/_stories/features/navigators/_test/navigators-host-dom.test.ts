// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The host-drawn navigator's DOM (the HostDrawn story's `wireHostNavigator`) in jsdom on a fake clock: its elements
// painted from the model, and a mouse drag over its window, zoomed in first as a reader does, panning the view through
// the host's pointer events alone.

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { HostSection, hostDrawnSpec } from '../compose.js';
import { wireHostNavigator } from '../host-navigator-dom.js';

const raf = DomHost.requestAnimationFrame;
let viewer: SeqViewer | undefined;
let unwire: (() => void) | undefined;
let nav: HTMLElement | undefined;

beforeEach(() => {
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  unwire?.();
  viewer?.retire();
  nav?.remove();
  unwire = viewer = nav = undefined;
  DomHost.requestAnimationFrame = raf;
  vi.useRealTimers();
  vi.restoreAllMocks();
});

/** A pointer event as a mouse sends it (jsdom has no PointerEvent: a MouseEvent of the type, with its pointer id). */
function pointer(type: string, x: number, buttons = 1): MouseEvent {
  const e = new MouseEvent(type, { clientX: x, clientY: 18, button: 0, buttons, bubbles: true, cancelable: true });
  return Object.assign(e, { pointerId: 1, pointerType: 'mouse' });
}

const windowEl = () => nav!.querySelector<HTMLElement>('[data-window]')!;
const left = (el: HTMLElement) => Number.parseFloat(el.style.left);

describe('the host-drawn navigator in a DOM', () => {
  it('a mouse drag over its window pans the view, after a zoom', async () => {
    const v = (viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined }));
    const model = v.navigator({ id: 'host', of: HostSection, content: 'axis' });
    nav = document.createElement('div');
    document.body.append(nav);
    unwire = wireHostNavigator(nav, v, model, 36);
    v.resize(720);
    void v.load(hostDrawnSpec());
    await vi.advanceTimersByTimeAsync(2000);
    v.view.setRange(HostSection, 159, 207);
    v.draw();
    await vi.advanceTimersByTimeAsync(600);
    // the host's elements: six domains and the window, placed through the model
    expect(nav.querySelectorAll('[data-block]').length).toBe(6);
    const w = windowEl();
    expect(left(w)).toBeCloseTo((159 / 393) * 600, 0);
    const mid = left(w) + Number.parseFloat(w.style.width) / 2;
    // press, move, release over the window: the view moves with the hand
    nav.dispatchEvent(pointer('pointerdown', mid));
    nav.dispatchEvent(pointer('pointermove', mid + (30 / 393) * 600));
    nav.dispatchEvent(pointer('pointerup', mid + (30 / 393) * 600, 0));
    nav.dispatchEvent(pointer('click', mid + (30 / 393) * 600, 0));
    expect(v.view.get(HostSection)!.start).toBeCloseTo(189, 0);
    v.draw();
    await vi.advanceTimersByTimeAsync(600);
    expect(left(windowEl())).toBeCloseTo((189 / 393) * 600, 0);
  });
});

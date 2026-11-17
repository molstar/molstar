// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { viewerItem } from '../config.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';
import { proteinOnGene } from './fixture.js';

// The mounted canvas stack in a DOM (jsdom): real events through the container's listeners, and what remains of a
// retired viewer in the document, its listeners and its observers.

let frames: Frames;
let viewer: SeqViewer | undefined;
const listeners = { added: 0, removed: 0 };
const observers = { observed: 0, disconnected: 0 };

class FakeResizeObserver {
  observe() {
    observers.observed++;
  }
  unobserve() {}
  disconnect() {
    observers.disconnected++;
  }
}

beforeEach(() => {
  frames = installFrames();
  listeners.added = listeners.removed = 0;
  observers.observed = observers.disconnected = 0;
  const add = EventTarget.prototype.addEventListener;
  const remove = EventTarget.prototype.removeEventListener;
  vi.spyOn(EventTarget.prototype, 'addEventListener').mockImplementation(function (
    this: EventTarget,
    ...a: [any, any, any]
  ) {
    if (this instanceof HTMLElement) listeners.added++;
    return add.apply(this, a);
  });
  vi.spyOn(EventTarget.prototype, 'removeEventListener').mockImplementation(function (
    this: EventTarget,
    ...a: [any, any, any]
  ) {
    if (this instanceof HTMLElement) listeners.removed++;
    return remove.apply(this, a);
  });
  vi.stubGlobal('ResizeObserver', FakeResizeObserver);
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    const ctx = ((this as any).__ctx ??= new RecordingContext());
    return ctx;
  });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

async function mount() {
  const target = document.createElement('div');
  document.body.appendChild(target);
  viewer = createSeqViewer(
    { ...DefaultSeqViewerSpec(), config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] },
    target,
  );
  await viewer.load(proteinOnGene());
  await viewer.settle();
  frames.flush();
  return { target, viewer, container: target.firstElementChild as HTMLElement };
}

const at = (type: string, x: number, y: number, init: MouseEventInit = {}) =>
  new MouseEvent(type, { clientX: x, clientY: y, bubbles: true, cancelable: true, ...init });

describe('the mounted canvas stack', () => {
  it('makes a focusable container with the base canvas and the overlay above it, sized to the layout', async () => {
    const { container, viewer } = await mount();
    const [base, overlay] = Array.from(container.querySelectorAll('canvas'));
    expect(container.tabIndex).toBe(0);
    expect(overlay.style.pointerEvents).toBe('none');
    expect(base.width).toBe(600 * viewer.dpr);
    expect(base.style.height).toBe(`${viewer.height}px`);
    expect((base as any).__ctx.calls('fillRect').length).toBeGreaterThan(3);
  });

  it('drives the viewer with real events: a Ctrl+wheel zooms, a pointer move hovers, a drag pans', async () => {
    const { container, viewer } = await mount();
    viewer.input.mac = false;
    const wheel = new WheelEvent('wheel', {
      clientX: 300,
      clientY: 18,
      deltaY: -100,
      deltaMode: 1,
      ctrlKey: true,
      bubbles: true,
      cancelable: true,
    });
    container.dispatchEvent(wheel);
    // a wheel a binding took is not the page's to scroll
    expect(wheel.defaultPrevented).toBe(true);
    frames.flush();
    const v = viewer.view.get('/views/main/protein')!;
    expect(v.end - v.start).toBeLessThan(4);
    container.dispatchEvent(at('pointermove', 300, 18));
    expect(viewer.hover.active).toBe(true);
    container.dispatchEvent(at('pointerleave', 300, 18));
    expect(viewer.hover.active).toBe(false);
    // a drag from the label column
    const start = viewer.view.get('/views/main/protein')!.start;
    container.dispatchEvent(at('pointerdown', 40, 18, { button: 0 }));
    container.dispatchEvent(at('pointermove', 100, 18));
    container.dispatchEvent(at('pointerup', 100, 18));
    frames.flush();
    expect(viewer.view.get('/views/main/protein')!.start).toBeLessThan(start);
  });

  it('keys run commands; Space held makes a drag pan on data', async () => {
    const { container, viewer } = await mount();
    viewer.input.mac = false;
    const key = new KeyboardEvent('keydown', { key: '=', bubbles: true, cancelable: true });
    container.dispatchEvent(key);
    expect(key.defaultPrevented).toBe(true);
    frames.flush();
    const v = viewer.view.get('/views/main/protein')!;
    expect(v.end - v.start).toBeLessThan(4);
    expect(viewer.input.down(300, 18)).toBe(false);
    container.dispatchEvent(new KeyboardEvent('keydown', { key: ' ', bubbles: true }));
    expect(viewer.input.down(300, 18)).toBe(true);
    viewer.input.up();
    container.dispatchEvent(new KeyboardEvent('keyup', { key: ' ', bubbles: true }));
    expect(viewer.input.down(300, 18)).toBe(false);
  });

  it('retire takes away every listener it added, the observer and the nodes, and no frame or timer is left', async () => {
    const { target, container, viewer } = await mount();
    container.dispatchEvent(at('pointermove', 300, 18));
    container.dispatchEvent(
      new WheelEvent('wheel', {
        clientX: 300,
        clientY: 18,
        deltaY: -100,
        deltaMode: 1,
        ctrlKey: true,
        cancelable: true,
      }),
    );
    expect(listeners.added).toBeGreaterThan(5);
    expect(observers.observed).toBe(1);
    viewer.retire();
    expect(listeners.removed).toBe(listeners.added);
    expect(observers.disconnected).toBe(1);
    expect(target.childElementCount).toBe(0);
    const drawn = viewer.counts.frames;
    frames.flush();
    expect(viewer.counts.frames).toBe(drawn);
    expect(frames.timers.live).toBe(0);
    // an event on the removed container does nothing
    const version = viewer.hover.version.value;
    container.dispatchEvent(at('pointermove', 100, 18));
    expect(viewer.hover.version.value).toBe(version);
  });

  it('a spec without a UI leaves the target empty', () => {
    const target = document.createElement('div');
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined }, target);
    expect(target.childElementCount).toBe(0);
    v.retire();
  });
});

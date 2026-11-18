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
/** The media queries asked for and their change listeners, as `matchMedia` made them. */
let queries: { media: string; listeners: Set<() => void> }[] = [];

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
  queries = [];
  // jsdom has none: defined on the document's window, which the mount reads
  const matchMedia = (media: string) => {
    const q = { media, listeners: new Set<() => void>() };
    queries.push(q);
    return {
      media,
      matches: true,
      addEventListener: (_: string, fn: () => void) => q.listeners.add(fn),
      removeEventListener: (_: string, fn: () => void) => q.listeners.delete(fn),
    } as unknown as MediaQueryList;
  };
  Object.defineProperty(document.defaultView!, 'matchMedia', { configurable: true, writable: true, value: matchMedia });
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
  Reflect.deleteProperty(document.defaultView!, 'matchMedia');
  Object.defineProperty(document.defaultView!, 'devicePixelRatio', { configurable: true, value: 1 });
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
  it('has the standalone viewer’s history when mounted from the default spec; none headless', async () => {
    const { viewer: v } = await mount();
    expect(v.spec.history).toBeDefined();
    const headless = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    expect(headless.spec.history).toBeUndefined();
    headless.retire();
  });

  it('makes a focusable container with the base canvas and the overlay above it, sized to the layout', async () => {
    const { container, viewer } = await mount();
    const [base, overlay, labels] = Array.from(container.querySelectorAll('canvas'));
    expect(container.tabIndex).toBe(0);
    expect(overlay.style.pointerEvents).toBe('none');
    expect(labels.style.pointerEvents).toBe('none');
    expect(labels.width).toBe(base.width);
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
    // on data a drag selects; with Space held it pans
    const mods = (space: boolean) => ({ shift: false, mod: false, alt: false, space });
    expect(viewer.bindings.resolve('drag', mods(false), 'cell')).toBe('drag-select');
    container.dispatchEvent(new KeyboardEvent('keydown', { key: ' ', bubbles: true }));
    const start = viewer.view.get('/views/main/protein')!.start;
    expect(viewer.input.down(300, 18)).toBe(true);
    viewer.input.move(340, 18);
    viewer.input.up();
    expect(viewer.view.get('/views/main/protein')!.start).not.toBe(start);
    // a Space that panned is no tap: it toggles nothing
    container.dispatchEvent(new KeyboardEvent('keyup', { key: ' ', bubbles: true }));
    expect(viewer.selection.value.value.shown.parts).toEqual([]);
  });

  it('a drag-select survives the platform’s real sequence: up, capture lost, then a click', async () => {
    const { container, viewer } = await mount();
    for (const [type, x] of [
      ['pointerdown', 300],
      ['gotpointercapture', 300],
      ['pointermove', 340],
      ['pointermove', 420],
      ['pointerup', 420],
      ['lostpointercapture', 420],
    ] as const)
      container.dispatchEvent(at(type, x, 18, { button: 0 }));
    const dragged = [...viewer.selection.value.value.shown.parts[0].ranges];
    expect(dragged[1] - dragged[0]).toBeGreaterThan(1);
    container.dispatchEvent(at('click', 420, 18, { button: 0, detail: 1 }));
    expect([...viewer.selection.value.value.shown.parts[0].ranges]).toEqual(dragged);
    // the next press without a move is a click again
    container.dispatchEvent(at('pointerdown', 300, 18, { button: 0 }));
    container.dispatchEvent(at('pointerup', 300, 18, { button: 0 }));
    container.dispatchEvent(at('lostpointercapture', 300, 18));
    container.dispatchEvent(at('click', 300, 18, { button: 0, detail: 1 }));
    expect(viewer.selection.value.value.shown.parts[0].ranges.length).toBe(2);
    expect(
      viewer.selection.value.value.shown.parts[0].ranges[1] - viewer.selection.value.value.shown.parts[0].ranges[0],
    ).toBe(1);
  });

  it('a drag ends when the pointer is cancelled or its capture is lost', async () => {
    const { container, viewer } = await mount();
    for (const end of ['pointercancel', 'lostpointercapture']) {
      container.dispatchEvent(at('pointerdown', 40, 18, { button: 0 }));
      expect(viewer.input.dragging).toBe(true);
      container.dispatchEvent(at(end, 40, 18));
      expect(viewer.input.dragging).toBe(false);
    }
  });

  it('Space does not scroll the page while it is the pan modifier and the canvas has focus', async () => {
    const { container, viewer } = await mount();
    expect(viewer.input.spaceIsModifier).toBe(true);
    const space = new KeyboardEvent('keydown', { key: ' ', bubbles: true, cancelable: true });
    container.dispatchEvent(space);
    expect(space.defaultPrevented).toBe(true);
    // a key no command has is the page's
    const other = new KeyboardEvent('keydown', { key: 'q', bubbles: true, cancelable: true });
    container.dispatchEvent(other);
    expect(other.defaultPrevented).toBe(false);
  });

  it('a change of the device pixel ratio asks for a frame and listens for the next ratio', async () => {
    const { viewer, container } = await mount();
    expect(queries.length).toBe(1);
    expect(queries[0].media).toBe('(resolution: 1dppx)');
    expect(queries[0].listeners.size).toBe(1);
    frames.flush();
    const drawn = viewer.counts.frames;
    // the ratio changes: the listener of the old query goes, a query for the new ratio is made, a frame follows
    vi.stubGlobal('devicePixelRatio', 2);
    Object.defineProperty(document.defaultView!, 'devicePixelRatio', { configurable: true, value: 2 });
    for (const fn of [...queries[0].listeners]) fn();
    expect(queries[0].listeners.size).toBe(0);
    expect(queries.length).toBe(2);
    expect(queries[1].media).toBe('(resolution: 2dppx)');
    expect(queries[1].listeners.size).toBe(1);
    frames.flush();
    expect(viewer.counts.frames).toBe(drawn + 1);
    expect(viewer.dpr).toBe(2);
    expect(container.querySelector('canvas')!.width).toBe(viewer.width * 2);
    viewer.retire();
    expect(queries[1].listeners.size).toBe(0);
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

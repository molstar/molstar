/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DomHost, pendingFrameCallbacks, render, ViewScope } from '../../../dom/index.js';
import { fakeContext, stubPath2D } from '../../logo/_test/fake-canvas.js';
import { Splash } from '../splash.js';

// The splash (stage §4): it assembles when it shows, draws only while something moves, none at rest, off screen or
// hidden; under reduced motion one still frame.

describe('Splash', () => {
  const host = DomHost.requestAnimationFrame;
  let queue: ((t: number) => void)[] = [];
  let now = 0;
  const frames = (ms: number) => {
    for (let t = 0; t < ms && queue.length; t += 16) {
      now += 16;
      const run = queue;
      queue = [];
      for (const f of run) f(now);
    }
  };
  let counts: ReturnType<typeof fakeContext>['counts'];
  let reduced = false;
  /** The observed element's visibility, as an IntersectionObserver reports it. */
  let intersect: ((on: boolean) => void) | undefined;

  beforeEach(() => {
    stubPath2D();
    reduced = false;
    vi.stubGlobal('matchMedia', () => ({ matches: reduced }));
    vi.stubGlobal(
      'IntersectionObserver',
      class {
        constructor(fn: IntersectionObserverCallback) {
          intersect = (on) => fn([{ isIntersecting: on } as IntersectionObserverEntry], this as never);
        }
        observe() {}
        disconnect() {
          intersect = undefined;
        }
      },
    );
    const fake = fakeContext();
    counts = fake.counts;
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockReturnValue(fake.ctx as never);
    DomHost.requestAnimationFrame = (f) => {
      queue.push(f);
      return 0;
    };
  });
  afterEach(() => {
    document.body.innerHTML = '';
    frames(160);
    DomHost.requestAnimationFrame = host;
    queue = [];
  });

  function mount(scope = new ViewScope()) {
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <Splash scale={2} />, container, scope);
    return { container, dispose, scope };
  }

  it('assembles when it shows, then rests with no frames; again when shown again', () => {
    const { container, dispose, scope } = mount();
    expect(container.querySelector('[role="img"]')?.getAttribute('aria-label')).toBe('Mol*');
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(0);
    expect(pendingFrameCallbacks()).toBeGreaterThan(0);
    // The splash's gathering: 2.4 s and its stagger.
    frames(4000);
    expect(pendingFrameCallbacks()).toBe(0);
    const drawn = counts.drawImage;
    frames(500);
    expect(counts.drawImage).toBe(drawn);
    scope.setActive(false);
    expect(pendingFrameCallbacks()).toBe(0);
    scope.setActive(true);
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(drawn);
    dispose();
  });

  it('off screen, asks for no frames; back on screen, goes on', () => {
    const { dispose } = mount();
    frames(48);
    intersect?.(false);
    expect(pendingFrameCallbacks()).toBe(0);
    const drawn = counts.drawImage;
    frames(500);
    expect(counts.drawImage).toBe(drawn);
    intersect?.(true);
    expect(pendingFrameCallbacks()).toBeGreaterThan(0);
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(drawn);
    dispose();
  });

  it('off screen before its first frame: it assembles when it comes on screen', () => {
    const { dispose } = mount();
    intersect?.(false);
    frames(500);
    expect(counts.drawImage).toBe(0);
    intersect?.(true);
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(0);
    dispose();
  });

  it('a press shatters it; it gathers again and rests', () => {
    const { container, dispose } = mount();
    frames(4000);
    const word = container.querySelector('[role="img"]')!;
    word.dispatchEvent(new PointerEvent('pointerdown', { clientX: 40, clientY: 20, cancelable: true }));
    expect(pendingFrameCallbacks()).toBeGreaterThan(0);
    frames(3000);
    expect(pendingFrameCallbacks()).toBe(0);
    dispose();
  });

  it('under reduced motion: one still frame', () => {
    reduced = true;
    const { container, dispose } = mount();
    frames(16);
    expect(pendingFrameCallbacks()).toBe(0);
    const drawn = counts.drawImage + counts.stroke;
    expect(drawn).toBeGreaterThan(0);
    container.querySelector('[role="img"]')!.dispatchEvent(new PointerEvent('pointerdown', { clientX: 4 }));
    expect(pendingFrameCallbacks()).toBe(0);
    dispose();
  });
});

/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DomHost, pendingFrameCallbacks, render, ViewScope } from '../../../dom/index.js';
import { GlassWordmark } from '../glass-wordmark.js';
import { fakeContext, stubPath2D } from './fake-canvas.js';

// The heading's glass (glass-wordmark.tsx): it assembles each time it shows, asks for frames only while something moves,
// none while hidden or at rest.

describe('GlassWordmark', () => {
  const host = DomHost.requestAnimationFrame;
  let queue: ((t: number) => void)[] = [];
  let now = 0;
  /** Runs the page's frames for `ms` (60 a second), as the browser would. */
  const frames = (ms: number) => {
    for (let t = 0; t < ms && queue.length; t += 16) {
      now += 16;
      const run = queue;
      queue = [];
      for (const f of run) f(now);
    }
  };
  let counts: ReturnType<typeof fakeContext>['counts'];

  beforeEach(() => {
    stubPath2D();
    vi.stubGlobal('matchMedia', () => ({ matches: false }));
    const fake = fakeContext();
    counts = fake.counts;
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockReturnValue(fake.ctx as any);
    DomHost.requestAnimationFrame = (f) => void queue.push(f);
  });
  afterEach(() => {
    // The frames asked for run (not dropped: the runtime would wait for a tick that never comes), then the page's own.
    document.body.innerHTML = '';
    frames(160);
    DomHost.requestAnimationFrame = host;
    queue = [];
  });

  it('assembles when it shows, rests (no frames), none while hidden, assembles again when shown', () => {
    const scope = new ViewScope();
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <GlassWordmark />, container, scope);
    expect(container.querySelector('[role="img"]')?.getAttribute('aria-label')).toBe('Mol*');
    // Shown: its first frame scatters it, the next draw the pieces coming together.
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(0);
    expect(pendingFrameCallbacks()).toBeGreaterThan(0);
    // Assembled (0.9 s and its stagger): whole, at rest, nothing asked for.
    frames(2000);
    expect(pendingFrameCallbacks()).toBe(0);
    const drawn = counts.drawImage;
    frames(500);
    expect(counts.drawImage).toBe(drawn);
    // Hidden (its view folded): nothing; shown again: it assembles again.
    scope.setActive(false);
    expect(pendingFrameCallbacks()).toBe(0);
    scope.setActive(true);
    frames(48);
    expect(counts.drawImage).toBeGreaterThan(drawn);
    dispose();
  });

  it('a pointer resting over it asks for no frames once the pieces settle', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <GlassWordmark />, container, new ViewScope());
    frames(2000);
    const word = container.querySelector('[role="img"]')!;
    word.dispatchEvent(new PointerEvent('pointerenter'));
    word.dispatchEvent(new PointerEvent('pointermove', { clientX: 20, clientY: 10 }));
    expect(pendingFrameCallbacks()).toBeGreaterThan(0);
    // The pieces lean and settle (under a second); held still under the pointer, nothing asks for a frame.
    frames(2000);
    expect(pendingFrameCallbacks()).toBe(0);
    dispose();
  });

  it('a press on it starts no text selection (default prevented) and still reaches what listens around it', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <GlassWordmark />, container, new ViewScope());
    let heard = false;
    container.addEventListener('pointerdown', () => (heard = true));
    const press = new PointerEvent('pointerdown', { bubbles: true, cancelable: true, clientX: 20, clientY: 10 });
    container.querySelector('[role="img"]')!.dispatchEvent(press);
    expect(press.defaultPrevented).toBe(true);
    expect(heard).toBe(true);
    dispose();
  });
});

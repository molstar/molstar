/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { EmptyCanvas } from '../empty-canvas.js';

describe('EmptyCanvas', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  function mount(initial: boolean) {
    const empty = observableState(initial);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} />, container);
    const mark = () => (container.querySelector('[role=img]')?.parentElement ?? null) as HTMLElement | null;
    return { empty, mark };
  }

  it('shows the mark only after the canvas has been empty a while, then fades it out when something loads', () => {
    const { empty, mark } = mount(true);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(400);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(400);
    expect(mark()?.dataset.phase).toBe('in');

    empty.set(false);
    expect(mark()?.dataset.phase).toBe('out');
    vi.advanceTimersByTime(300);
    expect(mark()).toBeNull();
  });

  it('never shows when something loads before the delay (clearing a scene to load the next)', () => {
    const { empty, mark } = mount(true);
    vi.advanceTimersByTime(300);
    empty.set(false);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(1000);
    expect(mark()).toBeNull();
  });

  it('shows the splash under the mark with it, once the viewer is idle on the empty canvas', () => {
    const empty = observableState(true);
    const busy = observableState(false);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} busy={busy.asReadonly} />, container);
    const splash = () => container.querySelector('canvas');
    // Working (a load starting): the mark shows, the splash waits.
    busy.set(true);
    vi.advanceTimersByTime(800);
    expect(container.querySelector('[role=img]')).not.toBeNull();
    expect(splash()).toBeNull();
    vi.advanceTimersByTime(5000);
    expect(splash()).toBeNull();
    // Idle: with the mark.
    busy.set(false);
    vi.advanceTimersByTime(1);
    expect(splash()).not.toBeNull();
    // Something loaded: it goes with the mark.
    empty.set(false);
    expect(splash()).not.toBeNull();
    vi.advanceTimersByTime(300);
    expect(splash()).toBeNull();
  });

  it('shows no splash when the spec turns it off', () => {
    const empty = observableState(true);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} splash={false} />, container);
    vi.advanceTimersByTime(10_000);
    expect(container.querySelector('[role=img]')).not.toBeNull();
    expect(container.querySelector('canvas')).toBeNull();
  });
});

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
import { runInScope, Scope } from '../../../dom/scope.js';
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

  it('holds its wait while its scope is paused, and does not show a mark from hidden UI', () => {
    const empty = observableState(true);
    const scope = new Scope();
    const el = runInScope(scope, () => <EmptyCanvas empty={empty.asReadonly} />) as HTMLElement;
    const mark = () => el.querySelector('[role=img]');
    vi.advanceTimersByTime(300);
    scope.setActive(false);
    vi.advanceTimersByTime(5000);
    expect(mark()).toBeNull();
    scope.setActive(true);
    vi.advanceTimersByTime(650);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(100);
    expect(mark()).not.toBeNull();
    scope.dispose();
  });

  it('never shows when something loads before the delay (clearing a scene to load the next)', () => {
    const { empty, mark } = mount(true);
    vi.advanceTimersByTime(300);
    empty.set(false);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(1000);
    expect(mark()).toBeNull();
  });

  it('waits for the viewer to be idle: work restarts the wait, and once shown work alone does not hide it', () => {
    const empty = observableState(true);
    const busy = observableState(true);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} busy={busy.asReadonly} />, container);
    const mark = () => container.querySelector('[role=img]');
    const splash = () => container.querySelector('canvas');
    // Empty but working (a long load): nothing, however long.
    vi.advanceTimersByTime(5000);
    expect(mark()).toBeNull();
    // Idle a moment, then working again: the wait starts over.
    busy.set(false);
    vi.advanceTimersByTime(500);
    busy.set(true);
    vi.advanceTimersByTime(100);
    busy.set(false);
    vi.advanceTimersByTime(500);
    expect(mark()).toBeNull();
    vi.advanceTimersByTime(250);
    // Idle the whole delay: the mark, and the splash with it.
    expect(mark()).not.toBeNull();
    expect(splash()).not.toBeNull();
    // Working while shown: it stays.
    busy.set(true);
    vi.advanceTimersByTime(2000);
    expect(splash()).not.toBeNull();
    // Something loaded: both go.
    empty.set(false);
    expect(splash()).not.toBeNull();
    vi.advanceTimersByTime(300);
    expect(splash()).toBeNull();
    expect(mark()).toBeNull();
  });

  it('shows no splash when the spec turns it off', () => {
    const empty = observableState(true);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} splash={false} />, container);
    vi.advanceTimersByTime(10_000);
    expect(container.querySelector('[role=img]')).not.toBeNull();
    expect(container.querySelector('canvas')).toBeNull();
  });

  it('shows particles one time and waves another: a fresh draw each time it shows', () => {
    const random = vi.spyOn(Math, 'random');
    const empty = observableState(true);
    const container = document.createElement('div');
    render(() => <EmptyCanvas empty={empty.asReadonly} />, container);
    random.mockClear();
    vi.advanceTimersByTime(800);
    expect(container.querySelector('canvas')).not.toBeNull();
    const draws = random.mock.calls.length;
    expect(draws).toBeGreaterThan(0);
    // Gone and back: drawn again.
    empty.set(false);
    vi.advanceTimersByTime(300);
    empty.set(true);
    vi.advanceTimersByTime(800);
    expect(container.querySelector('canvas')).not.toBeNull();
    expect(random.mock.calls.length).toBeGreaterThan(draws);
    random.mockRestore();
  });

  it('runs either splash when asked for it', () => {
    for (const field of ['light', 'quantum'] as const) {
      const empty = observableState(true);
      const container = document.createElement('div');
      render(() => <EmptyCanvas empty={empty.asReadonly} field={field} />, container);
      vi.advanceTimersByTime(800);
      expect(container.querySelector('canvas')).not.toBeNull();
    }
  });

  it('quiets the canvas while the mark shows, and shows it whenever something is in it', () => {
    const empty = observableState(true);
    const busy = observableState(false);
    const onCanvasShown = vi.fn();
    const container = document.createElement('div');
    render(
      () => (
        <EmptyCanvas empty={empty.asReadonly} busy={busy.asReadonly} splash={false} onCanvasShown={onCanvasShown} />
      ),
      container,
    );
    // Empty at the start: nothing (the host starts it quiet).
    expect(onCanvasShown).not.toHaveBeenCalled();
    // Cleared to load the next at once: shown, never quieted.
    vi.advanceTimersByTime(300);
    empty.set(false);
    // Empty while working (a long load), then loaded: still never quieted.
    empty.set(true);
    busy.set(true);
    vi.advanceTimersByTime(3000);
    empty.set(false);
    busy.set(false);
    expect(onCanvasShown.mock.calls.every(([shown]) => shown)).toBe(true);
    // Empty and idle the whole delay: the mark, the canvas quiet under it; what loads now comes in.
    onCanvasShown.mockClear();
    empty.set(true);
    vi.advanceTimersByTime(800);
    expect(container.querySelector('[role=img]')).not.toBeNull();
    expect(onCanvasShown.mock.calls).toEqual([[false]]);
    empty.set(false);
    expect(onCanvasShown.mock.calls.at(-1)).toEqual([true]);
  });
});

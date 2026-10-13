/** @jest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { EmptyCanvas } from '../empty-canvas.js';

describe('EmptyCanvas', () => {
  beforeEach(() => jest.useFakeTimers());
  afterEach(() => jest.useRealTimers());

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
    jest.advanceTimersByTime(400);
    expect(mark()).toBeNull();
    jest.advanceTimersByTime(400);
    expect(mark()?.dataset.phase).toBe('in');

    empty.set(false);
    expect(mark()?.dataset.phase).toBe('out');
    jest.advanceTimersByTime(300);
    expect(mark()).toBeNull();
  });

  it('never shows when something loads before the delay (clearing a scene to load the next)', () => {
    const { empty, mark } = mount(true);
    jest.advanceTimersByTime(300);
    empty.set(false);
    expect(mark()).toBeNull();
    jest.advanceTimersByTime(1000);
    expect(mark()).toBeNull();
  });
});

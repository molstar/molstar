/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import type { SnapshotTransitionModel } from '../../../model/playback.js';
import { observableState } from '../../../reactive/observable-state.js';
import { mockSnapshots } from '../../../testing/playback.js';
import { PlaybackDeck, SnapshotRow, TransitionRow } from '../playback.js';

// The canvas's playback card: the snapshots' row from two snapshots on, the transition's while the snapshot has one,
// its scrubber seeking and its button playing from the start or stopping.

afterEach(() => {
  document.body.innerHTML = '';
});

function transition(duration: number) {
  const time = observableState(0);
  const playing = observableState(false);
  const length = observableState(duration);
  const toggles: boolean[] = [];
  const model: SnapshotTransitionModel = {
    duration: length.asReadonly,
    step: observableState(100).asReadonly,
    time: time.asReadonly,
    isPlaying: playing.asReadonly,
    seek: (t) => time.set(t),
    togglePlay: () => {
      toggles.push(!playing.value);
      playing.set(!playing.value);
    },
  };
  return { model, time, length, toggles };
}

function mount(count: number, duration: number) {
  const snapshots = mockSnapshots(() => {}, count);
  const t = transition(duration);
  const container = document.createElement('div');
  document.body.append(container);
  render(
    () => (
      <PlaybackDeck>
        <SnapshotRow model={snapshots} />
        <TransitionRow model={t.model} />
      </PlaybackDeck>
    ),
    container,
  );
  const rows = () => [...container.querySelectorAll('[data-row]')];
  return { container, snapshots, t, rows };
}

describe('the playback card', () => {
  it('has the snapshots row from two snapshots on, and the transition row while there is one', () => {
    expect(mount(1, 0).rows()).toHaveLength(0);
    document.body.innerHTML = '';
    const m = mount(3, 0);
    expect(m.rows()).toHaveLength(1);
    expect(m.container.querySelector('[aria-label="Snapshots"] button')?.textContent).toBe('1 / 3');
    m.t.length.set(800);
    expect(m.rows()).toHaveLength(2);
  });

  it('steps the snapshots, seeks the transition and plays it', () => {
    const m = mount(3, 800);
    m.container.querySelector<HTMLButtonElement>('[data-tip="Next snapshot"]')!.click();
    expect(m.snapshots.index.value).toBe(1);
    const scrubber = m.container.querySelector<HTMLInputElement>('[aria-label="Transition time"]')!;
    expect(scrubber.max).toBe('800');
    scrubber.value = '400';
    scrubber.dispatchEvent(new Event('input', { bubbles: true }));
    expect(m.t.time.value).toBe(400);
    expect(m.container.textContent).toContain('0.4 / 0.8 s');
    m.container.querySelector<HTMLButtonElement>('[data-tip="Play the transition"]')!.click();
    expect(m.t.toggles).toEqual([true]);
    expect(m.container.querySelector('[data-tip="Stop the transition"]')).not.toBeNull();
  });
});

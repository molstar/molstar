/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { lostReleases, stubPointerCapture } from '@molstar/prism-ui/testing/gestures';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import type { AudioModel } from '../../../features/snapshots/audio.js';
import { AudioBar } from '../audio-bar.js';

// The audio bar's seek line: a press seeks where it is, a drag seeks once a frame and where it is let go.

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
  vi.restoreAllMocks();
});

describe('the audio bar', () => {
  it('seeks once a frame while dragging along the line, and where it is let go', () => {
    const frames: FrameRequestCallback[] = [];
    vi.spyOn(globalThis, 'requestAnimationFrame').mockImplementation((f) => frames.push(f));
    vi.spyOn(globalThis, 'cancelAnimationFrame').mockImplementation((id) => void (frames[id - 1] = () => {}));
    const seeks: number[] = [];
    const model = {
      time: observableState(0).asReadonly,
      duration: observableState(100).asReadonly,
      playing: observableState(false).asReadonly,
      seek: (f: number) => seeks.push(f),
    } as unknown as AudioModel;
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <AudioBar model={model} />, container);
    const track = container.querySelector<HTMLElement>('[role="slider"]')!;
    track.getBoundingClientRect = () => ({ left: 0, width: 200 }) as DOMRect;
    const fire = (type: string, clientX: number) =>
      track.dispatchEvent(new MouseEvent(type, { bubbles: true, clientX, buttons: type === 'pointerup' ? 0 : 1 }));

    fire('pointerdown', 50);
    expect(seeks).toEqual([0.25]);
    fire('pointermove', 60);
    fire('pointermove', 80);
    fire('pointermove', 100);
    // Nothing yet: the frame has not come.
    expect(seeks).toEqual([0.25]);
    frames.at(-1)!(0);
    expect(seeks).toEqual([0.25, 0.5]);
    fire('pointermove', 120);
    fire('pointerup', 140);
    expect(seeks.at(-1)).toBe(0.7);
    expect(seeks).toHaveLength(3);
  });

  it('a drag that loses its release leaves no listener, and the next one seeks where it is let go', () => {
    const seeks: number[] = [];
    const model = {
      time: observableState(0).asReadonly,
      duration: observableState(100).asReadonly,
      playing: observableState(false).asReadonly,
      seek: (f: number) => seeks.push(f),
    } as unknown as AudioModel;
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <AudioBar model={model} />, container);
    const track = container.querySelector<HTMLElement>('[role="slider"]')!;
    track.getBoundingClientRect = () => ({ left: 0, width: 200 }) as DOMRect;
    const r = lostReleases({
      target: () => track,
      from: { x: 50, y: 0 },
      to: { x: 150, y: 0 },
      reset: () => (seeks.length = 0),
      worked: () => seeks.at(-1) === 0.75,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});

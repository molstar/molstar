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
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { AudioModel } from '../../../features/snapshots/audio.js';
import { AudioBar } from '../audio-bar.js';

// The audio bar's seek line: a press seeks where it is, a drag seeks once a frame and where it is let go.

afterEach(() => {
  document.body.innerHTML = '';
  vi.restoreAllMocks();
});

describe('the audio bar', () => {
  it('seeks once a frame while dragging along the line, and where it is let go', () => {
    const frames: FrameRequestCallback[] = [];
    vi.spyOn(globalThis, 'requestAnimationFrame').mockImplementation((f) => frames.push(f));
    vi.spyOn(globalThis, 'cancelAnimationFrame').mockImplementation((id) => void (frames[id - 1] = () => {}));
    HTMLElement.prototype.setPointerCapture = () => undefined;
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
      track.dispatchEvent(new MouseEvent(type, { bubbles: true, clientX }));

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
});

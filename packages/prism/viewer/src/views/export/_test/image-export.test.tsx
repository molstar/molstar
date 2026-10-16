/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { mockImageExport } from '../../../testing/image-export.js';
import { ImageExport } from '../image-export.js';

// The image export's crop box (v5's screenshot crop): dragged to move, by a handle to resize, set once per frame.

afterEach(() => {
  document.body.innerHTML = '';
  vi.unstubAllGlobals();
});

function setup() {
  const frames: FrameRequestCallback[] = [];
  vi.stubGlobal('requestAnimationFrame', (cb: FrameRequestCallback) => frames.push(cb));
  vi.stubGlobal('cancelAnimationFrame', () => {});
  const host = document.createElement('div');
  document.body.append(host);
  const model = mockImageExport(() => {});
  const setCrop = vi.spyOn(model, 'setCrop');
  render(() => <ImageExport model={model} />, host);
  const crop = host.querySelector('[role="group"]') as HTMLElement;
  // The preview is 200 × 100 px.
  crop.parentElement!.getBoundingClientRect = () => ({ width: 200, height: 100 }) as DOMRect;
  const pointer = (type: string, target: Element, x: number, y: number) =>
    target.dispatchEvent(new MouseEvent(type, { bubbles: true, button: 0, clientX: x, clientY: y }));
  const flush = () => {
    for (const cb of frames.splice(0)) cb(0);
  };
  return { model, setCrop, crop, pointer, flush };
}

describe('the image crop', () => {
  it('moves with a drag, once per frame, and turns auto-crop off', () => {
    const { model, setCrop, crop, pointer, flush } = setup();
    expect(model.autoCrop.value).toBe(true);
    // The automatic box of the mock: 0.17, 0.1, 0.66 × 0.8.
    pointer('pointerdown', crop, 100, 50);
    pointer('pointermove', crop, 110, 55);
    pointer('pointermove', crop, 120, 60);
    expect(setCrop).not.toHaveBeenCalled();
    flush();
    expect(setCrop).toHaveBeenCalledTimes(1);
    const box = setCrop.mock.calls[0][0];
    expect(box.x).toBeCloseTo(0.27);
    // 10 px down is 0.1, as far as the box can go (its height is 0.8).
    expect(box.y).toBeCloseTo(0.2);
    expect(box.width).toBeCloseTo(0.66);
    expect(model.autoCrop.value).toBe(false);
    pointer('pointerup', crop, 120, 60);
    expect(crop.hasAttribute('data-dragging')).toBe(false);
  });

  it('resizes from a corner, kept inside the image', () => {
    const { setCrop, crop, pointer, flush } = setup();
    const se = crop.querySelector('[data-edge="se"]')!;
    pointer('pointerdown', se, 0, 0);
    pointer('pointermove', crop, 1000, 1000);
    flush();
    const box = setCrop.mock.calls[0][0];
    expect(box.x).toBeCloseTo(0.17);
    expect(box.x + box.width).toBeCloseTo(1);
    expect(box.y + box.height).toBeCloseTo(1);
  });

  it('moves with the arrow keys and resizes with ⇧', () => {
    const { model, crop } = setup();
    crop.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    expect(model.crop.value.x).toBeCloseTo(0.18);
    crop.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', shiftKey: true, bubbles: true }));
    expect(model.crop.value.height).toBeCloseTo(0.81);
  });
});

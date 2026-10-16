/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type FramingCamera, frameBelowTop, paddedFocus } from '../framing.js';

/** A camera whose focus is the sphere's center, a radius-proportional distance away, looking down -z, y up. */
const camera = (height: number): FramingCamera => ({
  viewport: { height },
  getFocus: (center, radius) => ({
    target: [center[0], center[1], center[2]],
    position: [center[0], center[1], center[2] + 3 * radius],
    up: [0, 1, 0],
    radius,
  }),
});

describe('the framing below the band', () => {
  it("is Mol*'s own focus without a band", () => {
    const f = paddedFocus(camera(800), [1, 2, 3], 10, 0);
    expect(f.radius).toBe(10);
    expect(f.target).toEqual([1, 2, 3]);
  });

  it('fits the sphere in the canvas below the band and moves it down by half of it', () => {
    const f = paddedFocus(camera(800), [0, 0, 0], 10, 200);
    // The canvas is a quarter shorter for the structure: its radius is a third larger.
    expect(f.radius).toBeCloseTo(10 * (800 / 600), 6);
    // The camera is higher by half the band's px: 100 px of the framed sphere's 2r over 800 px.
    const r = 10 * (800 / 600);
    expect((f.target as number[])[1]).toBeCloseTo((200 * r) / 800, 6);
    expect((f.position as number[])[1]).toBeCloseTo((200 * r) / 800, 6);
    expect((f.target as number[])[0]).toBe(0);
  });

  it('leaves a band that would take the canvas alone', () => {
    expect(paddedFocus(camera(300), [0, 0, 0], 10, 200).radius).toBe(10);
  });

  it('makes resets without a framing of their own frame below the band, and leaves the others alone', () => {
    const calls: { snapshot?: unknown }[] = [];
    const canvas3d = { requestCameraReset: (o?: { snapshot?: unknown }) => void calls.push(o ?? {}) };
    const element = { clientHeight: 500 } as HTMLElement;
    let band = 120;
    const undo = frameBelowTop(canvas3d as never, element, () => band);

    // A reset of its own (a snapshot, a session's camera): as it is.
    const own = { radius: 3 };
    canvas3d.requestCameraReset({ snapshot: own });
    expect(calls[0].snapshot).toBe(own);

    // Reset view: framed below the band, as paddedFocus does (the canvas 500 CSS px, its viewport 1000 px: band × 2).
    canvas3d.requestCameraReset();
    const provider = calls[1].snapshot as (scene: unknown, camera: FramingCamera) => unknown;
    const cam = camera(1000);
    const scene = { boundingSphereVisible: { center: [0, 0, 0], radius: 10 } };
    expect(provider(scene, cam)).toEqual(paddedFocus(camera(1000), [0, 0, 0], 10, 240));
    // No band: Mol*'s own focus (nothing added).
    band = 0;
    expect(provider(scene, cam)).toEqual({});

    undo();
    canvas3d.requestCameraReset();
    expect(calls[2].snapshot).toBeUndefined();
  });
});

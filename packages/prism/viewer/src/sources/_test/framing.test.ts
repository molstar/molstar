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
    // The canvas is a quarter shorter for the structure: its radius is a third larger (and the margin).
    const r = 10 * (800 / 600) * 1.04;
    expect(f.radius).toBeCloseTo(r, 6);
    // The camera is higher by half the band's px: 100 px of the framed sphere's 2r over 800 px.
    expect((f.target as number[])[1]).toBeCloseTo((200 * r) / 800, 6);
    expect((f.position as number[])[1]).toBeCloseTo((200 * r) / 800, 6);
    expect((f.target as number[])[0]).toBe(0);
  });

  it("takes in a perspective's own margin: as close as Mol*'s reset when the margin holds the band, never closer", () => {
    const fov = Math.PI / 4;
    const persp = { ...camera(800), state: { mode: 'perspective' as const, fov } };
    // A band of 3%: Mol*'s own margin (10 / cos(fov / 2) of half-height) holds it, the framing is Mol*'s.
    expect(paddedFocus(persp, [0, 0, 0], 10, 24).radius).toBe(10);
    // A band of 10%: the sphere below it needs 10 / 0.9 · 1.04 of half-height, a radius of that times cos(fov / 2).
    expect(paddedFocus(persp, [0, 0, 0], 10, 80).radius).toBeCloseTo((10 / 0.9) * 1.04 * Math.cos(fov / 2), 6);
    // A band of 25%: farther than Mol*'s, less than the orthographic rule would.
    const big = paddedFocus(persp, [0, 0, 0], 10, 200);
    const r = 10 * (800 / 600) * 1.04 * Math.cos(fov / 2);
    expect(big.radius).toBeCloseTo(r, 6);
    expect((big.target as number[])[1]).toBeCloseTo((200 * r) / (Math.cos(fov / 2) * 800), 6);
  });

  it('on a portrait canvas: the width bounds the sphere, the shift is in the half-height Mol* shows there', () => {
    // 400 × 800 px, orthographic: Mol* shows a half-height of r / (400 / 800) = 2r. A band of 120 px.
    const tall = { ...camera(800), viewport: { height: 800, width: 400 } };
    const f = paddedFocus(tall, [0, 0, 0], 10, 120);
    // The width needs a half-height of 10 · 800 / 400 = 20 (more than the 10 · 800 / 680 below the band), with margin.
    const r = 20 * 1.04 * 0.5;
    expect(f.radius).toBeCloseTo(r, 6);
    // Half the band (60 px) of the 2r / 0.5 world units over 800 px.
    expect((f.target as number[])[1]).toBeCloseTo((120 * r) / (0.5 * 800), 6);
  });

  it('leaves a band that would take the canvas alone', () => {
    expect(paddedFocus(camera(300), [0, 0, 0], 10, 200).radius).toBe(10);
  });

  it('makes resets without a framing of their own frame below the band, and leaves the others alone', () => {
    const calls: { snapshot?: unknown }[] = [];
    const canvas3d = { requestCameraReset: (o?: { snapshot?: unknown }) => void calls.push(o ?? {}) };
    const element = { clientHeight: 500 } as HTMLElement;
    let band = 120;
    const { undo } = frameBelowTop(canvas3d as never, element, () => band);

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
    band = 120;

    undo();
    canvas3d.requestCameraReset();
    expect(calls[2].snapshot).toBeUndefined();
  });

  it("frames Mol*'s first framing of a scene below the band, and again when the band comes, unless the camera moved", () => {
    const calls: { snapshot?: unknown; durationMs?: number }[] = [];
    let n = 0;
    const listeners: ((n: number) => void)[] = [];
    const state = { target: [0, 0, 0], position: [0, 0, 30] };
    const canvas3d = {
      requestCameraReset: (o?: { snapshot?: unknown; durationMs?: number }) => void calls.push(o ?? {}),
      reprCount: {
        get value() {
          return n;
        },
        subscribe: (fn: (n: number) => void) => {
          listeners.push(fn);
          fn(n);
          return { unsubscribe: () => listeners.splice(listeners.indexOf(fn), 1) };
        },
      },
      camera: { state, transition: { inTransition: false } },
    };
    const element = { clientHeight: 1000 } as HTMLElement;
    let band = 0;
    const framing = frameBelowTop(canvas3d as never, element, () => band);
    const resolve = (i: number) => {
      // what Canvas3D does with the reset asked: the snapshot read, the camera put there
      const provider = calls[i].snapshot as (
        scene: unknown,
        camera: FramingCamera,
      ) => { target?: number[]; position?: number[] };
      const cam = camera(1000);
      const own = provider({ boundingSphereVisible: { center: [0, 0, 0], radius: 10 } }, cam);
      const f = own.target ? own : cam.getFocus([0, 0, 0], 10);
      state.target = [...(f.target as number[])];
      state.position = [...(f.position as number[])];
    };

    // the scene from empty to drawn: a reset joins Mol*'s first framing, instantly
    n = 3;
    for (const l of listeners) l(n);
    expect(calls).toHaveLength(1);
    expect(calls[0].durationMs).toBe(0);
    resolve(0);

    // the band comes: framed again, the camera not moved since
    band = 200;
    framing.reframe();
    expect(calls).toHaveLength(2);
    resolve(1);
    expect(state.target[1]).toBeGreaterThan(0);

    // the person turned the camera: the band going leaves it
    state.position = [5, 5, 5];
    band = 0;
    framing.reframe();
    expect(calls).toHaveLength(2);
    framing.undo();
    expect(listeners).toHaveLength(0);
  });

  it("frames again a scene Mol* framed before the hold (or before the band), the camera where Mol*'s framing put it", () => {
    const calls: unknown[] = [];
    const cam = camera(1000);
    const own = cam.getFocus([0, 0, 0], 10);
    const state = { target: [...(own.target as number[])], position: [...(own.position as number[])] };
    const canvas3d = {
      requestCameraReset: (o?: unknown) => void calls.push(o ?? {}),
      boundingSphereVisible: { center: [0, 0, 0], radius: 10 },
      camera: { state, transition: { inTransition: false }, getFocus: cam.getFocus },
    };
    const framing = frameBelowTop(canvas3d as never, { clientHeight: 1000 } as HTMLElement, () => 200);
    framing.reframe();
    expect(calls).toHaveLength(1);
    // turned by the person: left alone
    state.position = [1, 2, 3];
    framing.reframe();
    expect(calls).toHaveLength(1);
    framing.undo();
  });

  it('joins each commit’s reset until one has framed a drawn scene; a camera of its own (a session’s) is left', () => {
    const calls: { snapshot?: unknown; durationMs?: number }[] = [];
    let n = 0;
    const listeners: ((n: number) => void)[] = [];
    const canvas3d = {
      requestCameraReset: (o?: { snapshot?: unknown; durationMs?: number }) => void calls.push(o ?? {}),
      reprCount: {
        get value() {
          return n;
        },
        subscribe: (fn: (n: number) => void) => {
          listeners.push(fn);
          return { unsubscribe: () => listeners.splice(listeners.indexOf(fn), 1) };
        },
      },
    };
    const commit = (count: number) => {
      n = count;
      for (const l of listeners) l(n);
    };
    const framing = frameBelowTop(canvas3d as never, { clientHeight: 1000 } as HTMLElement, () => 0);
    // the first commit's reset met nothing visible: its snapshot never read (Canvas3D frames nothing then)
    commit(1);
    expect(calls).toHaveLength(1);
    // the next commit asks again; this one frames
    commit(4);
    expect(calls).toHaveLength(2);
    (calls[1].snapshot as (scene: unknown, camera: FramingCamera) => unknown)(
      { boundingSphereVisible: { center: [0, 0, 0], radius: 10 } },
      camera(1000),
    );
    commit(5);
    expect(calls).toHaveLength(2);

    // emptied, and a session's camera restored with the next scene: its own, not framed again
    commit(0);
    canvas3d.requestCameraReset({ snapshot: { radius: 3 } });
    commit(3);
    expect(calls).toHaveLength(3);
    framing.undo();
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The camera's default framing with the top of the canvas taken (sequence.md §8): the sequence band covers the
// top `topPx` of the canvas, so the structure is framed in the region below it. Mol*'s reset takes a snapshot provider
// (`requestCameraReset({ snapshot })`); this is one. The user's camera is never touched by it: only a reset (or the
// first framing) asks for it.

import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { Canvas3D } from '@molstar/graphics/canvas3d/canvas3d';

type V3 = [number, number, number] | Float32Array | number[];

/** What the framing reads of Mol*'s camera: its viewport, its projection, and the focus on a sphere (`Camera.getFocus`). */
export interface FramingCamera {
  readonly viewport: { readonly height: number; readonly width?: number };
  /** Its projection: a perspective framing leaves a margin of its own (absent: none, as orthographic). */
  readonly state?: { readonly mode: 'perspective' | 'orthographic'; readonly fov: number };
  getFocus(center: V3, radius: number): { target?: V3; position?: V3; up?: V3 } & Record<string, unknown>;
}

/** Room left around the sphere below the band: a little, where Mol*'s own framing leaves its perspective's margin. */
const Margin = 1.04;

/**
 * The focus on a sphere that fits it in the canvas below the top `topPx`, centered there, as tight as Mol*'s own reset
 * or tighter, never closer. Mol*'s focus on radius r shows a half-height of r / (k·a) (`Camera.targetDistance`): k is
 * cos(fov / 2) in perspective (the sphere's margin, about as much as a band at the default field of view), 1 in
 * orthographic; a is the aspect w / h on a portrait canvas, 1 on a landscape one. Below the band the sphere needs a
 * half-height of r·h / (h − topPx), and a half-width of r (a portrait canvas: h / w of half-height), each with `Margin`;
 * the radius asked is the larger of them times k·a, and the camera goes up by half the band in the half-height shown
 * (the sphere moves down to the middle of what is below it).
 */
export function paddedFocus(camera: FramingCamera, center: V3, radius: number, topPx: number) {
  const h = camera.viewport.height;
  if (!(topPx > 0) || !(h > 2 * topPx)) return camera.getFocus(center, radius);
  const w = camera.viewport.width ?? h;
  const k = camera.state?.mode === 'perspective' ? Math.cos(camera.state.fov / 2) : 1;
  const a = h > w ? w / h : 1;
  // The half-height the sphere needs: below the band, and across the width (a portrait canvas).
  const need = Math.max((radius * h) / (h - topPx), (radius * h) / w) * Margin;
  // Never closer than Mol*'s own framing (its radius also sets the camera's clipping).
  const r = Math.max(radius, need * k * a);
  const focus = camera.getFocus(center, r);
  const { target, position, up } = focus;
  if (!target || !position || !up) return focus;
  // The half-height shown is r / (k·a) over h / 2 px: the camera goes up half the band's px of it.
  const shift = (topPx * r) / (k * a * h);
  for (let i = 0; i < 3; i++) {
    target[i] += up[i] * shift;
    position[i] += up[i] * shift;
  }
  return focus;
}

/** The framing's hold on the plugin's camera resets. */
export interface FrameBelowTop {
  /**
   * Frames the scene again below the band as it is now (the band shown or gone), over `durationMs` (0: at once), if the
   * camera is where the last framing put it: the user's camera is never moved.
   */
  reframe(durationMs?: number): void;
  undo(): void;
}

/** How long the scene glides to its new frame when the band is shown or closed (not at its first coming: counted in). */
export const FrameGlideMs = 250;

/** How far (relative to the framed radius) the camera may be from the last framing and still count as not moved. */
const Untouched = 1e-3;

/**
 * Makes the plugin's camera resets that bring no framing of their own (Reset view, the sequence's reset) frame the scene
 * below the top `topPx()` CSS px of `canvas` (the sequence band; 0: the whole canvas), by `paddedFocus`. A reset with
 * its own snapshot (a snapshot, a session's camera) is left as it is. Mol*'s own first framing of a scene (the scene
 * from empty to drawn) is framed so too, with no change to the graphics package (README decision 28): Canvas3D asks for
 * that reset in its commit and tells `reprCount` after, before the reset is resolved, so a reset asked there joins it
 * (each commit until one has framed a drawn scene).
 */
export function frameBelowTop(canvas3d: Canvas3D, canvas: HTMLElement, topPx: () => number): FrameBelowTop {
  const original = canvas3d.requestCameraReset;
  // Where the last framing put the camera (its target and position), and the radius it framed.
  let framed: { target: number[]; position: number[]; radius: number } | undefined;
  // The scene drawn has had its first framing (ours, or a camera of its own: a session's, a snapshot's)
  let firstDone = false;
  const reset = (options?: Parameters<Canvas3D['requestCameraReset']>[0]) => {
    if (options?.snapshot) {
      framed = undefined;
      firstDone = true;
      return original.call(canvas3d, options);
    }
    original.call(canvas3d, {
      ...options,
      snapshot: (scene, camera) => {
        const { center, radius } = scene.boundingSphereVisible;
        const top = topPx();
        // The band in the camera's pixels (its viewport is the drawing buffer's).
        const scale = camera.viewport.height / (canvas.clientHeight || camera.viewport.height);
        const focus = (
          top > 0 ? paddedFocus(camera, center, radius, top * scale) : camera.getFocus(center, radius)
        ) as Partial<Camera.Snapshot>;
        firstDone = true;
        framed =
          focus.target && focus.position
            ? { target: [...focus.target], position: [...focus.position], radius: Math.max(radius, 1e-6) }
            : undefined;
        return top > 0 ? focus : {};
      },
    });
  };
  canvas3d.requestCameraReset = reset;
  // The scene from empty to drawn: Mol*'s first framing, asked in a commit, takes the framing below the band. Asked
  // again at each commit until a framing has run on a drawn scene: the first commit's reset may meet a scene with
  // nothing visible yet (no framing, its options dropped), and Mol*'s next one would be its own. The scene emptied,
  // the next one is a first framing again.
  const first = canvas3d.reprCount?.subscribe((n) => {
    if (n === 0) {
      framed = undefined;
      firstDone = false;
    } else if (!firstDone) reset({ durationMs: 0 });
  });
  const near = (a: ArrayLike<number>, b: ArrayLike<number>, d: number) =>
    Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) <= d;
  return {
    reframe(durationMs = 0) {
      const camera = canvas3d.camera;
      if (!camera || camera.transition?.inTransition) return;
      const at = (f: { target?: ArrayLike<number>; position?: ArrayLike<number> } | undefined, radius: number) =>
        !!f?.target &&
        !!f.position &&
        near(camera.state.target, f.target, radius * Untouched) &&
        near(camera.state.position, f.position, radius * Untouched);
      // Where the last framing put it, or where Mol*'s own framing of the scene does (one made before this hold, or
      // before the band was known): not moved by the person since.
      const sphere = canvas3d.boundingSphereVisible;
      const own = sphere && sphere.radius > 0 ? camera.getFocus(sphere.center, sphere.radius) : undefined;
      if (!(framed && at(framed, framed.radius)) && !(own && at(own, sphere.radius))) return;
      reset({ durationMs });
    },
    undo() {
      first?.unsubscribe();
      canvas3d.requestCameraReset = original;
    },
  };
}

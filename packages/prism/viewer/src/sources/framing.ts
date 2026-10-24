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

/**
 * Makes the plugin's camera resets that bring no framing of their own (Reset view, the sequence's reset) frame the scene
 * below the top `topPx()` CSS px of `canvas` (the sequence band; 0: the whole canvas), by `paddedFocus`. A reset with
 * its own snapshot (a snapshot, a session's camera) is left as it is; so is Mol*'s own first framing of a scene, which
 * has no hook short of changing the graphics package (README decision 28). Returns the undo.
 */
export function frameBelowTop(canvas3d: Canvas3D, canvas: HTMLElement, topPx: () => number): () => void {
  const original = canvas3d.requestCameraReset;
  canvas3d.requestCameraReset = (options) => {
    if (options?.snapshot) return original.call(canvas3d, options);
    original.call(canvas3d, {
      ...options,
      snapshot: (scene, camera) => {
        const top = topPx();
        if (!(top > 0)) return {};
        // The band in the camera's pixels (its viewport is the drawing buffer's).
        const scale = camera.viewport.height / (canvas.clientHeight || camera.viewport.height);
        const { center, radius } = scene.boundingSphereVisible;
        return paddedFocus(camera, center, radius, top * scale) as Partial<Camera.Snapshot>;
      },
    });
  };
  return () => {
    canvas3d.requestCameraReset = original;
  };
}

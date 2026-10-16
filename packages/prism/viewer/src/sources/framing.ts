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

/** What the framing reads of Mol*'s camera: its viewport, and the focus on a sphere (`Camera.getFocus`). */
export interface FramingCamera {
  readonly viewport: { readonly height: number };
  getFocus(center: V3, radius: number): { target?: V3; position?: V3; up?: V3 } & Record<string, unknown>;
}

/**
 * The focus on a sphere that fits it in the canvas below the top `topPx`: the sphere framed as if the canvas were that
 * much shorter (a larger radius), then moved down by half the band, as the center of the lower region is, by moving the
 * camera up the same amount (a world length per px from the framed radius).
 */
export function paddedFocus(camera: FramingCamera, center: V3, radius: number, topPx: number) {
  const h = camera.viewport.height;
  if (!(topPx > 0) || !(h > 2 * topPx)) return camera.getFocus(center, radius);
  const r = (radius * h) / (h - topPx);
  const focus = camera.getFocus(center, r);
  const { target, position, up } = focus;
  if (!target || !position || !up) return focus;
  // The framed sphere spans the canvas' height, 2r over h px: the camera goes up half the band's px of that.
  const shift = (topPx * r) / h;
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

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The camera's default framing with the top of the canvas taken (prism-sequence.md §4): the sequence band covers the
// top `topPx` of the canvas, so the structure is framed in the region below it. Mol*'s reset takes a snapshot provider
// (`requestCameraReset({ snapshot })`); this is one. The user's camera is never touched by it: only a reset (or the
// first framing) asks for it.

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

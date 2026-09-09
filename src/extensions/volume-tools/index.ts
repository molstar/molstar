/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Interactive volume tools for Mol*: `mask` carves a single soft-edged mask out of a volume,
 * `segmentor` splits one into several bodies. Both are driven by polygons drawn over the
 * viewport, so they share `ViewMask`, its projection and the in-place volume operations.
 *
 * For a full interactive UI, see `src/examples/volume-tools/`.
 */

export { VolumeMaskBehavior } from './mask';
export { VolumeSegmentorBehavior } from './segmentor';
export type { Point2D, ViewMask } from './types';

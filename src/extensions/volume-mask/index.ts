/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Interactive Volume Mask Creator extension for Mol*.
 *
 * Usage
 * -----
 * 1. Add `VolumeMaskBehavior` to `PluginSpec.behaviors`.
 * 2. Create a `VolumeMaskController` (pass the PluginContext).
 * 3. Render `<MaskCreatorPanel plugin={...} controller={...} />` in your UI.
 *
 * The panel renders a drawing overlay canvas into `.msp-viewport` via a React
 * portal whenever drawing mode is active.
 */

export { VolumeMaskBehavior, VolumeMaskController, MASK_OVERLAY_COLOR } from './behavior';
export { MaskVolumeFromSource } from './transformers';
export { MaskCreatorPanel } from './ui/mask-panel';
export { DrawingCanvas } from './ui/drawing-canvas';
export type { ViewMask, MaskCreatorState, MaskSource, Point2D } from './types';

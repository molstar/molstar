/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Volume Mask Creator extension for Mol*
 *
 * For a full interactive UI, see `src/examples/volume-tools/`.
 */

export { VolumeMaskBehavior } from './behavior';
export { MaskVolumeFromSource } from './transformers';
export { MaskSelection } from './selection';
export type { SelectionStore } from './selection';
export { MaskSelectionColorThemeProvider, MaskSelectionColorThemeParams } from './theme';
export type { ViewMask, MaskCreatorState, MaskSource, Point2D } from './types';

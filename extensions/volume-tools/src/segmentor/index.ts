/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Volume Segmentor extension for Mol*: interactively split a volume into several bodies and
 * export one soft-edged MRC mask per body.
 *
 * For a full interactive UI, see `src/examples/volume-tools/`.
 */

export { VolumeSegmentorBehavior } from './behavior.js';
export { VolumeSegmentorManager, isBodyMaskCell } from './manager.js';
export type { VolumeSegmentorState, VolumeSegmentorStats, PreviewMode } from './manager.js';
export { BodyMaskFromLabels, BodyMaskFromLabelsTag } from './transformers.js';
export { BodyLabelColorThemeProvider, BodyLabelColorThemeParams } from './theme.js';
export { BodyLabels } from './labels.js';
export { computeBodyMask, buildCroppedMaskVolume, scatterToFullBox } from './internal/mask-compute.js';
export { computeCandidates } from '@molstar/volume-tools-extension/candidates';
export { flipHandednessInPlace, removeDustInPlace, recomputeStats } from './internal/volume-edit.js';
export { exportBodyMasks, bodyMaskFileName, downloadVolumeMrc, maskBaseName, orderBodiesBySize, writeBodyMaskMrc } from './internal/export.js';
export type { BodyId, BodyInfo, LabelStore, AssignMode, BodyMaskParams, BodyMaskResult, GridBox, ViewMask, Point2D } from './types.js';

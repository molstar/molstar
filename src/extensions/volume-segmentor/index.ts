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

export { VolumeSegmentorBehavior } from './behavior';
export { VolumeSegmentorManager, isBodyMaskCell } from './manager';
export type { VolumeSegmentorState, VolumeSegmentorStats, PreviewMode } from './manager';
export { BodyMaskFromLabels, BodyMaskFromLabelsTag } from './transformers';
export { BodyLabelColorThemeProvider, BodyLabelColorThemeParams } from './theme';
export { BodyLabels } from './labels';
export { computeBodyMask, buildCroppedMaskVolume, scatterToFullBox, softValue } from './internal/mask-compute';
export { computeCandidates } from './internal/candidates';
export { flipHandednessInPlace, removeDustInPlace, recomputeStats } from './internal/volume-edit';
export { exportBodyMasks, bodyMaskFileName, downloadVolumeMrc, maskBaseName, orderBodiesBySize, writeBodyMaskMrc } from './internal/export';
export type { BodyId, BodyInfo, LabelStore, AssignMode, BodyMaskParams, BodyMaskResult, GridBox, ViewMask, Point2D } from './types';

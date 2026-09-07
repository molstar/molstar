/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Volume Bodies extension for Mol*: interactively split a volume into several bodies and
 * export one soft-edged MRC mask per body.
 *
 * For a full interactive UI, see `src/examples/volume-bodies/`.
 */

export { VolumeBodiesBehavior } from './behavior';
export { VolumeBodiesManager, isBodyMaskCell } from './manager';
export type { VolumeBodiesState, VolumeBodiesStats, PreviewMode } from './manager';
export { BodyMaskFromLabels, BodyMaskFromLabelsTag } from './transformers';
export { BodyLabelColorThemeProvider, BodyLabelColorThemeParams } from './theme';
export { BodyLabels } from './labels';
export { computeBodyMask, buildCroppedMaskVolume, scatterToFullBox, softValue } from './internal/mask-compute';
export { computeCandidates } from './internal/candidates';
export { removeDustInPlace, recomputeStats } from './internal/volume-edit';
export { exportBodyMasks, bodyMaskFileName, maskBaseName, orderBodiesBySize, writeBodyMaskMrc } from './internal/export';
export type { BodyId, BodyInfo, LabelStore, AssignMode, BodyMaskParams, BodyMaskResult, GridBox, ViewMask, Point2D } from './types';

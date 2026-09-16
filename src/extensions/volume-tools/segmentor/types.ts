/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Color } from '../../../mol-util/color';
import type { GridBox, ViewMask } from '../types';

export type { ViewMask, Point2D, GridBox } from '../types';

/** Body identifier, 1..255. Ids are never reused within a session; 0 means "unassigned". */
export type BodyId = number;

export const MaxBodyId = 255;

/**
 * Definition of one body. Its voxels are the candidates (above threshold) whose projection
 * lies inside every view polygon (inverted views exclude instead); a remainder body takes
 * every candidate no other body claimed. Bodies are resolved in list order, so where polygons
 * overlap the body listed first wins.
 */
export interface BodyInfo {
    id: BodyId
    name: string
    color: Color
    /** View polygons (with camera snapshots) intersected to define the body. */
    views: ViewMask[]
    /** Takes all voxels above the threshold not claimed by other bodies. */
    remainder: boolean
    /** Number of voxels currently labelled with this body. */
    voxelCount: number
    /** Per-body override of the global extend (voxels). */
    extend?: number
    /** Per-body override of the global soft-edge width (voxels). */
    softEdge?: number
}

/**
 * Voxel labels of one volume, derived from the body definitions. Attached to the source
 * `Volume` as an ad-hoc property (see `BodyLabels`), so the color theme and the mask
 * transformer can read it without pushing large arrays through state-tree params.
 */
export interface LabelStore {
    /** One label per voxel in the memory order of `volume.grid.cells.data`; 0 = unassigned. */
    labels: Uint8Array
    /** Bumped after every change; mirrored into theme/transformer params to trigger updates. */
    version: number
    /** Body definitions in priority order. */
    bodies: BodyInfo[]
    nextId: BodyId
}

export type AssignMode = 'replace' | 'unassigned-only' | 'erase';

export interface BodyMaskParams {
    /** Dilate the binary body by this many voxels before the soft edge. */
    extend: number
    /** Width of the raised-cosine soft edge in voxels. */
    softEdge: number
    /** Exclude labelled voxels whose density is below the current threshold. */
    pruneBelowThreshold: boolean
}

/** Axis-aligned box in grid (voxel index) space; `min` inclusive, `dims` voxel counts. */
export interface BodyMaskResult {
    box: GridBox
    /** Mask values in [0, 1], canonical local order (x fastest) within `box`. */
    data: Float32Array
    /** Number of binary body voxels the mask was derived from. */
    voxelCount: number
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { Volume } from '../../mol-model/volume';
import { StateTransform } from '../../mol-state';

export type Point2D = [number, number];

/**
 * A 2D polygon drawn on top of the viewport at a specific camera orientation.
 * Everything needed to project any voxel into this view is captured at draw time.
 */
export interface ViewMask {
    id: string;
    label: string;
    /** Polygon vertices in CSS pixel coords on the overlay canvas (y=0 at top). */
    polygon: Point2D[];
    /** CSS pixel dimensions of the overlay canvas at draw time. */
    canvasWidth: number;
    canvasHeight: number;
    /** Physical pixel dimensions of the WebGL canvas at draw time (may differ by dpr). */
    viewportWidth: number;
    viewportHeight: number;
    /** Full camera state frozen at draw time. */
    cameraSnapshot: Camera.Snapshot;
    /** When true, voxels OUTSIDE this polygon are selected instead of inside. */
    inverted?: boolean;
}

export type MaskSource = 'volume' | 'structure';

export interface MaskCreatorState {
    isDrawing: boolean;
    targetVolumeRef: StateTransform.Ref | undefined;
    targetStructureRef: StateTransform.Ref | undefined;
    /** All available auth chain IDs in the current target structure. */
    availableChainIds: string[];
    viewMasks: ViewMask[];
    threshold: Volume.IsoValue;
    dilation: number;
    softEdge: number;
    maskSource: MaskSource;
    proteinRadius: number;
    maskVolumeRef: StateTransform.Ref | undefined;
    /** Float32 mask data (0–1) kept for MRC export; has soft edge values when softEdge > 0. */
    maskData: Float32Array | undefined;
    /** When true, the displayed/exported mask is inverted (1 − original). */
    maskInverted: boolean;
    /** Opacity (0–1) applied to the target volume's representations. */
    volumeOpacity: number;
    /** Auth chain IDs (e.g. 'A', 'B') used for structure-based masks. Empty = no chains selected. */
    selectedChainIds: string[];
}

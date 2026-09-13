/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { StateTransform } from '../../../mol-state';
import type { ViewMask } from '../types';

export type { Point2D, ViewMask } from '../types';

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

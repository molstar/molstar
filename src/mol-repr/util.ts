/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { defaults } from '../mol-util';
import { Structure } from '../mol-model/structure';
import { VisualQuality } from '../mol-geo/geometry/base';
import { Box3D, SpacegroupCell } from '../mol-math/geometry';
import { ModelSymmetry } from '../mol-model-formats/structure/property/symmetry';
import { Volume } from '../mol-model/volume';
import { Location } from '../mol-model/location';

export interface VisualUpdateState {
    updateTransform: boolean
    updateMatrix: boolean
    updateColor: boolean
    updateSize: boolean
    createGeometry: boolean
    createNew: boolean

    /** holds contextual info, is not reset  */
    info: { [k: string]: unknown }
}
export namespace VisualUpdateState {
    export function create(): VisualUpdateState {
        return {
            updateTransform: false,
            updateMatrix: false,
            updateColor: false,
            updateSize: false,
            createGeometry: false,
            createNew: false,

            info: {}
        };
    }
    export function reset(state: VisualUpdateState) {
        state.updateTransform = false;
        state.updateMatrix = false;
        state.updateColor = false;
        state.updateSize = false;
        state.createGeometry = false;
        state.createNew = false;
    }
}

export type LocationCallback = (loc: Location, isSecondary: boolean) => void

//

export interface QualityProps {
    quality: VisualQuality
    detail: number
    radialSegments: number
    linearSegments: number
    resolution: number
    probePositions: number
    doubleSided: boolean
    xrayShaded: boolean | 'inverted'
    alpha: number
    transparentBackfaces: 'off' | 'on' | 'opaque'
}

export const DefaultQualityThresholds = {
    lowestElementCount: 1_000_000,
    lowerElementCount: 500_000,
    lowElementCount: 100_000,
    mediumElementCount: 20_000,
    highElementCount: 2_000,
    coarseGrainedFactor: 10,

    elementCountFactor: 1
};
export type QualityThresholds = typeof DefaultQualityThresholds

export function getStructureQuality(structure: Structure, tresholds: Partial<QualityThresholds> = {}): VisualQuality {
    const t = { ...DefaultQualityThresholds, ...tresholds };
    let score = structure.elementCount * t.elementCountFactor;
    if (structure.isCoarseGrained || structure.isCoarse) score *= t.coarseGrainedFactor;
    if (score > t.lowestElementCount) {
        return 'lowest';
    } else if (score > t.lowerElementCount) {
        return 'lower';
    } else if (score > t.lowElementCount) {
        return 'low';
    } else if (score > t.mediumElementCount) {
        return 'medium';
    } else if (score > t.highElementCount) {
        return 'high';
    } else {
        return 'higher';
    }
}

/**
 * Uses cell volume to avoid costly boundary calculation if
 * - single model
 * - non-empty 'P 1' spacegroup
 */
function getRootVolume(structure: Structure) {
    if (structure.root.models.length === 1) {
        const sym = ModelSymmetry.Provider.get(structure.root.model);
        if (sym && sym.spacegroup.name === 'P 1' && !SpacegroupCell.isZero(sym.spacegroup.cell)) {
            return sym.spacegroup.cell.volume;
        }
    }
    return Box3D.volume(structure.root.boundary.box);
}

export function getQualityProps(props: Partial<QualityProps>, data?: any) {
    let quality = defaults(props.quality, 'auto' as VisualQuality);
    let detail = defaults(props.detail, 1);
    let radialSegments = defaults(props.radialSegments, 12);
    let linearSegments = defaults(props.linearSegments, 8);
    let resolution = defaults(props.resolution, 2);
    let probePositions = defaults(props.probePositions, 12);
    let doubleSided = defaults(props.doubleSided, true);

    let volume = 0;
    if (quality === 'auto') {
        if (data instanceof Structure) {
            quality = getStructureQuality(data.root);
            volume = getRootVolume(data);
        } else if (Volume.is(data)) {
            const [x, y, z] = data.grid.cells.space.dimensions;
            volume = x * y * z;
            quality = volume < 10_000_000 ? 'medium' : 'low';
        }
    }

    switch (quality) {
        case 'highest':
            detail = 3;
            radialSegments = 36;
            linearSegments = 18;
            resolution = 0.1;
            probePositions = 72;
            doubleSided = true;
            break;
        case 'higher':
            detail = 3;
            radialSegments = 28;
            linearSegments = 14;
            resolution = 0.3;
            probePositions = 48;
            doubleSided = true;
            break;
        case 'high':
            detail = 2;
            radialSegments = 20;
            linearSegments = 10;
            resolution = 0.5;
            probePositions = 36;
            doubleSided = true;
            break;
        case 'medium':
            detail = 1;
            radialSegments = 12;
            linearSegments = 8;
            resolution = 0.8;
            probePositions = 24;
            doubleSided = true;
            break;
        case 'low':
            detail = 0;
            radialSegments = 8;
            linearSegments = 3;
            resolution = 1.3;
            probePositions = 24;
            doubleSided = false;
            break;
        case 'lower':
            detail = 0;
            radialSegments = 4;
            linearSegments = 2;
            resolution = 3;
            probePositions = 12;
            doubleSided = false;
            break;
        case 'lowest':
            detail = 0;
            radialSegments = 2;
            linearSegments = 1;
            resolution = 8;
            probePositions = 12;
            doubleSided = false;
            break;
        case 'custom':
            // use defaults or given props as set above
            break;
    }

    // max resolution based on volume (for 'auto' quality)
    resolution = Math.max(resolution, volume / 500_000_000);
    resolution = Math.min(resolution, 20);

    if (props.transparentBackfaces === 'off' && ((props.alpha !== undefined && props.alpha < 1) || !!props.xrayShaded)) {
        doubleSided = false;
    }

    return {
        detail,
        radialSegments,
        linearSegments,
        resolution,
        probePositions,
        doubleSided
    };
}
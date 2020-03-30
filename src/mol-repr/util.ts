/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { defaults } from '../mol-util';
import { Structure } from '../mol-model/structure';
import { VisualQuality } from '../mol-geo/geometry/base';

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
        }
    }
    export function reset(state: VisualUpdateState) {
        state.updateTransform = false
        state.updateMatrix = false
        state.updateColor = false
        state.updateSize = false
        state.createGeometry = false
        state.createNew = false
    }
}

//

export interface QualityProps {
    quality: VisualQuality
    detail: number
    radialSegments: number
    linearSegments: number
    resolution: number
    doubleSided: boolean
}

export const DefaultQualityThresholds = {
    lowestElementCount: 500_000,
    lowerElementCount: 200_000,
    lowElementCount: 60_000,
    mediumElementCount: 20_000,
    highElementCount: 2_000,
    coarseGrainedFactor: 10,

    elementCountFactor: 1
}
export type QualityThresholds = typeof DefaultQualityThresholds

export function getStructureQuality(structure: Structure, tresholds: Partial<QualityThresholds> = {}): VisualQuality {
    const t = { ...DefaultQualityThresholds, ...tresholds }
    let score = structure.elementCount * t.elementCountFactor
    if (structure.isCoarseGrained) score *= t.coarseGrainedFactor
    if (score > t.lowestElementCount) {
        return 'lowest'
    } else if (score > t.lowerElementCount) {
        return 'lower'
    } else if (score > t.lowElementCount) {
        return 'low'
    } else if (score > t.mediumElementCount) {
        return 'medium'
    } else if (score > t.highElementCount) {
        return 'high'
    } else {
        return 'higher'
    }
}

export function getQualityProps(props: Partial<QualityProps>, data?: any) {
    let quality = defaults(props.quality, 'auto' as VisualQuality)
    let detail = defaults(props.detail, 1)
    let radialSegments = defaults(props.radialSegments, 12)
    let linearSegments = defaults(props.linearSegments, 8)
    let resolution = defaults(props.resolution, 2)
    let doubleSided = defaults(props.doubleSided, true)

    if (quality === 'auto' && data instanceof Structure) {
        quality = getStructureQuality(data.root)
    }

    switch (quality) {
        case 'highest':
            detail = 3
            radialSegments = 36
            linearSegments = 18
            resolution = 0.1
            doubleSided = true
            break
        case 'higher':
            detail = 3
            radialSegments = 28
            linearSegments = 14
            resolution = 0.3
            doubleSided = true
            break
        case 'high':
            detail = 2
            radialSegments = 20
            linearSegments = 10
            resolution = 0.5
            doubleSided = true
            break
        case 'medium':
            detail = 1
            radialSegments = 12
            linearSegments = 8
            resolution = 1
            doubleSided = true
            break
        case 'low':
            detail = 0
            radialSegments = 8
            linearSegments = 3
            resolution = 2
            doubleSided = false
            break
        case 'lower':
            detail = 0
            radialSegments = 4
            linearSegments = 2
            resolution = 4
            doubleSided = false
            break
        case 'lowest':
            detail = 0
            radialSegments = 2
            linearSegments = 1
            resolution = 8
            doubleSided = false
            break
        case 'custom':
            // use defaults or given props as set above
            break
    }

    return {
        detail,
        radialSegments,
        linearSegments,
        resolution,
        doubleSided
    }
}
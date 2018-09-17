/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { defaults } from 'mol-util';
import { Structure } from 'mol-model/structure';
import { VisualQuality } from '../geometry/geometry';

export interface QualityProps {
    quality: VisualQuality
    detail: number
    radialSegments: number
    linearSegments: number
}

export function getQualityProps(props: Partial<QualityProps>, structure?: Structure) {
    let quality = defaults(props.quality, 'auto' as VisualQuality)
    let detail = defaults(props.detail, 1)
    let radialSegments = defaults(props.radialSegments, 12)
    let linearSegments = defaults(props.linearSegments, 8)

    if (quality === 'auto' && structure) {
        let score = structure.elementCount
        if (structure.isCoarse) score *= 10
        if (score > 500_000) {
            quality = 'lowest'
        } else if (score > 100_000) {
            quality = 'low'
        } else if (score > 30_000) {
            quality = 'medium'
        } else {
            quality = 'high'
        }
    }

    switch (quality) {
        case 'highest':
            detail = 3
            radialSegments = 36
            linearSegments = 18
            break
        case 'high':
            detail = 2
            radialSegments = 24
            linearSegments = 12
            break
        case 'medium':
            detail = 1
            radialSegments = 12
            linearSegments = 8
            break
        case 'low':
            detail = 0
            radialSegments = 8
            linearSegments = 3
            break
        case 'lowest':
            detail = 0
            radialSegments = 4
            linearSegments = 2
            break
        case 'custom':
            // use defaults or given props as set above
            break
    }

    return {
        detail,
        radialSegments,
        linearSegments
    }
}
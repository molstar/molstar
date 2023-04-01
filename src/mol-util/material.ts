/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from './type-helpers';
import { ParamDefinition as PD } from './param-definition';

export interface Material {
    /** Normalized to [0, 1] range */
    metalness: number,
    /** Normalized to [0, 1] range */
    roughness: number
    /** Normalized to [0, 1] range */
    bumpiness: number
}

export function Material(values?: Partial<Material>) {
    return { ...Material.Zero, ...values };
}

export namespace Material {
    export const Zero: Material = { metalness: 0, roughness: 0, bumpiness: 0 };

    export function toArray<T extends NumberArray>(material: Material, array: T, offset: number) {
        array[offset] = material.metalness * 255;
        array[offset + 1] = material.roughness * 255;
        array[offset + 2] = material.bumpiness * 255;
        return array;
    }

    export function toString({ metalness, roughness, bumpiness }: Material) {
        return `M ${metalness.toFixed(2)} | R ${roughness.toFixed(2)} | B ${bumpiness.toFixed(2)}`;
    }

    export function getParam(info?: { isExpanded?: boolean, isFlat?: boolean }) {
        return PD.Group({
            metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
            roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
            bumpiness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        }, {
            ...info,
            presets: [
                [{ metalness: 0, roughness: 1, bumpiness: 0 }, 'Matte'],
                [{ metalness: 0, roughness: 0.2, bumpiness: 0 }, 'Plastic'],
                [{ metalness: 0, roughness: 0.6, bumpiness: 0 }, 'Glossy'],
                [{ metalness: 1.0, roughness: 0.6, bumpiness: 0 }, 'Metallic'],
            ]
        });
    }
}

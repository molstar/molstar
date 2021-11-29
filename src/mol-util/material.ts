/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from './type-helpers';
import { ParamDefinition as PD } from './param-definition';
import { toFixed } from './number';

/** Material properties expressed as a single number */
export type Material = { readonly '@type': 'material' } & number

export function Material(hex: number) { return hex as Material; }

export namespace Material {
    export function fromNormalized(metalness: number, roughness: number): Material {
        return (((metalness * 255) << 16) | ((roughness * 255) << 8)) as Material;
    }

    export function fromObjectNormalized(v: { metalness: number, roughness: number }): Material {
        return fromNormalized(v.metalness, v.roughness);
    }

    export function toObjectNormalized(material: Material, fractionDigits?: number) {
        const metalness = (material >> 16 & 255) / 255;
        const roughness = (material >> 8 & 255) / 255;
        return {
            metalness: fractionDigits ? toFixed(metalness, fractionDigits) : metalness,
            roughness: fractionDigits ? toFixed(roughness, fractionDigits) : roughness
        };
    }

    export function toArray(material: Material, array: NumberArray, offset: number) {
        array[offset] = (material >> 16 & 255);
        array[offset + 1] = (material >> 8 & 255);
        return array;
    }

    export function toString(material: Material) {
        const metalness = (material >> 16 & 255) / 255;
        const roughness = (material >> 8 & 255) / 255;
        return `M ${metalness} | R ${roughness}`;
    }

    export function getParam(info?: { isExpanded?: boolean, isFlat?: boolean }) {
        return PD.Converted(
            (v: Material) => toObjectNormalized(v, 2),
            (v: { metalness: number, roughness: number }) => fromObjectNormalized(v),
            PD.Group({
                metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
                roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
            }, {
                ...info,
                presets: [
                    [{ metalness: 0, roughness: 1 }, 'Matte'],
                    [{ metalness: 0.5, roughness: 0.5 }, 'Metallic'],
                    [{ metalness: 0, roughness: 0.25 }, 'Plastic'],
                ]
            })
        );
    }
}

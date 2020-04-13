/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../mol-util/type-helpers';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Hcl } from './spaces/hcl';
import { Lab } from './spaces/lab';

/** RGB color triplet expressed as a single number */
export type Color = { readonly '@type': 'color' } & number

export function Color(hex: number) { return hex as Color; }

export namespace Color {
    export function toStyle(hexColor: Color) {
        return `rgb(${hexColor >> 16 & 255}, ${hexColor >> 8 & 255}, ${hexColor & 255})`;
    }

    export function toHexString(hexColor: Color) {
        return '0x' + ('000000' + hexColor.toString(16)).slice(-6);
    }

    export function toRgbString(hexColor: Color) {
        return `RGB: ${Color.toRgb(hexColor).join(', ')}`;
    }

    export function toRgb(hexColor: Color): [number, number, number] {
        return [ hexColor >> 16 & 255, hexColor >> 8 & 255, hexColor & 255 ];
    }

    export function toRgbNormalized(hexColor: Color): [number, number, number] {
        return [ (hexColor >> 16 & 255) / 255, (hexColor >> 8 & 255) / 255, (hexColor & 255) / 255 ];
    }

    export function fromRgb(r: number, g: number, b: number): Color {
        return ((r << 16) | (g << 8) | b) as Color;
    }

    export function fromNormalizedRgb(r: number, g: number, b: number): Color {
        return (((r * 255) << 16) | ((g * 255) << 8) | (b * 255)) as Color;
    }

    export function fromArray(array: NumberArray, offset: number): Color {
        return fromRgb(array[offset], array[offset + 1], array[offset + 2]);
    }

    export function fromNormalizedArray(array: NumberArray, offset: number): Color {
        return fromNormalizedRgb(array[offset], array[offset + 1], array[offset + 2]);
    }

    /** Copies hex color to rgb array */
    export function toArray(hexColor: Color, array: NumberArray, offset: number) {
        array[ offset ] = (hexColor >> 16 & 255);
        array[ offset + 1 ] = (hexColor >> 8 & 255);
        array[ offset + 2 ] = (hexColor & 255);
        return array;
    }

    /** Copies normalized (0 to 1) hex color to rgb array */
    export function toArrayNormalized<T extends NumberArray>(hexColor: Color, array: T, offset: number) {
        array[ offset ] = (hexColor >> 16 & 255) / 255;
        array[ offset + 1 ] = (hexColor >> 8 & 255) / 255;
        array[ offset + 2 ] = (hexColor & 255) / 255;
        return array;
    }

    /** Copies hex color to rgb vec3 */
    export function toVec3(out: Vec3, hexColor: Color) {
        out[0] = (hexColor >> 16 & 255);
        out[1] = (hexColor >> 8 & 255);
        out[2] = (hexColor & 255);
        return out;
    }

    /** Copies normalized (0 to 1) hex color to rgb vec3 */
    export function toVec3Normalized(out: Vec3, hexColor: Color) {
        out[0] = (hexColor >> 16 & 255) / 255;
        out[1] = (hexColor >> 8 & 255) / 255;
        out[2] = (hexColor & 255) / 255;
        return out;
    }

    /** Linear interpolation between two colors in rgb space */
    export function interpolate(c1: Color, c2: Color, t: number): Color {
        const r1 = c1 >> 16 & 255;
        const g1 = c1 >> 8 & 255;
        const b1 = c1 & 255;
        const r2 = c2 >> 16 & 255;
        const g2 = c2 >> 8 & 255;
        const b2 = c2 & 255;

        const r = r1 + (r2 - r1) * t;
        const g = g1 + (g2 - g1) * t;
        const b = b1 + (b2 - b1) * t;

        return ((r << 16) | (g << 8) | b) as Color;
    }

    const tmpSaturateHcl = [0, 0, 0] as Hcl;
    export function saturate(c: Color, amount: number): Color {
        Hcl.fromColor(tmpSaturateHcl, c);
        return Hcl.toColor(Hcl.saturate(tmpSaturateHcl, tmpSaturateHcl, amount));
    }

    export function desaturate(c: Color, amount: number): Color {
        return saturate(c, -amount);
    }

    const tmpDarkenLab = [0, 0, 0] as Lab;
    export function darken(c: Color, amount: number): Color {
        Lab.fromColor(tmpDarkenLab, c);
        return Lab.toColor(Lab.darken(tmpDarkenLab, tmpDarkenLab, amount));
    }

    export function lighten(c: Color, amount: number): Color {
        return darken(c, -amount);
    }
}

export interface ColorList {
    label: string
    description: string
    list: Color[]
    type: 'sequential' | 'diverging' | 'qualitative'
}
export function ColorList(label: string, type: 'sequential' | 'diverging' | 'qualitative', description: string, list: number[]): ColorList {
    return { label, description, list: list as Color[], type };
}

export type ColorTable<T extends { [k: string]: number[] }> = { [k in keyof T]: Color[] }
export function ColorTable<T extends { [k: string]: number[] }>(o: T) { return o as unknown as ColorTable<T>; }

export type ColorMap<T extends { [k: string]: number }> = { [k in keyof T]: Color }
export function ColorMap<T extends { [k: string]: number }>(o: T) { return o as unknown as ColorMap<T>; }
export function getAdjustedColorMap<T extends { [k: string]: number }>(map: ColorMap<T>, saturation: number, lightness: number) {
    const adjustedMap: { [k: string]: Color } = {};
    for (const e in map) {
        let c = map[e];
        c = Color.saturate(c, saturation);
        c = Color.darken(c, -lightness);
        adjustedMap[e] = c;
    }
    return adjustedMap as ColorMap<T>;
}

export type ColorSwatch = [string, Color][]
export function ColorSwatch(l: [string, number][]) { return l as unknown as ColorSwatch; }
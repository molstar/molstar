/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../mol-util/type-helpers';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Hcl } from './spaces/hcl';
import { Lab } from './spaces/lab';
import { Hsl } from './spaces/hsl';

/** RGB color triplet expressed as a single number */
export type Color = { readonly '@type': 'color' } & number

export function Color(hex: number) { return hex as Color; }

export namespace Color {
    export function toStyle(hexColor: Color): string {
        return `rgb(${hexColor >> 16 & 255}, ${hexColor >> 8 & 255}, ${hexColor & 255})`;
    }

    export function toHexStyle(hexColor: Color): string {
        return '#' + ('000000' + hexColor.toString(16)).slice(-6);
    }

    export function toHexString(hexColor: Color): string {
        return '0x' + ('000000' + hexColor.toString(16)).slice(-6);
    }

    export function toRgbString(hexColor: Color): string {
        return `RGB: ${Color.toRgb(hexColor).join(', ')}`;
    }

    export function toRgb(hexColor: Color): [number, number, number] {
        return [hexColor >> 16 & 255, hexColor >> 8 & 255, hexColor & 255];
    }

    export function toRgbNormalized(hexColor: Color): [number, number, number] {
        return [(hexColor >> 16 & 255) / 255, (hexColor >> 8 & 255) / 255, (hexColor & 255) / 255];
    }

    export function fromHexStyle(s: string): Color {
        return parseInt(s.replace('#', '0x')) as Color;
    }

    export function fromHexString(s: string): Color {
        return parseInt(s) as Color;
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

    export function fromColorListEntry(entry: ColorListEntry): Color {
        if (typeof entry === 'number') return entry;
        else return entry[0];
    }

    /** Copies hex color to rgb array */
    export function toArray<T extends NumberArray>(hexColor: Color, array: T, offset: number): T {
        array[offset] = (hexColor >> 16 & 255);
        array[offset + 1] = (hexColor >> 8 & 255);
        array[offset + 2] = (hexColor & 255);
        return array;
    }

    /** Copies normalized (0 to 1) hex color to rgb array */
    export function toArrayNormalized<T extends NumberArray>(hexColor: Color, array: T, offset: number): T {
        array[offset] = (hexColor >> 16 & 255) / 255;
        array[offset + 1] = (hexColor >> 8 & 255) / 255;
        array[offset + 2] = (hexColor & 255) / 255;
        return array;
    }

    /** Copies hex color to rgb vec3 */
    export function toVec3(out: Vec3, hexColor: Color): Vec3 {
        out[0] = (hexColor >> 16 & 255);
        out[1] = (hexColor >> 8 & 255);
        out[2] = (hexColor & 255);
        return out;
    }

    /** Copies normalized (0 to 1) hex color to rgb vec3 */
    export function toVec3Normalized(out: Vec3, hexColor: Color): Vec3 {
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

    const _interpolateHsl1 = Hsl.zero();
    const _interpolateHsl2 = Hsl.zero();

    /** Linear interpolation between two colors in HSL space */
    export function interpolateHsl(c1: Color, c2: Color, t: number): Color {
        const hsl1 = Hsl.fromColor(_interpolateHsl1, c1);
        const hsl2 = Hsl.fromColor(_interpolateHsl2, c2);
        Hsl.interpolate(hsl1, hsl1, hsl2, t);
        return Hsl.toColor(hsl1);
    }

    export function hasHue(c: Color): boolean {
        const r = c >> 16 & 255;
        const g = c >> 8 & 255;
        const b = c & 255;
        return r !== g || r !== b;
    }

    const tmpSaturateHcl = [0, 0, 0] as unknown as Hcl;
    export function saturate(c: Color, amount: number): Color {
        if (!hasHue(c)) return c;
        Hcl.fromColor(tmpSaturateHcl, c);
        return Hcl.toColor(Hcl.saturate(tmpSaturateHcl, tmpSaturateHcl, amount));
    }

    export function desaturate(c: Color, amount: number): Color {
        return saturate(c, -amount);
    }

    const tmpDarkenLab = [0, 0, 0] as unknown as Lab;
    export function darken(c: Color, amount: number): Color {
        Lab.fromColor(tmpDarkenLab, c);
        return Lab.toColor(Lab.darken(tmpDarkenLab, tmpDarkenLab, amount));
    }

    export function lighten(c: Color, amount: number): Color {
        return darken(c, -amount);
    }

    function _luminance(x: number): number {
        return x <= 0.03928 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4);
    }

    /**
     * Relative luminance
     * http://www.w3.org/TR/2008/REC-WCAG20-20081211/#relativeluminancedef
     */
    export function luminance(c: Color): number {
        const r = _luminance((c >> 16 & 255) / 255);
        const g = _luminance((c >> 8 & 255) / 255);
        const b = _luminance((c & 255) / 255);
        return 0.2126 * r + 0.7152 * g + 0.0722 * b;
    }

    /**
     * WCAG contrast ratio
     * http://www.w3.org/TR/2008/REC-WCAG20-20081211/#contrast-ratiodef
     */
    export function contrast(a: Color, b: Color): number {
        const l1 = luminance(a);
        const l2 = luminance(b);
        return l1 > l2 ? (l1 + 0.05) / (l2 + 0.05) : (l2 + 0.05) / (l1 + 0.05);
    };

    //

    function _sRGBToLinear(c: number): number {
        return (c < 0.04045) ? c * 0.0773993808 : Math.pow(c * 0.9478672986 + 0.0521327014, 2.4);
    }

    export function sRGBToLinear(c: Color): Color {
        return fromNormalizedRgb(
            _sRGBToLinear((c >> 16 & 255) / 255),
            _sRGBToLinear((c >> 8 & 255) / 255),
            _sRGBToLinear((c & 255) / 255)
        );
    }

    function _linearToSRGB(c: number): number {
        return (c < 0.0031308) ? c * 12.92 : 1.055 * (Math.pow(c, 0.41666)) - 0.055;
    }

    export function linearToSRGB(c: Color): Color {
        return fromNormalizedRgb(
            _linearToSRGB((c >> 16 & 255) / 255),
            _linearToSRGB((c >> 8 & 255) / 255),
            _linearToSRGB((c & 255) / 255)
        );
    }
}

export type ColorListEntry = Color | [Color, number /** normalized value from 0 to 1 */]

export interface ColorList {
    label: string
    description: string
    list: ColorListEntry[]
    type: 'sequential' | 'diverging' | 'cyclical' | 'qualitative'
}
export function ColorList(label: string, type: ColorList['type'], description: string, list: (number | [number, number])[]): ColorList {
    return { label, description, list: list as ColorListEntry[], type };
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
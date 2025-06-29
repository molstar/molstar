/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color, ColorListEntry } from './color';
import { ColorNames } from './names';

const hexColorRegex = /^#([0-9A-F]{3}){1,2}$/i;
const rgbColorRegex = /^rgb\(\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*\)$/i;

export function decodeColor(colorString: string | undefined | null): Color | undefined {
    if (colorString === undefined || colorString === null) return undefined;
    let result: Color | undefined;
    if (hexColorRegex.test(colorString)) {
        if (colorString.length === 4) {
            // convert short form to full form (#f0f -> #ff00ff)
            colorString = `#${colorString[1]}${colorString[1]}${colorString[2]}${colorString[2]}${colorString[3]}${colorString[3]}`;
        }
        result = Color.fromHexStyle(colorString);
        if (result !== undefined && !isNaN(result)) return result;
    }

    result = ColorNames[colorString.toLowerCase() as keyof typeof ColorNames];
    if (result !== undefined) return result;

    const rgbMatch = rgbColorRegex.exec(colorString);
    if (rgbMatch) {
        const r = parseInt(rgbMatch[1], 10);
        const g = parseInt(rgbMatch[2], 10);
        const b = parseInt(rgbMatch[3], 10);
        if (r >= 0 && r <= 255 && g >= 0 && g <= 255 && b >= 0 && b <= 255) {
            return Color.fromRgb(r, g, b);
        }
    }

    return undefined;
}

export function getColorGradientBanded(colors: ColorListEntry[]) {
    const n = colors.length;
    const styles: string[] = [];

    const hasOffsets = colors.every(c => Array.isArray(c));
    if (hasOffsets) {
        const off = [...colors] as [Color, number][];
        // 0 colors present
        if (!off[0]) {
            return 'linear-gradient(to right, #000 0%, #000 100%)';
        }
        off.sort((a, b) => a[1] - b[1]);
        styles.push(`${Color.toStyle(off[0][0])} ${(100 * off[0][1]).toFixed(2)}%`);
        for (let i = 0, il = off.length - 1; i < il; ++i) {
            const [c0, o0] = off[i];
            const [c1, o1] = off[i + 1];
            const o = o0 + (o1 - o0) / 2;
            styles.push(
                `${Color.toStyle(c0)} ${(100 * o).toFixed(2)}%`,
                `${Color.toStyle(c1)} ${(100 * o).toFixed(2)}%`
            );
        }
        styles.push(`${Color.toStyle(off[off.length - 1][0])} ${(100 * off[off.length - 1][1]).toFixed(2)}%`);
    } else {
        styles.push(`${colorEntryToStyle(colors[0])} ${100 * (1 / n)}%`);
        for (let i = 1, il = n - 1; i < il; ++i) {
            styles.push(
                `${colorEntryToStyle(colors[i])} ${100 * (i / n)}%`,
                `${colorEntryToStyle(colors[i])} ${100 * ((i + 1) / n)}%`
            );
        }
        styles.push(`${colorEntryToStyle(colors[n - 1])} ${100 * ((n - 1) / n)}%`);
    }

    return `linear-gradient(to right, ${styles.join(', ')})`;
}

export function getColorGradient(colors: ColorListEntry[]) {
    if (colors.length === 0) return 'linear-gradient(to right, #000 0%, #000 100%)';

    const hasOffsets = colors.every(c => Array.isArray(c));
    let styles;

    if (hasOffsets) {
        const off = [...colors] as [Color, number][];
        off.sort((a, b) => a[1] - b[1]);
        styles = off.map(c => colorEntryToStyle(c, true));
    } else {
        styles = colors.map(c => colorEntryToStyle(c));
    }

    return `linear-gradient(to right, ${styles.join(', ')})`;
}

function colorEntryToStyle(e: ColorListEntry, includeOffset = false) {
    if (Array.isArray(e)) {
        if (includeOffset) return `${Color.toStyle(e[0])} ${(100 * e[1]).toFixed(2)}%`;
        return Color.toStyle(e[0]);
    }
    return Color.toStyle(e);
}

export function parseColorList(input: string, separator: RegExp = /,/): ColorListEntry[] {
    const ret: ColorListEntry[] = [];
    const trimmed = input.replace(/\s+/g, '');
    let tokenStart = 0;
    let bracketLevel = 0;
    for (let i = 0, il = trimmed.length; i < il; ++i) {
        const c = trimmed[i];
        if (c === '(') {
            bracketLevel++;
            continue;
        } else if (c === ')') {
            if (bracketLevel > 0) {
                bracketLevel--;
            }
            continue;
        }

        if (bracketLevel > 0) continue;

        if (!separator.test(c)) {
            continue;
        }

        const color = trimmed.substring(tokenStart, i);
        tokenStart = i + 1;
        const decoded = decodeColor(color);
        if (decoded !== undefined) {
            ret.push(decoded);
        }
    }

    if (tokenStart < trimmed.length) {
        const color = trimmed.substring(tokenStart);
        const decoded = decodeColor(color);
        console.log(color, decoded);
        if (decoded !== undefined) {
            ret.push(decoded);
        }
    }

    return ret;
}
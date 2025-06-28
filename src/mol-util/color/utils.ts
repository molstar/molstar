/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color, ColorListEntry } from './color';

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
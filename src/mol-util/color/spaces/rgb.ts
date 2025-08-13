/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Color } from '../color';

export { Rgb };

interface Rgb extends Array<number> { [d: number]: number, '@type': 'normalized-rgb', length: 3 }

function Rgb() {
    return Rgb.zero();
}

namespace Rgb {
    export function zero(): Rgb {
        const out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out as Rgb;
    }

    export function fromColor(out: Rgb, hexColor: Color) {
        out[0] = (hexColor >> 16 & 255) / 255;
        out[1] = (hexColor >> 8 & 255) / 255;
        out[2] = (hexColor & 255) / 255;
        return out;
    }

    export function toColor(rgb: Rgb): Color {
        return (((rgb[0] * 255) << 16) | ((rgb[1] * 255) << 8) | (rgb[2] * 255)) as Color;
    }
}

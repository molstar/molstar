/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Color conversion and interpolation code adapted from chroma.js (https://github.com/gka/chroma.js)
 * Copyright (c) 2011-2018, Gregor Aisch, BSD license
 */

import type { Color } from '../color';
import { Rgb } from './rgb';

export { Hsl };

interface Hsl extends Array<number> { [d: number]: number, '@type': 'hsl', length: 3 }

function Hsl() {
    return Hsl.zero();
}

namespace Hsl {
    export function zero(): Hsl {
        const out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out as Hsl;
    }

    const _rgb = Rgb();
    export function fromColor(out: Hsl, color: Color) {
        Rgb.fromColor(_rgb, color);
        return Hsl.fromRgb(out, _rgb);
    }

    export function toColor(hsl: Hsl): Color {
        toRgb(_rgb, hsl);
        return Rgb.toColor(_rgb);
    }

    export function fromRgb(out: Hsl, rgb: Rgb) {
        const [r, g, b] = rgb;

        const minRgb = Math.min(r, g, b);
        const maxRgb = Math.max(r, g, b);

        const l = (maxRgb + minRgb) / 2;
        let s: number = 0, h: number = 0;

        if (maxRgb === minRgb) {
            s = 0;
            h = Number.NaN;
        } else {
            s = l < 0.5
                ? (maxRgb - minRgb) / (maxRgb + minRgb)
                : (maxRgb - minRgb) / (2 - maxRgb - minRgb);
        }

        if (r === maxRgb) h = (g - b) / (maxRgb - minRgb);
        else if (g === maxRgb) h = 2 + (b - r) / (maxRgb - minRgb);
        else if (b === maxRgb) h = 4 + (r - g) / (maxRgb - minRgb);

        h *= 60;
        if (h < 0) h += 360;
        out[0] = h;
        out[1] = s;
        out[2] = l;
        return out;
    }

    const _t3 = [0, 0, 0];
    const _c = [0, 0, 0];
    export function toRgb(out: Rgb, hsl: Hsl) {
        const [h, s, l] = hsl;
        let r: number, g: number, b: number;
        if (s === 0) {
            r = g = b = l;
        } else {
            const t3 = _t3;
            const c = _c;
            const t2 = l < 0.5 ? l * (1 + s) : l + s - l * s;
            const t1 = 2 * l - t2;
            const h_ = h / 360;
            t3[0] = h_ + 1 / 3;
            t3[1] = h_;
            t3[2] = h_ - 1 / 3;
            for (let i = 0; i < 3; i++) {
                if (t3[i] < 0) t3[i] += 1;
                if (t3[i] > 1) t3[i] -= 1;
                if (6 * t3[i] < 1) c[i] = t1 + (t2 - t1) * 6 * t3[i];
                else if (2 * t3[i] < 1) c[i] = t2;
                else if (3 * t3[i] < 2) c[i] = t1 + (t2 - t1) * (2 / 3 - t3[i]) * 6;
                else c[i] = t1;
            }
            r = c[0];
            g = c[1];
            b = c[2];
        }
        out[0] = r;
        out[1] = g;
        out[2] = b;
        return out;
    }

    export function interpolate(out: Hsl, col1: Hsl, col2: Hsl, t: number) {
        const xyz0 = col1, xyz1 = col2;

        const [hue0, sat0, lbv0] = xyz0;
        const [hue1, sat1, lbv1] = xyz1;

        let sat, hue, dh;

        if (!isNaN(hue0) && !isNaN(hue1)) {
            // both colors have hue
            if (hue1 > hue0 && hue1 - hue0 > 180) {
                dh = hue1 - (hue0 + 360);
            } else if (hue1 < hue0 && hue0 - hue1 > 180) {
                dh = hue1 + 360 - hue0;
            } else {
                dh = hue1 - hue0;
            }
            hue = hue0 + t * dh;
        } else if (!isNaN(hue0)) {
            hue = hue0;
            if ((lbv1 === 1 || lbv1 === 0)) sat = sat0;
        } else if (!isNaN(hue1)) {
            hue = hue1;
            if ((lbv0 === 1 || lbv0 === 0)) sat = sat1;
        } else {
            hue = Number.NaN;
        }

        if (sat === undefined) sat = sat0 + t * (sat1 - sat0);
        const lbv = lbv0 + t * (lbv1 - lbv0);
        out[0] = hue;
        out[1] = sat;
        out[2] = lbv;
        return out;
    }
}

/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Color conversion code adapted from chroma.js (https://github.com/gka/chroma.js)
 * Copyright (c) 2011-2018, Gregor Aisch, BSD license
 */

import { Color } from '../color';
import { Hcl } from './hcl';
import { radToDeg } from '../../../mol-math/misc';
import { clamp } from '../../../mol-math/interpolate';

export { Lab };

interface Lab extends Array<number> { [d: number]: number, '@type': 'lab', length: 3 }

/**
 * CIE LAB color
 *
 * - L* [0..100] - lightness from black to white
 * - a [-100..100] - green (-) to red (+)
 * - b [-100..100] - blue (-) to yellow (+)
 *
 * see https://en.wikipedia.org/wiki/CIELAB_color_space
 */
function Lab() {
    return Lab.zero();
}

namespace Lab {
    export function zero(): Lab {
        const out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out as Lab;
    }

    export function create(l: number, a: number, b: number): Lab {
        const out = zero();
        out[0] = l;
        out[1] = a;
        out[2] = b;
        return out;
    }

    export function fromColor(out: Lab, color: Color): Lab {
        const [r, g, b] = Color.toRgb(color);
        const [x, y, z] = rgbToXyz(r, g, b);
        const l = 116 * y - 16;
        out[0] = l < 0 ? 0 : l;
        out[1] = 500 * (x - y);
        out[2] = 200 * (y - z);
        return out;
    }

    export function fromHcl(out: Lab, hcl: Hcl): Lab {
        return Hcl.toLab(out, hcl);
    }

    export function toColor(lab: Lab): Color {
        let y = (lab[0] + 16) / 116;
        let x = isNaN(lab[1]) ? y : y + lab[1] / 500;
        let z = isNaN(lab[2]) ? y : y - lab[2] / 200;

        y = Yn * lab_xyz(y);
        x = Xn * lab_xyz(x);
        z = Zn * lab_xyz(z);

        const r = xyz_rgb(3.2404542 * x - 1.5371385 * y - 0.4985314 * z);  // D65 -> sRGB
        const g = xyz_rgb(-0.9692660 * x + 1.8760108 * y + 0.0415560 * z);
        const b = xyz_rgb(0.0556434 * x - 0.2040259 * y + 1.0572252 * z);

        return Color.fromRgb(
            Math.round(clamp(r, 0, 255)),
            Math.round(clamp(g, 0, 255)),
            Math.round(clamp(b, 0, 255))
        );
    }

    export function toHcl(out: Hcl, lab: Lab): Hcl {
        const [l, a, b] = lab;
        const c = Math.sqrt(a * a + b * b);
        let h = (radToDeg(Math.atan2(b, a)) + 360) % 360;
        if (Math.round(c * 10000) === 0) h = Number.NaN;
        out[0] = h;
        out[1] = c;
        out[2] = l;
        return out;
    }

    export function copy(out: Lab, c: Lab): Lab {
        out[0] = c[0];
        out[1] = c[1];
        out[2] = c[2];
        return out;
    }

    export function darken(out: Lab, c: Lab, amount: number): Lab {
        out[0] = c[0] - Kn * amount;
        out[1] = c[1];
        out[2] = c[2];
        return out;
    }

    export function lighten(out: Lab, c: Lab, amount: number): Lab {
        return darken(out, c, -amount);
    }

    const tmpSaturateHcl = [0, 0, 0] as Hcl;
    export function saturate(out: Lab, c: Lab, amount: number): Lab {
        toHcl(tmpSaturateHcl, c);
        return Hcl.toLab(out, Hcl.saturate(tmpSaturateHcl, tmpSaturateHcl, amount));
    }

    export function desaturate(out: Lab, c: Lab, amount: number): Lab {
        return saturate(out, c, -amount);
    }

    // Corresponds roughly to RGB brighter/darker
    const Kn = 18;

    /** D65 standard referent */
    const Xn = 0.950470;
    const Yn = 1;
    const Zn = 1.088830;

    const T0 = 0.137931034; // 4 / 29
    const T1 = 0.206896552; // 6 / 29
    const T2 = 0.12841855;  // 3 * t1 * t1
    const T3 = 0.008856452; // t1 * t1 * t1

    /** convert component from xyz to rgb */
    function xyz_rgb(c: number) {
        return 255 * (c <= 0.00304 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055);
    }

    /** convert component from lab to xyz */
    function lab_xyz(t: number) {
        return t > T1 ? t * t * t : T2 * (t - T0);
    }

    /** convert component from rgb to xyz */
    function rgb_xyz(c: number) {
        if ((c /= 255) <= 0.04045) return c / 12.92;
        return Math.pow((c + 0.055) / 1.055, 2.4);
    }

    /** convert component from xyz to lab */
    function xyz_lab(t: number) {
        if (t > T3) return Math.pow(t, 1 / 3);
        return t / T2 + T0;
    }

    function rgbToXyz(r: number, g: number, b: number) {
        r = rgb_xyz(r);
        g = rgb_xyz(g);
        b = rgb_xyz(b);
        const x = xyz_lab((0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / Xn);
        const y = xyz_lab((0.2126729 * r + 0.7151522 * g + 0.0721750 * b) / Yn);
        const z = xyz_lab((0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / Zn);
        return [x, y, z];
    }
}
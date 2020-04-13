/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Color conversion code adapted from chroma.js (https://github.com/gka/chroma.js)
 * Copyright (c) 2011-2018, Gregor Aisch, BSD license
 */

import { Color } from '../color';
import { degToRad } from '../../../mol-math/misc';
import { Lab } from './lab';

export { Hcl };

interface Hcl extends Array<number> { [d: number]: number, '@type': 'hcl', length: 3 }

/**
 * CIE HCL (Hue-Chroma-Luminance) color
 *
 * - H [0..360]
 * - C [0..100]
 * - L [0..100]
 *
 * Cylindrical representation of CIELUV (see https://en.wikipedia.org/wiki/CIELUV)
 */
function Hcl() {
    return Hcl.zero();
}

namespace Hcl {
    export function zero(): Hcl {
        const out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out as Hcl;
    }

    export function create(h: number, c: number, l: number): Hcl {
        const out = zero();
        out[0] = h;
        out[1] = c;
        out[2] = l;
        return out;
    }

    const tmpFromColorLab = [0, 0, 0] as Lab;
    export function fromColor(out: Hcl, color: Color): Hcl {
        return Lab.toHcl(out, Lab.fromColor(tmpFromColorLab, color));
    }

    export function fromLab(hcl: Hcl, lab: Lab): Hcl {
        return Lab.toHcl(hcl, lab);
    }

    const tmpToColorLab = [0, 0, 0] as Lab;
    export function toColor(hcl: Hcl): Color {
        return Lab.toColor(toLab(tmpToColorLab, hcl));
    }

    /**
     * Convert from a qualitative parameter h and a quantitative parameter l to a 24-bit pixel.
     *
     * These formulas were invented by David Dalrymple to obtain maximum contrast without going
     * out of gamut if the parameters are in the range 0-1.
     * A saturation multiplier was added by Gregor Aisch
     */
    export function toLab(out: Lab, hcl: Hcl): Lab {
        let [h, c, l] = hcl;
        if (isNaN(h)) h = 0;
        h = degToRad(h);
        out[0] = l;
        out[1] = Math.cos(h) * c;
        out[2] = Math.sin(h) * c;
        return out;
    }

    export function copy(out: Hcl, c: Hcl): Hcl {
        out[0] = c[0];
        out[1] = c[1];
        out[2] = c[2];
        return out;
    }

    export function saturate(out: Hcl, c: Hcl, amount: number): Hcl {
        out[0] = c[0];
        out[1] = Math.max(0, c[1] + Kn * amount);
        out[2] = c[2];
        return out;
    }

    export function desaturate(out: Hcl, c: Hcl, amount: number): Hcl {
        return saturate(out, c, -amount);
    }

    const tmpDarkenLab = [0, 0, 0] as Lab;
    export function darken(out: Hcl, c: Hcl, amount: number): Hcl {
        toLab(tmpDarkenLab, c);
        return Lab.toHcl(out, Lab.darken(tmpDarkenLab, tmpDarkenLab, amount));
    }

    export function lighten(out: Hcl, c: Hcl, amount: number): Hcl {
        return darken(out, c, -amount);
    }

    // Corresponds roughly to RGB brighter/darker
    const Kn = 18;
}
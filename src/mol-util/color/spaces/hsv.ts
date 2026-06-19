/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 *
 * Adapted from kin-parser.ts file from the NGL project:
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * Adapted from hsl.ts in this same directory:
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Color } from '../color';
import { Rgb } from './rgb';

export { Hsv };

/** Hsv tuple: [h, s, v]
 *  - h in [0,360] degrees
 *  - s in [0,100] percent
 *  - v in [0,100] percent
 */
interface Hsv extends Array<number> { [d: number]: number, '@type': 'hsv', length: 3 }

function Hsv() {
    return Hsv.zero();
}

namespace Hsv {
    export function zero(): Hsv {
        const out = [0.0, 0.0, 0.0];
        out[0] = 0;
        return out as Hsv;
    }

    /** Copy values from an array-like 3-tuple into `out`. */
    export function fromArray(arr: ArrayLike<number>): Hsv {
      const out = Hsv.zero();
      out[0] = arr[0] ?? 0;
      out[1] = arr[1] ?? 0;
      out[2] = arr[2] ?? 0;
      return out;
    }

    const _rgb = Rgb();
    export function toColor(hsv: Hsv): Color {
        toRgb(_rgb, hsv);
        return Rgb.toColor(_rgb);
    }

    export function toRgb(out: Rgb, hsv: Hsv) {
        let [h, s, v] = hsv;
        h /= 360;
        s /= 100;
        v /= 100;
        let r = 0, g = 0, b = 0;
        const i = Math.floor(h * 6);
        const f = h * 6 - i;
        const p = v * (1 - s);
        const q = v * (1 - f * s);
        const t = v * (1 - (1 - f) * s);
        switch (i % 6) {
          case 0: r = v; g = t; b = p; break;
          case 1: r = q; g = v; b = p; break;
          case 2: r = p; g = v; b = t; break;
          case 3: r = p; g = q; b = v; break;
          case 4: r = t; g = p; b = v; break;
          case 5: r = v; g = p; b = q; break;
        }
        out[0] = r;
        out[1] = g;
        out[2] = b;
        return out;
    }
}

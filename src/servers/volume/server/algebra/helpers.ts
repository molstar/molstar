import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';

// to prevent floating point rounding errors
export function round(v: number) {
    return Math.round(10000000 * v) / 10000000;
}

export function transform(x: { [index: number]: number }, matrix: Mat4) {
    return Vec3.transformMat4(Vec3.zero(), x as Vec3, matrix);
}

export function snap(v: number, to: 'bottom' | 'top') {
    switch (to) {
        case 'bottom': return Math.floor(round(v)) | 0;
        case 'top': return Math.ceil(round(v)) | 0;
    }
}
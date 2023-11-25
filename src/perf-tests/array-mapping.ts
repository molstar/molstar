import * as B from 'benchmark';
import { SymmetryOperator } from '../mol-math/geometry/symmetry-operator';
import { Mat4 } from '../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../mol-math/linear-algebra/3d/vec3';

interface CoordinateMapper<T extends number> { (index: T, slot: Vec3): Vec3 }
interface ArrayMapping<T extends number> {
    readonly coordinates: Coordinates,
    readonly operator: SymmetryOperator,
    readonly invariantPosition: CoordinateMapper<T>,
    readonly position: CoordinateMapper<T>,
    x(index: T): number,
    y(index: T): number,
    z(index: T): number,
    r(index: T): number
}

interface Coordinates { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }

function _createMapping<T extends number>(operator: SymmetryOperator, coords: Coordinates, radius: ((index: T) => number)): ArrayMapping<T> {
    const invariantPosition = createCoordinateMapper(SymmetryOperator.Default, coords);
    const position = operator.isIdentity ? invariantPosition : createCoordinateMapper(operator, coords);
    const { x, y, z } = createProjections(operator, coords);
    return { operator, coordinates: coords, invariantPosition, position, x, y, z, r: radius };
}

function createMapping<T extends number>(operator: SymmetryOperator, coords: Coordinates, radius: ((index: T) => number) = _zeroRadius) {
    return _createMapping(operator, coords, radius);
}

function createCoordinateMapper<T extends number>(t: SymmetryOperator, coords: Coordinates): CoordinateMapper<T> {
    if (t.isIdentity) return identityPosition(coords);
    return generalPosition(t, coords);
}

function _zeroRadius(i: number) { return 0; }

interface Projections { x(index: number): number, y(index: number): number, z(index: number): number }

function createProjections(t: SymmetryOperator, coords: SymmetryOperator.Coordinates): Projections {
    if (t.isIdentity) return { x: projectCoord(coords.x), y: projectCoord(coords.y), z: projectCoord(coords.z) };
    return { x: projectX(t, coords), y: projectY(t, coords), z: projectZ(t, coords) };
}

function projectCoord(xs: ArrayLike<number>) {
    return function projectCoord(i: number) {
        return xs[i];
    };
}

function isW1(m: Mat4) {
    return m[3] === 0 && m[7] === 0 && m[11] === 0 && m[15] === 1;
}

function projectX({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[0], yy = m[4], zz = m[8], tx = m[12];

    if (isW1(m)) {
        // this should always be the case.
        return function projectX_W1(i: number) {
            return xx * xs[i] + yy * ys[i] + zz * zs[i] + tx;
        };
    }

    return function projectX(i: number) {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tx) / w;
    };
}

function projectY({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[1], yy = m[5], zz = m[9], ty = m[13];

    if (isW1(m)) {
        // this should always be the case.
        return function projectY_W1(i: number) {
            return xx * xs[i] + yy * ys[i] + zz * zs[i] + ty;
        };
    }

    return function projectY(i: number) {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ty) / w;
    };
}

function projectZ({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[2], yy = m[6], zz = m[10], tz = m[14];

    if (isW1(m)) {
        // this should always be the case.
        return function projectZ_W1(i: number) {
            return xx * xs[i] + yy * ys[i] + zz * zs[i] + tz;
        };
    }

    return function projectZ(i: number) {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tz) / w;
    };
}

function identityPosition<T extends number>({ x, y, z }: SymmetryOperator.Coordinates): SymmetryOperator.CoordinateMapper<T> {
    return function identityPosition(i: T, s: Vec3): Vec3 {
        s[0] = x[i];
        s[1] = y[i];
        s[2] = z[i];
        return s;
    };
}

function generalPosition<T extends number>({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    if (isW1(m)) {
        // this should always be the case.
        return function generalPosition_W1(i: T, r: Vec3): Vec3 {
            const x = xs[i], y = ys[i], z = zs[i];
            r[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
            r[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
            r[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
            return r;
        };
    }
    return function generalPosition(i: T, r: Vec3): Vec3 {
        r[0] = xs[i];
        r[1] = ys[i];
        r[2] = zs[i];
        Vec3.transformMat4(r, r, m);
        return r;
    };
}

//

function getData(size: number) {
    const op = SymmetryOperator.create('id', Mat4.identity());
    const coords = {
        x: new Float32Array(size),
        y: new Float32Array(size),
        z: new Float32Array(size),
    };
    return { op, coords };
}

function getArrayMappingClass(op: SymmetryOperator, coords: Coordinates) {
    return SymmetryOperator.createMapping(op, coords);
}

function getArrayMappingClosures(op: SymmetryOperator, coords: Coordinates) {
    return createMapping(op, coords);
}

function testXYZ(mapping: SymmetryOperator.ArrayMapping<number>) {
    const n = mapping.coordinates.x.length;
    let sum = 0;

    for (let i = 0; i < n; i++) {
        sum += mapping.x(i);
        sum += mapping.y(i);
        sum += mapping.z(i);
    }

    return sum /= n;
}

function testPos(mapping: SymmetryOperator.ArrayMapping<number>) {
    const n = mapping.coordinates.x.length;
    let sum = 0;
    const p = Vec3();

    for (let i = 0; i < n; i++) {
        mapping.position(i, p);
        sum += p[0];
        sum += p[1];
        sum += p[2];
    }

    return sum /= n;
}

function testInvPos(mapping: SymmetryOperator.ArrayMapping<number>) {
    const n = mapping.coordinates.x.length;
    let sum = 0;
    const p = Vec3();

    for (let i = 0; i < n; i++) {
        mapping.invariantPosition(i, p);
        sum += p[0];
        sum += p[1];
        sum += p[2];
    }

    return sum /= n;
}

function run(size: number) {
    const suite = new B.Suite();

    const { op, coords } = getData(size);
    const mappingClass = getArrayMappingClass(op, coords);
    const mappingClosures = getArrayMappingClosures(op, coords);

    suite
        .add('class testXYZ', () => testXYZ(mappingClass))
        .add('class testPos', () => testPos(mappingClass))
        .add('class testInvPos', () => testInvPos(mappingClass))
        .add('closures testXYZ', () => testXYZ(mappingClosures))
        .add('closures testPos', () => testPos(mappingClosures))
        .add('closures testInvPos', () => testInvPos(mappingClosures))
        .on('cycle', (e: any) => {
            console.log(size, String(e.target));
        })
        .run();
}

run(10000);
run(100000);
run(1000000);
run(10000000);

// 10000 class testXYZ x 48,383 ops/sec ±0.39% (94 runs sampled)
// 10000 class testPos x 56,421 ops/sec ±0.50% (98 runs sampled)
// 10000 class testInvPos x 59,019 ops/sec ±0.15% (100 runs sampled)
// 10000 closures testXYZ x 7,788 ops/sec ±6.39% (91 runs sampled)
// 10000 closures testPos x 35,736 ops/sec ±11.59% (98 runs sampled)
// 10000 closures testInvPos x 35,501 ops/sec ±11.78% (95 runs sampled)
// 100000 class testXYZ x 603 ops/sec ±0.40% (95 runs sampled)
// 100000 class testPos x 1,391 ops/sec ±0.45% (98 runs sampled)
// 100000 class testInvPos x 1,650 ops/sec ±0.73% (92 runs sampled)
// 100000 closures testXYZ x 724 ops/sec ±0.39% (95 runs sampled)
// 100000 closures testPos x 1,612 ops/sec ±0.47% (94 runs sampled)
// 100000 closures testInvPos x 1,626 ops/sec ±0.32% (95 runs sampled)
// 1000000 class testXYZ x 55.15 ops/sec ±0.49% (59 runs sampled)
// 1000000 class testPos x 132 ops/sec ±0.52% (77 runs sampled)
// 1000000 class testInvPos x 181 ops/sec ±0.85% (85 runs sampled)
// 1000000 closures testXYZ x 70.03 ops/sec ±0.54% (73 runs sampled)
// 1000000 closures testPos x 163 ops/sec ±0.53% (84 runs sampled)
// 1000000 closures testInvPos x 156 ops/sec ±0.42% (81 runs sampled)
// 10000000 class testXYZ x 6.09 ops/sec ±1.32% (20 runs sampled)
// 10000000 class testPos x 14.15 ops/sec ±0.49% (40 runs sampled)
// 10000000 class testInvPos x 18.02 ops/sec ±0.59% (49 runs sampled)
// 10000000 closures testXYZ x 7.30 ops/sec ±0.69% (23 runs sampled)
// 10000000 closures testPos x 16.39 ops/sec ±0.57% (45 runs sampled)
// 10000000 closures testInvPos x 15.23 ops/sec ±0.53% (42 runs sampled)

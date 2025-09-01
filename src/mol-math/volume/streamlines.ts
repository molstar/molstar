// mol-math/volume/streamlines.ts
import { Tensor, Vec3, Mat4 } from '../../mol-math/linear-algebra';

export type StreamlinePoint = { x: number; y: number; z: number; value: number };
export type Streamline = StreamlinePoint[];
export type StreamlineSet = Streamline[];

export type StreamlineMode = 'simple' | 'advanced';

export interface StreamlineParams {
  seedDensity: number;
  maxSteps: number;
  stepSize: number; // Å
  minSpeed: number; // 1/Å  (|E|)
  minLevel: number;
  maxLevel: number;
  mid_interpolation: boolean; // RK2 vs RK4 for simple
  writeStride: number; // tracer subsampling
}

export interface GridInfo {
  space: Tensor.Space; // volume.grid.cells.space
  data: ArrayLike<number>; // volume.grid.cells.data
  dims: Vec3; // [nx, ny, nz]
  cellSize: Vec3; // [hx, hy, hz] in Å
  gridToCartn: Mat4; // transformation matrix from grid to Cartesian coordinates
}

/** Pure function: returns polylines in INDEX space (no transforms, no builders). */
export function collectStreamlines(
  grid: GridInfo,
  params: StreamlineParams,
  mode: StreamlineMode,
  out: StreamlineSet = []
): StreamlineSet {
  return (mode === 'advanced')
    ? computeAdvancedStreamlines(grid, params, out)
    : computeSimpleStreamlines(grid, params, out);
}

function getCellSize(gridToCartn: Mat4) {
    const b = Mat4.extractBasis(gridToCartn);
    return { hx: Vec3.magnitude(b.x), hy: Vec3.magnitude(b.y), hz: Vec3.magnitude(b.z) };
}

// ---------------- Simple tracer (RK2/RK4) in index space ----------------

function traceOneDirection(
    out: StreamlinePoint[],
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, dsWorld: number, eps: number,
    hx: number, hy: number, hz: number, dirSign: 1 | -1, skipFirst: boolean,
    S: { p: Vec3, g: Vec3, g2: Vec3, g3: Vec3, g4: Vec3, v1: Vec3, v2: Vec3, v3: Vec3, v4: Vec3, k1: Vec3, k2: Vec3, k3: Vec3, k4: Vec3, dp: Vec3, t1: Vec3, t2: Vec3, t3: Vec3 },
    writeStride = 2,
    order: 2 | 4
): number {
    const [nx, ny, nz] = space.dimensions as Vec3;
    const p = S.p; Vec3.copy(p, seed);

    // convert a world step to index steps per-axis
    const sIx = dirSign * (dsWorld / hx), sIy = dirSign * (dsWorld / hy), sIz = dirSign * (dsWorld / hz);

    let written = 0;
    for (let step = 0; step < maxSteps; step++) {
        if (p[0] < 1 || p[0] > (nx as number) - 2 || p[1] < 1 || p[1] > (ny as number) - 2 || p[2] < 1 || p[2] > (nz as number) - 2) break;

        // gW = -∇φ in WORLD units (1/Å)
        gradientAtP_world(S.g, space, data, p, hx, hy, hz, S.t3);
        const m = Vec3.magnitude(S.g);
        if (!(m > eps)) break;
        Vec3.scale(S.v1, S.g, 1 / m); // unit world

        if (order === 2) {
            // RK2 midpoint in world, converted to index steps
            S.t1[0] = p[0] + 0.5 * sIx * S.v1[0];
            S.t1[1] = p[1] + 0.5 * sIy * S.v1[1];
            S.t1[2] = p[2] + 0.5 * sIz * S.v1[2];

            gradientAtP_world(S.g2, space, data, S.t1, hx, hy, hz, S.t3);
            const m2 = Math.max(Vec3.magnitude(S.g2), eps);
            S.v2[0] = S.g2[0] / m2; S.v2[1] = S.g2[1] / m2; S.v2[2] = S.g2[2] / m2;

            S.dp[0] = sIx * S.v2[0]; S.dp[1] = sIy * S.v2[1]; S.dp[2] = sIz * S.v2[2];
        } else {
            // RK4 in world, converted to index
            Vec3.set(S.k1, sIx * S.v1[0], sIy * S.v1[1], sIz * S.v1[2]);

            S.t1[0] = p[0] + 0.5 * S.k1[0]; S.t1[1] = p[1] + 0.5 * S.k1[1]; S.t1[2] = p[2] + 0.5 * S.k1[2];
            gradientAtP_world(S.g2, space, data, S.t1, hx, hy, hz, S.t3);
            const m2 = Math.max(Vec3.magnitude(S.g2), eps);
            S.v2[0] = S.g2[0] / m2; S.v2[1] = S.g2[1] / m2; S.v2[2] = S.g2[2] / m2;
            Vec3.set(S.k2, sIx * S.v2[0], sIy * S.v2[1], sIz * S.v2[2]);

            S.t1[0] = p[0] + 0.5 * S.k2[0]; S.t1[1] = p[1] + 0.5 * S.k2[1]; S.t1[2] = p[2] + 0.5 * S.k2[2];
            gradientAtP_world(S.g3, space, data, S.t1, hx, hy, hz, S.t3);
            const m3 = Math.max(Vec3.magnitude(S.g3), eps);
            S.v3[0] = S.g3[0] / m3; S.v3[1] = S.g3[1] / m3; S.v3[2] = S.g3[2] / m3;
            Vec3.set(S.k3, sIx * S.v3[0], sIy * S.v3[1], sIz * S.v3[2]);

            S.t1[0] = p[0] + S.k3[0]; S.t1[1] = p[1] + S.k3[1]; S.t1[2] = p[2] + S.k3[2];
            gradientAtP_world(S.g4, space, data, S.t1, hx, hy, hz, S.t3);
            const m4 = Math.max(Vec3.magnitude(S.g4), eps);
            S.v4[0] = S.g4[0] / m4; S.v4[1] = S.g4[1] / m4; S.v4[2] = S.g4[2] / m4;
            Vec3.set(S.k4, sIx * S.v4[0], sIy * S.v4[1], sIz * S.v4[2]);

            // dp = (k1 + 2*k2 + 2*k3 + k4) / 6
            Vec3.copy(S.dp, S.k1);
            Vec3.scaleAndAdd(S.dp, S.dp, S.k2, 2);
            Vec3.scaleAndAdd(S.dp, S.dp, S.k3, 2);
            Vec3.add(S.dp, S.dp, S.k4);
            Vec3.scale(S.dp, S.dp, 1 / 6);
        }

        if (!(skipFirst && step === 0) && (step % writeStride === 0)) {
            const value = getInterpolatedValue(space, data, p);
            out.push({ x: p[0], y: p[1], z: p[2], value });
            written++;
        }
        Vec3.add(p, p, S.dp);
    }
    return written;
}

function reverseSegmentInPlace<T>(arr: T[], start: number, end: number) {
    for (let i = start, j = end; i < j; i++, j--) { const tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp; }
}

// safe to hoist once at module level
const Scratch = {
  p: Vec3(), g: Vec3(), g2: Vec3(), g3: Vec3(), g4: Vec3(),
  v1: Vec3(), v2: Vec3(), v3: Vec3(), v4: Vec3(),
  k1: Vec3(), k2: Vec3(), k3: Vec3(), k4: Vec3(),
  dp: Vec3(), t1: Vec3(), t2: Vec3(), t3: Vec3()
};

function traceStreamlineBothDirs(space: Tensor.Space, data: ArrayLike<number>,
    seed: Vec3, maxSteps: number, dsWorld: number, eps: number, hx: number,
    hy: number, hz: number, writeStride: number, order: 2 | 4): StreamlinePoint[] {
    const line: StreamlinePoint[] = [];
    const S = Scratch;

    const nBack = traceOneDirection(line, space, data, seed, maxSteps, dsWorld, eps, hx, hy, hz, -1, false, S, writeStride, order);
    if (nBack > 1) reverseSegmentInPlace(line, 0, nBack - 1);
    traceOneDirection(line, space, data, seed, maxSteps, dsWorld, eps, hx, hy, hz, +1, true, S, writeStride, order);
    return line;
}

// --- put your optimized simple/advanced tracer implementation here ---
// Important: no imports from mol-model/* or mol-repr/*; only math & types.
function computeSimpleStreamlines(grid: GridInfo, props: StreamlineParams, out: StreamlineSet): StreamlineSet {
    const { space, data, gridToCartn } = grid;

    const [nx, ny, nz] = space.dimensions as Vec3;
    const { hx, hy, hz } = getCellSize(gridToCartn);

    const seedDensity = props.seedDensity ?? 8;
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    const maxSteps = props.maxSteps ?? 2000;
    const dsWorld = props.stepSize ?? 0.35; // Å
    const eps = props.minSpeed ?? 1e-6; // |E| threshold in 1/Å
    const writeStride = Math.max(1, props.writeStride | 0);
    const order: 2 | 4 = props.mid_interpolation ? 2 : 4;

    const pos = Vec3(), g = Vec3();

    // bounds in index space avoiding edges
    const xStart = 1, xEnd = nx - 2;
    const yStart = 1, yEnd = ny - 2;
    const zStart = 1, zEnd = nz - 2;

    // seeds on a coarse lattice, shuffled
    const seeds: number[] = [];
    for (let z = zStart; z <= zEnd; z += seedStep)
        for (let y = yStart; y <= yEnd; y += seedStep)
            for (let x = xStart; x <= xEnd; x += seedStep) { seeds.push(x, y, z); }
    for (let i = seeds.length - 3; i > 0; i -= 3) {
        const j = (Math.floor(Math.random() * (i / 3 + 1)) * 3) | 0;
        const tx = seeds[i], ty = seeds[i+1], tz = seeds[i+2];
        seeds[i] = seeds[j]; seeds[i+1] = seeds[j+1]; seeds[i+2] = seeds[j+2];
        seeds[j] = tx; seeds[j+1] = ty; seeds[j+2] = tz;
    }

    for (let i = 0; i < seeds.length; i += 3) {
        Vec3.set(pos, seeds[i], seeds[i + 1], seeds[i + 2]);

        const phi = getInterpolatedValue(space, data, pos);
        gradientAtP_world(g, space, data, pos, hx, hy, hz);
        const mag = Vec3.magnitude(g);
        if ((phi < props.minLevel || phi > props.maxLevel) || mag < eps) continue;

        const line = traceStreamlineBothDirs(space, data, pos, maxSteps, dsWorld, eps, hx, hy, hz, writeStride, order);
        if (line.length >= 2) out.push(line);
    }
    return out;
}

// ---------------- Advanced (PyMOL-style) tracer ----------------
// Algorithm based on PyMOL’s isosurface gradient tracing (Schrödinger LLC),
// reimplemented independently.

const OBLATION_SPACING = 2; // cells

function computeAdvancedStreamlines(grid: GridInfo, props: StreamlineParams, out: StreamlineSet): StreamlineSet {
    const { space, data, gridToCartn } = grid;
    const [nx, ny, nz] = space.dimensions as Vec3;
    const { hx, hy, hz } = getCellSize(gridToCartn);

    // map range (full grid)
    const i0 = 0, j0 = 0, k0 = 0;
    const i1 = nx, j1 = ny, k1 = nz;

    // precompute WORLD gradients (units: 1/Å) on the whole grid
    const grad = buildGradientFieldWorld(space, data, hx, hy, hz);

    // index → world (we store index coords and transform at geometry stage)
    const toIndexPoint = (i:number, j:number, k:number, fx:number, fy:number, fz:number, outV:Vec3) => { Vec3.set(outV, i + fx, j + fy, k + fz); };

    // randomized order of cells
    const rx = i1 - i0, ry = j1 - j0, rz = k1 - k0;
    const rangeSize = rx * ry * rz;
    const order = new Int32Array(rangeSize * 3);
    {
        let t = 0;
        for (let k = k0; k < k1; k++) for (let j = j0; j < j1; j++) for (let i = i0; i < i1; i++) { order[t++] = i; order[t++] = j; order[t++] = k; }
        // deterministic but good shuffle
        let s = (nx * ny * nz) >>> 0; const rnd = () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; };
        for (let a = 0; a < rangeSize; a++) {
            const p = (Math.floor(rnd() * rangeSize) * 3) | 0;
            const q = (Math.floor(rnd() * rangeSize) * 3) | 0;
            const t0 = order[p], t1 = order[p+1], t2 = order[p+2]; order[p] = order[q]; order[p+1] = order[q+1]; order[p+2] = order[q+2]; order[q] = t0; order[q+1] = t1; order[q+2] = t2;
        }
    }

    // oblation flags
    const flag = new Uint8Array(rx * ry * rz);
    const strideX = 1, strideY = rx, strideZ = rx * ry;
    const flagIndex = (i:number, j:number, k:number) => ((i - i0) * strideX) + ((j - j0) * strideY) + ((k - k0) * strideZ);

    const minDot = 0.0; // disallow sharp reversals
    const maxSteps = Math.max(1, props.maxSteps | 0);
    const dsWorld = props.stepSize; // Å
    const writeStride = Math.max(1, props.writeStride | 0);
    const minSlope = Math.max(props.minSpeed, 1e-6); // 1/Å

    // helpers to keep (fx,fy,fz) in [0,1) while moving voxel indices
    function normalizeFrac(f:number, l:number, lo:number, hi:number) {
        while (f < 0) { f += 1; l--; }
        while (f >= 1) { f -= 1; l++; }
        return (l < lo || l > hi) ? { f, l, done: true } : { f, l, done: false };
    }

    // emit segments similar to PyMOL bookkeeping but directly into out as separate polylines
    const pushNewPolyline = () => { out.push([]); return out.length - 1; };
    // scratch (NO allocations in the inner loop)
    const prev = Vec3();
    const gW = Vec3();
    const pw = Vec3();
    const activeCells: number[] = [];

    for (let a = 0; a < rangeSize; a++) {
        activeCells.length = 0; // reuse
        // two passes: +gradient, -gradient (no closures)
        for (let pass = 0; pass < 2; pass++) {
            const sign: 1 | -1 = pass === 0 ? +1 : -1;
            let li = order[a*3], lj = order[a*3+1], lk = order[a*3+2];
            let fx = 0.5, fy = 0.5, fz = 0.5; // start at voxel center (stabler)
            let nVert = 0; let havePrev = false;
            const segIdx = pushNewPolyline();

            for (let step = 0; step < maxSteps; step++) {
                // keep within valid trilinear range [0..N-2]
                let r = normalizeFrac(fx, li, i0, i1 - 2); fx = r.f; li = r.l; if (r.done) break;
                r = normalizeFrac(fy, lj, j0, j1 - 2); fy = r.f; lj = r.l; if (r.done) break;
                r = normalizeFrac(fz, lk, k0, k1 - 2); fz = r.f; lk = r.l; if (r.done) break;

                if (flag[flagIndex(li, lj, lk)]) break;

                const level = interpolateValueFrac(space, data, li, lj, lk, fx, fy, fz);
                if (level < props.minLevel || level > props.maxLevel) break;

                // world gradient (unit dir with sign)
                interpolateGradientWorld(space, grad, li, lj, lk, fx, fy, fz, gW);
                const gm = Vec3.magnitude(gW);
                if (gm < minSlope) break;
                Vec3.scale(gW, gW, sign / gm); // unit world dir with sign

                // decimate writes in tracer: only record every Nth step
                if ((step % writeStride) === 0) {
                    // add point (index coords; geometry stage does gridToCartn)
                    toIndexPoint(li, lj, lk, fx, fy, fz, pw);
                    out[segIdx].push({ x: pw[0], y: pw[1], z: pw[2], value: level });
                    nVert++;
                }

                // record visited cell for oblation (store once per change)
                const n = activeCells.length;
                if (n < 3 || activeCells[n-3] !== li || activeCells[n-2] !== lj || activeCells[n-1] !== lk) {
                    activeCells.push(li, lj, lk);
                }

                // directional coherence (avoid flip-flop)
                if (havePrev) {
                    const dp = gW[0]*prev[0] + gW[1]*prev[1] + gW[2]*prev[2];
                    if (dp < minDot) break;
                }
                Vec3.copy(prev, gW); havePrev = true;

                // advance in INDEX coordinates using a world step
                fx += gW[0] * (dsWorld / hx);
                fy += gW[1] * (dsWorld / hy);
                fz += gW[2] * (dsWorld / hz);
            }

            // discard degenerate polylines
            if (nVert < 2) {
                out[segIdx] = [] as StreamlinePoint[]; // keep slot empty; harmless to geometry stage
            }
        }
        // oblation (flag cells in a sphere of radius R around visited cells)
        const R = OBLATION_SPACING|0; const R2 = R*R;
        for (let t = 0; t < activeCells.length; t += 3) {
            const ii = activeCells[t], jj = activeCells[t+1], kk = activeCells[t+2];
            for (let k = Math.max(k0, kk-R); k <= Math.min(k1-1, kk+R); k++) {
                const dz2 = (kk-k)*(kk-k);
                for (let j = Math.max(j0, jj-R); j <= Math.min(j1-1, jj+R); j++) {
                    const dy2 = (jj-j)*(jj-j) + dz2; if (dy2 > R2) continue;
                    for (let i = Math.max(i0, ii-R); i <= Math.min(i1-1, ii+R); i++) {
                        const d2 = (ii-i)*(ii-i) + dy2; if (d2 <= R2) flag[flagIndex(i,j,k)] = 1;
                    }
                }
            }
        }
    }
    return out;
}


// ---------------------------------------------------------------------------------------
// Field sampling (clamped + world-aware gradients)
// ---------------------------------------------------------------------------------------

// Trilinear interpolation from a Vec3 pos (index coords)
function getInterpolatedValue(space: Tensor.Space, data: ArrayLike<number>, pos: Vec3): number {
    const x = Math.floor(pos[0]), y = Math.floor(pos[1]), z = Math.floor(pos[2]);
    const fx = pos[0] - x, fy = pos[1] - y, fz = pos[2] - z;
    let acc = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        acc += data[space.dataOffset(x + dx, y + dy, z + dz)] * w;
    }
    return acc;
}

// Same but with explicit (i,j,k) + fractional offsets
function interpolateValueFrac(space: Tensor.Space, data: ArrayLike<number>, i: number, j: number, k: number, fx: number, fy: number, fz: number) {
    let acc = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        acc += data[space.dataOffset(i + dx, j + dy, k + dz)] * w;
    }
    return acc;
}

// Central-diff gradient of φ in WORLD units (1/Å), then E = -∇φ
function gradientAtP_world(
  out: Vec3,
  space: Tensor.Space,
  data: ArrayLike<number>,
  p: Vec3,
  hx: number, hy: number, hz: number
): Vec3 {
  // grid dims
  const nx = space.dimensions[0] as number;
  const ny = space.dimensions[1] as number;
  const nz = space.dimensions[2] as number;

  // base voxel (clamped)
  // |0 is a fast floor for non-negative numbers; we clamp right after anyway.
  let x0 = p[0] | 0; if (x0 < 0) x0 = 0; else if (x0 > nx - 1) x0 = nx - 1;
  let y0 = p[1] | 0; if (y0 < 0) y0 = 0; else if (y0 > ny - 1) y0 = ny - 1;
  let z0 = p[2] | 0; if (z0 < 0) z0 = 0; else if (z0 > nz - 1) z0 = nz - 1;

  // neighbor indices (clamped)
  const xp = x0 + 1 <= nx - 1 ? x0 + 1 : nx - 1;
  const xm = x0 - 1 >= 0 ? x0 - 1 : 0;
  const yp = y0 + 1 <= ny - 1 ? y0 + 1 : ny - 1;
  const ym = y0 - 1 >= 0 ? y0 - 1 : 0;
  const zp = z0 + 1 <= nz - 1 ? z0 + 1 : nz - 1;
  const zm = z0 - 1 >= 0 ? z0 - 1 : 0;

  // central differences of φ, then E = -∇φ (world units)
  const φxp = data[space.dataOffset(xp, y0, z0)];
  const φxm = data[space.dataOffset(xm, y0, z0)];
  const φyp = data[space.dataOffset(x0, yp, z0)];
  const φym = data[space.dataOffset(x0, ym, z0)];
  const φzp = data[space.dataOffset(x0, y0, zp)];
  const φzm = data[space.dataOffset(x0, y0, zm)];

  out[0] = - (φxp - φxm) / (2 * hx);
  out[1] = - (φyp - φym) / (2 * hy);
  out[2] = - (φzp - φzm) / (2 * hz);
  return out;
}

// Precompute WORLD gradient field (E = -∇φ) at voxel centers
function buildGradientFieldWorld(space: Tensor.Space, data: ArrayLike<number>, hx: number, hy: number, hz: number) {
    const nx = space.dimensions[0] as number, ny = space.dimensions[1] as number, nz = space.dimensions[2] as number;
    const gx = new Float32Array(nx * ny * nz);
    const gy = new Float32Array(nx * ny * nz);
    const gz = new Float32Array(nx * ny * nz);
    const idx = (x: number, y: number, z: number) => space.dataOffset(x, y, z);

    for (let x = 0; x < nx; x++) {
        const xm = Math.max(0, x - 1), xp = Math.min(nx - 1, x + 1);
        for (let y = 0; y < ny; y++) {
            const ym = Math.max(0, y - 1), yp = Math.min(ny - 1, y + 1);
            for (let z = 0; z < nz; z++) {
                const zm = Math.max(0, z - 1), zp = Math.min(nz - 1, z + 1);
                const off = idx(x,y,z);
                gx[off] = - (data[idx(xp, y, z)] - data[idx(xm, y, z)]) / (2 * hx);
                gy[off] = - (data[idx(x, yp, z)] - data[idx(x, ym, z)]) / (2 * hy);
                gz[off] = - (data[idx(x, y, zp)] - data[idx(x, y, zm)]) / (2 * hz);
            }
        }
    }
    return { gx, gy, gz };
}

// Trilinear interpolation of the WORLD gradient
function interpolateGradientWorld(space: Tensor.Space, grad: { gx: Float32Array, gy: Float32Array, gz: Float32Array }, i: number, j: number, k: number, fx: number, fy: number, fz: number, out: Vec3) {
    const { gx, gy, gz } = grad;
    let vx = 0, vy = 0, vz = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        const off = space.dataOffset(i + dx, j + dy, k + dz);
        vx += gx[off] * w; vy += gy[off] * w; vz += gz[off] * w;
    }
    out[0] = vx; out[1] = vy; out[2] = vz; return out;
}
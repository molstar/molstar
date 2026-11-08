/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { BufferUsage, MapMode, TextureUsage } from '@molstar/graphics/gl/webgpu/constants';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { getLodPixelBoundary, getStrides, type LodLevels } from '@molstar/graphics/geo/geometry/lod-levels';
import { Spheres as SpheresGeometry } from '@molstar/graphics/geo/geometry/spheres/spheres';
import {
  createTargetParticleTransform,
  ParticleTargetRepresentationParams,
} from '@molstar/graphics/repr/particles/representation/target/common';
import {
  getSpacefillParams,
  SpacefillRepresentationProvider,
} from '@molstar/graphics/repr/structure/representation/spacefill';
import { getPhysicalRadius } from '@molstar/graphics/theme/size/physical';
import { getParticleTargetGroups, type ParticleList } from '@molstar/model/model/particles/particle-list';
import type { ElementIndex, Structure, Unit } from '@molstar/model/model/structure';
import { HiZParams } from '@molstar/graphics/canvas3d/passes/hi-z';
import { structureRepr } from '../visual/scenes/structure.js';
import { type BenchmarkCase, BenchmarkEnv, type PathRun } from './framework.js';
import {
  getGlobule,
  getPostsyn,
  getBenchmarkEntry,
  getScene,
  LargeSize,
  postsynRepr,
  type RenderInput,
  renderFrames,
  rotateCamera,
  Size,
} from './render.js';

/**
 * Prototype of a compute raster for small opaque spheres (WebGPU with `atomic-vec2u-min-max`), measured against the
 * sphere impostors' opaque pass on the same scene and camera: a kernel per sphere writes each covered pixel's nearest
 * sphere into a visibility buffer with one 64-bit `atomicStoreMax` (depth in the high word `y`, id in `x`), a full-screen
 * resolve recomputes the hit and shades each pixel once. Variants: an early-out on a 32-bit depth buffer before the
 * 64-bit write, and two 32-bit passes (depth, then id where equal) as a fallback's cost. Spheres above `R` pixels in
 * radius are skipped (left to the impostors). The spheres are read from the data as the impostors draw them, instanced:
 * the atoms of each target once, a transform per instance (a particle's, or the identity), with the impostors' levels of
 * detail by distance (the atoms in level order, a level a prefix of them, its spheres scaled up): the level picked per
 * instance by the distance of its center (the impostors pick per sphere), whole instances culled by the frustum, a
 * workgroup per chunk of an instance's atoms. The compacted variant culls and picks levels in a pass per instance that
 * appends chunks of 64 spheres, then rasters the chunks in an indirect dispatch (no workgroups without spheres); on the
 * same chunks with exact pixel bounds, and with a workgroup tile; then with exact bounds and the instances also culled
 * by Hi-Z in two phases (the instances visible in the previous frame, a Hi-Z of their pixels, then every instance)
 */

interface Spheres {
  /** object-space centers and radii (`vec4` each), the targets' atoms one after the other */
  readonly atoms: Float32Array;
  /** a `mat4` per instance, the targets' instances one after the other */
  readonly transforms: Float32Array;
  readonly targets: readonly Target[];
  /** levels by distance, finest first (none: every atom at its size) */
  readonly levels: readonly { minDistance: number; maxDistance: number }[];
}

interface Target {
  atomOffset: number;
  atomCount: number;
  instanceOffset: number;
  instanceCount: number;
  /** the bounding sphere of the atoms (object space), its radius with the largest sphere at the coarsest scale */
  center: Vec3;
  radius: number;
  /** per level (finest first) the atoms drawn (a prefix) and their scale */
  counts: number[];
  scales: number[];
}

/**
 * `atoms` ordered into the levels of `lodLevels` (coarsest first) as the impostors order them (`orderPositionGroup`,
 * `getLodLevelsValue`): per level its count and scale, finest first; one level of every atom without levels
 */
function levelOrder(atoms: Float32Array, lodLevels: LodLevels, sizeFactor: number) {
  const n = atoms.length / 4;
  if (lodLevels.length === 0) return { atoms, counts: [n], scales: [1] };
  const strides = getStrides(lodLevels, sizeFactor);
  const levels = new Uint8Array(n);
  for (let j = 0; j < n; ++j) {
    let k = 0;
    while (k < strides.length && j % strides[k] !== 0) ++k;
    levels[j] = k;
  }
  const out = new Float32Array(atoms.length);
  const ends: number[] = [];
  let o = 0;
  for (let i = 0; i < strides.length; ++i) {
    for (let j = 0; j < n; ++j) {
      if (levels[j] === i) out.set(atoms.subarray(j * 4, j * 4 + 4), 4 * o++);
    }
    ends.push(o);
  }
  // finest first: level i draws the items up to the end of stride `strides[strides.length - 1 - i]`
  const counts = lodLevels.map((_, i) => ends[strides.length - 1 - i]);
  const scales = lodLevels.map((l, i) => Math.pow(Math.min(n, strides[strides.length - 1 - i]), 1 / l.scaleBias));
  return { atoms: out, counts, scales };
}

function sphereTarget(
  sphere: Sphere3D,
  atoms: Float32Array,
  lodLevels: LodLevels,
  sizeFactor: number,
  offsets: { atomOffset: number; instanceOffset: number; instanceCount: number },
): Target & { atoms: Float32Array } {
  const ordered = levelOrder(atoms, lodLevels, sizeFactor);
  let maxRadius = 0;
  for (let i = 3; i < atoms.length; i += 4) maxRadius = Math.max(maxRadius, atoms[i]);
  const { center, radius } = sphere;
  return {
    ...offsets,
    ...ordered,
    atomCount: atoms.length / 4,
    center: Vec3.clone(center),
    radius: radius + maxRadius * Math.max(...ordered.scales),
  };
}

interface RasterInput extends RenderInput {
  readonly spheres: () => Spheres;
  readonly label: string;
}

/** Atoms of `structure`, radii by `radius` (as the representation sizes them) */
function structureAtoms(structure: Structure, radius: (unit: Unit, element: ElementIndex) => number) {
  const out = new Float32Array(structure.elementCount * 4);
  const p = Vec3();
  let o = 0;
  for (const unit of structure.units) {
    const { elements, conformation } = unit;
    for (let i = 0, il = elements.length; i < il; ++i, o += 4) {
      conformation.position(elements[i], p);
      out[o] = p[0];
      out[o + 1] = p[1];
      out[o + 2] = p[2];
      out[o + 3] = radius(unit, elements[i]);
    }
  }
  return out;
}

function concat(arrays: Float32Array[]) {
  const out = new Float32Array(arrays.reduce((n, a) => n + a.length, 0));
  for (let i = 0, o = 0; i < arrays.length; o += arrays[i].length, ++i) out.set(arrays[i], o);
  return out;
}

/**
 * Spacefill of `structure`: physical radii times the representation's size factor (2 for coarse-grained structures);
 * without levels
 */
function structureSpheres(structure: Structure): Spheres {
  const sizeFactor = getSpacefillParams(undefined as any, structure).sizeFactor.defaultValue;
  const atoms = structureAtoms(structure, (u, e) => getPhysicalRadius(u, e) * sizeFactor);
  // a target per unit, its own instance (the identity): the chains culled and appended in parallel, as units
  const targets: Target[] = [];
  let atomOffset = 0;
  for (const unit of structure.units) {
    const n = unit.elements.length;
    const unitAtoms = atoms.subarray(atomOffset * 4, (atomOffset + n) * 4);
    const t = sphereTarget(atomsSphere(unitAtoms), unitAtoms, [], sizeFactor, {
      atomOffset,
      instanceOffset: targets.length,
      instanceCount: 1,
    });
    targets.push(t);
    atomOffset += n;
  }
  const transforms = new Float32Array(targets.length * 16);
  for (let i = 0; i < targets.length; ++i) transforms.set(Mat4.identity(), i * 16);
  return { atoms, transforms, targets, levels: [] };
}

/** The bounding sphere of `atoms`' centers (`vec4` each: the box's center, its farthest center) */
function atomsSphere(atoms: Float32Array): Sphere3D {
  const min = Vec3.create(Infinity, Infinity, Infinity);
  const max = Vec3.create(-Infinity, -Infinity, -Infinity);
  for (let i = 0; i < atoms.length; i += 4) {
    for (let k = 0; k < 3; ++k) {
      min[k] = Math.min(min[k], atoms[i + k]);
      max[k] = Math.max(max[k], atoms[i + k]);
    }
  }
  const center = Vec3.scale(Vec3(), Vec3.add(Vec3(), min, max), 0.5);
  let radius = 0;
  for (let i = 0; i < atoms.length; i += 4) {
    const dx = atoms[i] - center[0];
    const dy = atoms[i + 1] - center[1];
    const dz = atoms[i + 2] - center[2];
    radius = Math.max(radius, Math.sqrt(dx * dx + dy * dy + dz * dz));
  }
  return Sphere3D.create(center, radius);
}

/**
 * The atoms of every particle's structure target and the particles' transforms: the uniform size theme's radius (1, as
 * the PetWorld scene has it) times the target's size factor
 */
function particleSpheres(list: ParticleList, lodLevels: LodLevels): Spheres {
  const defaults = PD.getDefaultValues(ParticleTargetRepresentationParams);
  const { sizeFactor } = defaults.structure;
  const radius = 1 * sizeFactor;
  const groups = getParticleTargetGroups(list);
  const atoms: Float32Array[] = [];
  const transforms: Float32Array[] = [];
  const targets: Target[] = [];
  let atomOffset = 0;
  let instanceOffset = 0;
  for (let g = 0; g < groups.targetIds.length; ++g) {
    const target = list.targetMapping?.get(groups.targetIds[g]);
    if (target?.kind !== 'structure') continue;
    const { structure } = target;
    const indices = groups.sets[g];
    const a = structureAtoms(structure, () => radius);
    const instanceCount = OrderedSet.size(indices);
    const t = createTargetParticleTransform(
      list,
      indices,
      structure.boundary.sphere,
      defaults.cellSize,
      defaults.batchSize,
      false,
    ).aTransform.ref.value.subarray(0, instanceCount * 16);
    const k = sphereTarget(structure.boundary.sphere, a, lodLevels, sizeFactor, {
      atomOffset,
      instanceOffset,
      instanceCount,
    });
    atoms.push(k.atoms);
    transforms.push(t);
    targets.push(k);
    atomOffset += a.length / 4;
    instanceOffset += instanceCount;
  }
  return { atoms: concat(atoms), transforms: concat(transforms), targets, levels: lodLevels };
}

const CommonWgsl = /* wgsl */ `
struct Params {
    view: mat4x4f,
    proj: mat4x4f,
    size: vec2f,
    maxRadius: f32,
    near: f32,
    lodMin: vec4f,
    lodMax: vec4f,
    levels: u32,
    /** the previous frame's camera (one-phase occlusion: its pyramid's) */
    prevView: mat4x4f,
    prevProj: mat4x4f,
};
/** \`a\`: atom offset, atom count, instance offset, the first id; \`bounds\`: center (object space), radius */
struct Target {
    a: vec4u,
    bounds: vec4f,
    counts: vec4u,
    scales: vec4f,
};
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> atoms: array<vec4f>;
@group(0) @binding(2) var<storage, read> transforms: array<mat4x4f>;

/** the view-space direction through pixel \`x\`, \`y\` (perspective, the eye at the origin) */
fn pixelDirection(x: i32, y: i32) -> vec3f {
    let nx = (f32(x) + 0.5) / P.size.x * 2.0 - 1.0;
    let ny = 1.0 - (f32(y) + 0.5) / P.size.y * 2.0;
    return vec3f((nx + P.proj[2][0]) / P.proj[0][0], (ny + P.proj[2][1]) / P.proj[1][1], -1.0);
}

struct Level { count: u32, scale: f32 };

/** the level of instance transform \`m\` of target \`t\` by the distance of its center, count 0 outside the levels */
fn instanceLevel(t: Target, m: mat4x4f) -> Level {
    if (P.levels == 0u) { return Level(t.a.y, 1.0); }
    let d = -(P.view * (m * vec4f(t.bounds.xyz, 1.0))).z;
    for (var i = 0u; i < P.levels; i++) {
        if (d >= P.lodMin[i] && d < P.lodMax[i]) { return Level(t.counts[i], t.scales[i]); }
    }
    return Level(0u, 1.0);
}

/** the pixel rectangle (center, half size, +1) of the sphere at view-space \`c\` with radius \`r\` */
fn screenRect(c: vec3f, r: f32) -> array<vec2f, 2> {
    let dz = -c.z;
    let q = r / sqrt(dz * dz - r * r);
    let ndc = vec2f((P.proj[0][0] * c.x + P.proj[2][0] * c.z) / dz, (P.proj[1][1] * c.y + P.proj[2][1] * c.z) / dz);
    let center = vec2f((ndc.x * 0.5 + 0.5) * P.size.x, (0.5 - ndc.y * 0.5) * P.size.y);
    let rad = vec2f(q * P.proj[0][0] * 0.5 * P.size.x, q * P.proj[1][1] * 0.5 * P.size.y) + 1.0;
    return array(center, rad);
}

fn offScreen(rect: array<vec2f, 2>) -> bool {
    return any(rect[0] + rect[1] < vec2f(0.0)) || any(rect[0] - rect[1] > P.size);
}

struct Bounds { lo: vec2i, hi: vec2i, rad: vec2f };

/**
 * The pixels whose centers the sphere at view-space \`c\` with radius \`r\` (in front of the eye) can cover, from its
 * exact perspective bounds (the tangents from the eye: \`(c dz ± r sqrt(c² + dz² - r²)) / (dz² - r²)\` per axis);
 * empty where \`lo > hi\` (off screen, or between pixel centers); \`rad\` its half size in pixels
 */
fn pixelBounds(c: vec3f, r: f32) -> Bounds {
    return pixelBoundsWith(c, r, P.proj);
}

/** \`pixelBounds\` with projection \`proj\` (the previous frame's: one-phase occlusion) */
fn pixelBoundsWith(c: vec3f, r: f32, proj: mat4x4f) -> Bounds {
    let P = Params(P.view, proj, P.size, P.maxRadius, P.near, P.lodMin, P.lodMax, P.levels, P.prevView, P.prevProj);
    let dz = -c.z;
    let den = dz * dz - r * r;
    let t = sqrt(c.xy * c.xy + den);
    let u0 = (c.xy * dz - r * t) / den;
    let u1 = (c.xy * dz + r * t) / den;
    let s0 = vec2f(P.proj[0][0], P.proj[1][1]) * u0 - vec2f(P.proj[2][0], P.proj[2][1]);
    let s1 = vec2f(P.proj[0][0], P.proj[1][1]) * u1 - vec2f(P.proj[2][0], P.proj[2][1]);
    // screen y grows down: ndc y's maximum is the top
    let x = (vec2f(s0.x, s1.x) * 0.5 + 0.5) * P.size.x;
    let y = (0.5 - vec2f(s1.y, s0.y) * 0.5) * P.size.y;
    let lo = vec2i(i32(ceil(x.x - 0.5)), i32(ceil(y.x - 0.5)));
    let hi = vec2i(i32(floor(x.y - 0.5)), i32(floor(y.y - 0.5)));
    return Bounds(max(lo, vec2i(0)), min(hi, vec2i(P.size) - 1), vec2f(x.y - x.x, y.y - y.x) * 0.5);
}

/** the Hi-Z pyramid: level \`l\` the farthest depth key (the smallest, 0 uncovered) of tiles of \`8 << l\` pixels */
const HizLevels = 5u;
fn hizDims(level: u32) -> vec2u {
    let tile = 8u << level;
    return (vec2u(P.size) + tile - 1u) / tile;
}
fn hizOffset(level: u32) -> u32 {
    var o = 0u;
    for (var l = 0u; l < level; l++) {
        let d = hizDims(l);
        o += d.x * d.y;
    }
    return o;
}

/** the front hit's ray parameter along \`d\` of the sphere at view-space \`c\` with radius \`r\`, -1 on a miss */
fn hitT(d: vec3f, c: vec3f, r: f32) -> f32 {
    let b = dot(d, c);
    let dd = dot(d, d);
    let disc = b * b - dd * (dot(c, c) - r * r);
    if (disc < 0.0) { return -1.0; }
    return (b - sqrt(disc)) / dd;
}
`;

/**
 * A kernel per target (`T`): a workgroup per chunk of `ChunkAtoms` atoms of an instance (chunks along `x`,
 * instances along `y` and `z`), the instance's level and frustum test, then per sphere of the level `perSphere`
 * with `rad` (the projected radius in pixels, +1) and the loop over the sphere's pixel rectangle, `perPixel` with
 * `p` (the pixel), `key` (the depth) and `id` (1-based)
 */
function sphereKernel(bindings: string, perSphere: string, perPixel: string, enable = '') {
  return /* wgsl */ `
${enable}
${CommonWgsl}
@group(0) @binding(3) var<uniform> T: Target;
${bindings}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_index) lid: u32) {
    let instance = wid.y + wid.z * 65535u;
    let m = transforms[T.a.z + instance];
    let level = instanceLevel(T, m);
    let first = wid.x * ${ChunkAtoms}u;
    if (first >= level.count) { return; }
    let vb = (P.view * (m * vec4f(T.bounds.xyz, 1.0))).xyz;
    if (-vb.z - T.bounds.w > P.near && offScreen(screenRect(vb, T.bounds.w))) { return; }
    let last = min(first + ${ChunkAtoms}u, level.count);
    for (var atom = first + lid; atom < last; atom += 64u) {
        let s = atoms[T.a.x + atom];
        let scale = level.scale;
        let id = T.a.w + instance * T.a.y + atom + 1u;
        ${sphereBody(perSphere, perPixel, 'continue')}
    }
}
`;
}

/**
 * Sphere `s` (object space) of instance transform `m` at `scale` with `id`: frustum test, `perSphere`, the size limit,
 * then the loop over its pixel rectangle with `perPixel`; `exit` leaves the sphere (`continue` in a loop, or `return`)
 */
function sphereBody(perSphere: string, perPixel: string, exit: string) {
  return /* wgsl */ `
        let c = (P.view * (m * vec4f(s.xyz, 1.0))).xyz;
        let r = s.w * scale;
        if (-c.z - r < P.near) { ${exit}; }
        let rect = screenRect(c, r);
        if (offScreen(rect)) { ${exit}; }
        let center = rect[0];
        let rad = rect[1];
        ${perSphere}
        if (max(rad.x, rad.y) > P.maxRadius) { ${exit}; }
        let width = u32(P.size.x);
        let lo = max(vec2i(floor(center - rad)), vec2i(0));
        let hi = min(vec2i(ceil(center + rad)), vec2i(P.size) - 1);
        for (var y = lo.y; y <= hi.y; y++) {
            for (var x = lo.x; x <= hi.x; x++) {
                let t = hitT(pixelDirection(x, y), c, r);
                if (t <= 0.0) { continue; }
                // positive floats order as their bits: nearer is larger
                let key = bitcast<u32>(1.0 / t);
                let p = u32(y) * width + u32(x);
                ${perPixel}
            }
        }`;
}

/** the compacted variant's bindings: targets, each instance's target, the chunk counter, the chunks */
const CompactBindings = /* wgsl */ `
@group(0) @binding(8) var<storage, read> targets: array<Target>;
@group(0) @binding(9) var<storage, read> instanceTarget: array<u32>;
@group(0) @binding(11) var<storage, read_write> chunks: array<vec4u>;
`;

/**
 * The compacted variant's first pass: an invocation per instance (every target) culls it by the frustum, picks its
 * level and appends a chunk per 64 of the level's atoms (instance, first atom, end, scale). Occlusion in two phases
 * (the camera may move): `phase1` appends only the instances visible in the previous frame (`visible`); then the Hi-Z
 * is built from their pixels; `phase2` tests every instance against it (its farthest depth per tile, at the level
 * where the instance's rectangle spans at most 3 tiles a side), appends those newly visible and updates `visible`
 */
function cullKernel(mode: 'all' | 'one' | 'phase1' | 'phase2') {
  return /* wgsl */ `
${CommonWgsl}
${CompactBindings}
@group(0) @binding(10) var<storage, read_write> counter: array<atomic<u32>, 4>;
${mode === 'phase2' || mode === 'one' ? '@group(0) @binding(12) var<storage, read> hiz: array<u32>;' : ''}
${mode === 'phase1' || mode === 'phase2' ? `@group(0) @binding(14) var<storage, ${mode === 'phase2' ? 'read_write' : 'read'}> visible: array<u32>;` : ''}

/** whether the instance's sphere at view-space \`vb\` with radius \`r\` is behind the Hi-Z (its tiles' farthest depth) */
fn occluded(vb: vec3f, r: f32, b: Bounds) -> bool {
    ${
      mode === 'phase2' || mode === 'one'
        ? /* wgsl */ `
    let span = u32(max(b.hi.x - b.lo.x, b.hi.y - b.lo.y)) + 1u;
    var l = 0u;
    while (l < HizLevels - 1u && span > (16u << l)) { l++; }
    let tile = 8u << l;
    let tlo = vec2u(b.lo) / tile;
    let thi = vec2u(b.hi) / tile;
    if (any(thi - tlo >= vec2u(3u))) { return false; }
    let dims = hizDims(l);
    let offset = hizOffset(l);
    var far = 0xffffffffu;
    for (var y = tlo.y; y <= thi.y; y++) {
        for (var x = tlo.x; x <= thi.x; x++) { far = min(far, hiz[offset + y * dims.x + x]); }
    }
    return bitcast<u32>(1.0 / (-vb.z - r)) < far;`
        : 'return false;'
    }
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3u, @builtin(num_workgroups) groups: vec3u) {
    let instance = gid.x + gid.y * groups.x * 64u;
    if (instance >= arrayLength(&instanceTarget)) { return; }
    let t = targets[instanceTarget[instance]];
    let m = transforms[instance];
    let level = instanceLevel(t, m);
    let vb = (P.view * (m * vec4f(t.bounds.xyz, 1.0))).xyz;
    var inView = level.count > 0u;
    var occ = false;
    if (inView && -vb.z - t.bounds.w > P.near) {
        let b = pixelBounds(vb, t.bounds.w);
        inView = all(b.lo <= b.hi);
        ${
          mode === 'one'
            ? /* wgsl */ `
        // against the previous frame's pyramid, with its camera
        let vp = (P.prevView * (m * vec4f(t.bounds.xyz, 1.0))).xyz;
        if (inView && -vp.z - t.bounds.w > P.near) {
            let bp = pixelBoundsWith(vp, t.bounds.w, P.prevProj);
            if (all(bp.lo <= bp.hi)) { occ = occluded(vp, t.bounds.w, bp); }
        }`
            : 'if (inView) { occ = occluded(vb, t.bounds.w, b); }'
        }
    }
    ${
      mode === 'phase1'
        ? 'if (!inView || visible[instance] == 0u) { return; }'
        : mode === 'phase2'
          ? /* wgsl */ `
    let was = visible[instance];
    visible[instance] = select(0u, 1u, inView && !occ);
    if (!inView || occ || was != 0u) { return; }`
          : mode === 'one'
            ? 'if (!inView || occ) { return; }'
            : 'if (!inView) { return; }'
    }
    let n = (level.count + 63u) / 64u;
    let base = atomicAdd(&counter[0], n);
    let capacity = arrayLength(&chunks);
    for (var k = 0u; k < n && base + k < capacity; k++) {
        chunks[base + k] = vec4u(instance, k * 64u, min(k * 64u + 64u, level.count), bitcast<u32>(level.scale));
    }
}
`;
}

/** the Hi-Z pyramid's level `L.x`, from the visibility buffer (level 0) or the level below */
const HizBuildWgsl = /* wgsl */ `
${CommonWgsl}
@group(0) @binding(4) var<storage, read> vis: array<vec2u>;
@group(0) @binding(12) var<storage, read_write> hiz: array<u32>;
@group(0) @binding(13) var<uniform> L: vec4u;

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    let level = L.x;
    let dims = hizDims(level);
    if (gid.x >= dims.x || gid.y >= dims.y) { return; }
    var m = 0xffffffffu;
    if (level == 0u) {
        let width = u32(P.size.x);
        let lo = gid.xy * 8u;
        let hi = min(lo + 8u, vec2u(P.size));
        for (var y = lo.y; y < hi.y; y++) {
            for (var x = lo.x; x < hi.x; x++) { m = min(m, vis[y * width + x].y); }
        }
    } else {
        let below = hizDims(level - 1u);
        let offset = hizOffset(level - 1u);
        let lo = gid.xy * 2u;
        let hi = min(lo + 2u, below);
        for (var y = lo.y; y < hi.y; y++) {
            for (var x = lo.x; x < hi.x; x++) { m = min(m, hiz[offset + y * below.x + x]); }
        }
    }
    hiz[hizOffset(level) + gid.y * dims.x + gid.x] = m;
}
`;

/** the indirect dispatch of the chunks (65535 a row) */
const ArgsWgsl = /* wgsl */ `
@group(0) @binding(0) var<storage, read> counter: array<u32, 4>;
@group(0) @binding(1) var<storage, read_write> args: array<u32, 4>;
@group(0) @binding(2) var<uniform> capacity: vec4u;
@compute @workgroup_size(1)
fn main() {
    let n = min(counter[0], capacity.x);
    args[0] = min(n, 65535u);
    args[1] = (n + 65534u) / 65535u;
    args[2] = 1u;
}
`;

/** the compacted variant's second pass: a workgroup per chunk, an invocation per sphere */
function compactKernel(bindings: string, perPixel: string, enable = '') {
  return /* wgsl */ `
${enable}
${CommonWgsl}
${CompactBindings.replace('read_write> chunks', 'read> chunks')}
@group(0) @binding(10) var<storage, read> counter: array<u32, 4>;
${bindings}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_index) lid: u32) {
    let chunk = wid.x + wid.y * 65535u;
    if (chunk >= min(counter[0], arrayLength(&chunks))) { return; }
    let e = chunks[chunk];
    let atom = e.y + lid;
    if (atom >= e.z) { return; }
    let t = targets[instanceTarget[e.x]];
    let m = transforms[e.x];
    let s = atoms[t.a.x + atom];
    let scale = bitcast<f32>(e.w);
    let id = t.a.w + (e.x - t.a.z) * t.a.y + atom + 1u;
    ${sphereBody('', perPixel, 'return')}
}
`;
}

/**
 * The compacted second pass with the exact pixel bounds (`pixelBounds`) and, with `tile`, the chunk's depths first
 * reduced in a workgroup tile of `TileSize`² pixels at its instance's rectangle (where that fits): only each pixel's
 * nearest of the chunk writes the visibility buffer
 */
function tightKernel(tile: boolean) {
  const perPixel = /* wgsl */ `
                let t = hitT(pixelDirection(x, y), c, r);
                if (t <= 0.0) { continue; }
                let key = bitcast<u32>(1.0 / t);
                let p = u32(y) * width + u32(x);`;
  return /* wgsl */ `
${Enable64}
${CommonWgsl}
${CompactBindings.replace('read_write> chunks', 'read> chunks')}
@group(0) @binding(10) var<storage, read> counter: array<u32, 4>;
${Vis64}
${tile ? `var<workgroup> tile: array<atomic<u32>, ${TileSize * TileSize}>;` : ''}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_index) lid: u32) {
    ${tile ? `for (var k = lid; k < ${TileSize * TileSize}u; k += 64u) { atomicStore(&tile[k], 0u); }\n    workgroupBarrier();` : ''}
    let chunk = wid.x + wid.y * 65535u;
    var ok = chunk < min(counter[0], arrayLength(&chunks));
    var e = vec4u(0u);
    if (ok) { e = chunks[chunk]; }
    let atom = e.y + lid;
    ok = ok && atom < e.z;
    let t = targets[instanceTarget[e.x]];
    let m = transforms[e.x];
    var c = vec3f(0.0);
    var r = 0.0;
    var b = Bounds(vec2i(0), vec2i(-1), vec2f(0.0));
    if (ok) {
        let s = atoms[t.a.x + atom];
        c = (P.view * (m * vec4f(s.xyz, 1.0))).xyz;
        r = s.w * bitcast<f32>(e.w);
        ok = -c.z - r > P.near;
        if (ok) {
            b = pixelBounds(c, r);
            ok = all(b.lo <= b.hi) && max(b.rad.x, b.rad.y) <= P.maxRadius;
        }
    }
    let id = t.a.w + (e.x - t.a.z) * t.a.y + atom + 1u;
    let width = u32(P.size.x);
    ${
      tile
        ? /* wgsl */ `
    // the instance's rectangle: the tile's origin where it fits
    let vb = (P.view * (m * vec4f(t.bounds.xyz, 1.0))).xyz;
    var origin = vec2i(0);
    var useTile = false;
    if (-vb.z - t.bounds.w > P.near) {
        let ib = pixelBounds(vb, t.bounds.w);
        origin = ib.lo;
        useTile = all(ib.hi - ib.lo < vec2i(${TileSize}));
    }
    if (ok) {
        for (var y = b.lo.y; y <= b.hi.y; y++) {
            for (var x = b.lo.x; x <= b.hi.x; x++) {${perPixel}
                let l = vec2i(x, y) - origin;
                if (useTile && all(l >= vec2i(0)) && all(l < vec2i(${TileSize}))) {
                    atomicMax(&tile[l.y * ${TileSize} + l.x], key);
                } else {
                    atomicStoreMax(&vis[p], vec2u(id, key));
                }
            }
        }
    }
    workgroupBarrier();
    if (ok && useTile) {
        for (var y = b.lo.y; y <= b.hi.y; y++) {
            for (var x = b.lo.x; x <= b.hi.x; x++) {${perPixel}
                let l = vec2i(x, y) - origin;
                if (all(l >= vec2i(0)) && all(l < vec2i(${TileSize})) && atomicLoad(&tile[l.y * ${TileSize} + l.x]) == key) {
                    atomicStoreMax(&vis[p], vec2u(id, key));
                }
            }
        }
    }`
        : /* wgsl */ `
    if (ok) {
        for (var y = b.lo.y; y <= b.hi.y; y++) {
            for (var x = b.lo.x; x <= b.hi.x; x++) {${perPixel}
                atomicStoreMax(&vis[p], vec2u(id, key));
            }
        }
    }`
    }
}
`;
}

/** the workgroup tile's side in pixels */
const TileSize = 32;

/** atoms per workgroup of an instance (16 per invocation) */
const ChunkAtoms = 1024;

const Vis64 = '@group(0) @binding(4) var<storage, read_write> vis: array<atomic<vec2u>>;';
const Depth32 = '@group(0) @binding(5) var<storage, read_write> depth32: array<atomic<u32>>;';
const Id32 = '@group(0) @binding(6) var<storage, read_write> id32: array<u32>;';
const Bins = '@group(0) @binding(7) var<storage, read_write> bins: array<atomic<u32>, 4>;';
const Enable64 = 'enable atomic_vec2u_min_max;';

const Kernels = {
  raster64: sphereKernel(Vis64, '', 'atomicStoreMax(&vis[p], vec2u(id, key));', Enable64),
  raster64Early: sphereKernel(
    `${Vis64}\n${Depth32}`,
    '',
    'if (atomicLoad(&depth32[p]) > key) { continue; }\natomicMax(&depth32[p], key);\natomicStoreMax(&vis[p], vec2u(id, key));',
    Enable64,
  ),
  depth32: sphereKernel(Depth32, '', 'atomicMax(&depth32[p], key);'),
  id32: sphereKernel(`${Depth32}\n${Id32}`, '', 'if (atomicLoad(&depth32[p]) == key) { id32[p] = id; }'),
  raster64Compact: compactKernel(Vis64, 'atomicStoreMax(&vis[p], vec2u(id, key));', Enable64),
  rasterTight: tightKernel(false),
  rasterTile: tightKernel(true),
  cull: cullKernel('all'),
  cullPhase1: cullKernel('phase1'),
  cullPhase2: cullKernel('phase2'),
  cullOne: cullKernel('one'),
  hizBuild: HizBuildWgsl,
  args: ArgsWgsl,
  // spheres in the frustum by projected radius (+1) in pixels: ≤2, ≤4, ≤8, above
  count: sphereKernel(
    Bins,
    'let px = max(rad.x, rad.y);\natomicAdd(&bins[select(select(select(3u, 2u, px <= 8.0), 1u, px <= 4.0), 0u, px <= 2.0)], 1u);\ncontinue;',
    '',
  ),
};

const ClearWgsl = /* wgsl */ `
@group(0) @binding(0) var<storage, read_write> words: array<u32>;
@group(0) @binding(1) var<uniform> count: vec4u;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3u, @builtin(num_workgroups) n: vec3u) {
    let i = gid.x + gid.y * n.x * 64u;
    if (i < count.x) { words[i] = 0u; }
}
`;

const ResolveWgsl = /* wgsl */ `
${CommonWgsl}
@group(0) @binding(3) var<storage, read> targets: array<Target>;
@group(0) @binding(4) var<storage, read> vis: array<vec2u>;

@vertex
fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let p = array(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(p[i], 0.0, 1.0);
}

struct Out { @location(0) color: vec4f, @builtin(frag_depth) depth: f32 };

@fragment
fn fs(@builtin(position) pos: vec4f) -> Out {
    let x = i32(pos.x);
    let y = i32(pos.y);
    let v = vis[u32(y) * u32(P.size.x) + u32(x)];
    if (v.x == 0u) { discard; }
    let id = v.x - 1u;
    // the last target whose first id is at or below
    var lo = 0u;
    var hi = arrayLength(&targets) - 1u;
    while (lo < hi) {
        let mid = (lo + hi + 1u) / 2u;
        if (targets[mid].a.w <= id) { lo = mid; } else { hi = mid - 1u; }
    }
    let t = targets[lo];
    let local = id - t.a.w;
    let s = atoms[t.a.x + local % t.a.y];
    let m = transforms[t.a.z + local / t.a.y];
    let r = s.w * instanceLevel(t, m).scale;
    let c = (P.view * (m * vec4f(s.xyz, 1.0))).xyz;
    let d = pixelDirection(x, y);
    let hit = d * max(hitT(d, c, r), 0.0);
    let normal = normalize(hit - c);
    let shade = max(dot(normal, normalize(vec3f(0.3, 0.5, 1.0))), 0.0) * 0.8 + 0.2;
    let clip = P.proj * vec4f(hit, 1.0);
    return Out(vec4f(vec3f(shade), 1.0), clamp(clip.z / clip.w, 0.0, 1.0));
}
`;

/** the resolve's color (R unlimited) as a PNG data URL, beside the impostors' image (`--images`) */
async function textureImage(device: GPUDevice, texture: GPUTexture, width: number, height: number) {
  const bytesPerRow = Math.ceil((width * 4) / 256) * 256;
  const buffer = device.createBuffer({
    size: bytesPerRow * height,
    usage: BufferUsage.MAP_READ | BufferUsage.COPY_DST,
  });
  const encoder = device.createCommandEncoder();
  encoder.copyTextureToBuffer({ texture }, { buffer, bytesPerRow }, [width, height]);
  device.queue.submit([encoder.finish()]);
  await buffer.mapAsync(MapMode.READ);
  const src = new Uint8Array(buffer.getMappedRange());
  const image = new ImageData(width, height);
  for (let y = 0; y < height; ++y) {
    image.data.set(src.subarray(y * bytesPerRow, y * bytesPerRow + width * 4), y * width * 4);
  }
  buffer.unmap();
  buffer.destroy();
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  canvas.getContext('2d')!.putImageData(image, 0, 0);
  return canvas.toDataURL('image/png');
}

async function readBuffer(device: GPUDevice, src: GPUBuffer, size: number) {
  const buffer = device.createBuffer({ size, usage: BufferUsage.MAP_READ | BufferUsage.COPY_DST });
  const encoder = device.createCommandEncoder();
  encoder.copyBufferToBuffer(src, 0, buffer, 0, size);
  device.queue.submit([encoder.finish()]);
  await buffer.mapAsync(MapMode.READ);
  const out = new Uint32Array(buffer.getMappedRange().slice(0));
  buffer.unmap();
  buffer.destroy();
  return out;
}

/** GPU times (ms) of `passes` per frame over `frames`, by timestamp queries; medians */
async function timeFrames(
  device: GPUDevice,
  frames: number,
  passes: string[],
  encode: (encoder: GPUCommandEncoder, timestamps: (pass: number) => GPUComputePassTimestampWrites) => void,
) {
  const querySet = device.createQuerySet({ type: 'timestamp', count: passes.length * 2 });
  const resolve = device.createBuffer({
    size: passes.length * 16,
    usage: BufferUsage.QUERY_RESOLVE | BufferUsage.COPY_SRC,
  });
  const read = device.createBuffer({
    size: passes.length * 16,
    usage: BufferUsage.MAP_READ | BufferUsage.COPY_DST,
  });
  const times: number[][] = passes.map(() => []);
  for (let f = 0; f < frames; ++f) {
    const encoder = device.createCommandEncoder();
    encode(encoder, (k) => ({ querySet, beginningOfPassWriteIndex: 2 * k, endOfPassWriteIndex: 2 * k + 1 }));
    encoder.resolveQuerySet(querySet, 0, passes.length * 2, resolve, 0);
    encoder.copyBufferToBuffer(resolve, 0, read, 0, passes.length * 16);
    device.queue.submit([encoder.finish()]);
    await read.mapAsync(MapMode.READ);
    const t = new BigUint64Array(read.getMappedRange().slice(0));
    read.unmap();
    for (let k = 0; k < passes.length; ++k) times[k].push(Number(t[2 * k + 1] - t[2 * k]) / 1e6);
  }
  querySet.destroy();
  resolve.destroy();
  read.destroy();
  const median = (a: number[]) => a.sort((x, y) => x - y)[Math.floor(a.length / 2)];
  return Object.fromEntries(passes.map((p, k) => [p, median(times[k])]));
}

async function runPrototype(id: string, graphics: GraphicsContext, input: RasterInput): Promise<PathRun | string> {
  if (!isWebGPUContext(graphics)) return 'WebGPU only';
  const { device } = graphics;
  if (!device.features.has('atomic-vec2u-min-max' as GPUFeatureName)) {
    return 'no atomic-vec2u-min-max (Chrome 156: --browser-arg --enable-webgpu-developer-features)';
  }
  if (!device.features.has('timestamp-query')) return 'no timestamp-query';

  // the impostors on the same scene and camera (their opaque pass with --gpu-timing)
  // a profiler's run (`BenchmarkEnv.rasterProfile`): the impostors' frames alone, without end
  if (BenchmarkEnv.rasterProfile === 'impostors') for (;;) await renderFrames(id, graphics, input);
  const impostors: PathRun =
    BenchmarkEnv.rasterProfile === 'raster' ? { stages: [], triangles: 0 } : await renderFrames(id, graphics, input);
  const canvas3d = await getScene(id, graphics, input);
  const { camera } = canvas3d;
  // its matrices of the case's view (zoomed): updated by the impostors' frames, not drawn by a profiler's raster run
  camera.update();
  if (camera.state.mode !== 'perspective') return 'perspective cameras only';
  const direction = () => Vec3.normalize(Vec3(), Vec3.sub(Vec3(), camera.state.position, camera.state.target));
  const startDirection = direction();
  const { width, height } = input.size ?? Size;
  const pixels = width * height;

  const spheres = input.spheres();
  let count = 0;
  // a `Target` each (16 words): packed for the resolve, at the binding offset alignment as the kernels' uniforms
  const targetTable = new ArrayBuffer(spheres.targets.length * 64);
  const targetUniforms = new ArrayBuffer(spheres.targets.length * 256);
  spheres.targets.forEach((t, k) => {
    for (const [buffer, offset] of [
      [targetTable, k * 64],
      [targetUniforms, k * 256],
    ] as const) {
      const u = new Uint32Array(buffer, offset, 16);
      const f = new Float32Array(buffer, offset, 16);
      u.set([t.atomOffset, t.atomCount, t.instanceOffset, count]);
      f.set([...t.center, t.radius], 4);
      u.set(t.counts.slice(0, 4), 8);
      f.set(t.scales.slice(0, 4), 12);
    }
    count += t.atomCount * t.instanceCount;
  });
  if (count >= 2 ** 32 - 1) return `${count.toLocaleString('en')} spheres: above 32-bit ids`;
  const storage = (data: Float32Array | Uint32Array | ArrayBuffer, usage: number = BufferUsage.STORAGE) => {
    const b = device.createBuffer({ size: Math.max(16, data.byteLength), usage: usage | BufferUsage.COPY_DST });
    device.queue.writeBuffer(b, 0, data);
    return b;
  };
  const atomBuffer = storage(spheres.atoms);
  const transformBuffer = storage(spheres.transforms);
  const targetBuffer = storage(targetTable);
  const targetUniformBuffer = storage(targetUniforms, BufferUsage.UNIFORM);
  if (spheres.levels.length > 4) return 'more than 4 levels';
  const vis = device.createBuffer({ size: pixels * 8, usage: BufferUsage.STORAGE | BufferUsage.COPY_SRC });
  const depth32 = device.createBuffer({ size: pixels * 4, usage: BufferUsage.STORAGE });
  const id32 = device.createBuffer({ size: pixels * 4, usage: BufferUsage.STORAGE });
  const bins = device.createBuffer({ size: 16, usage: BufferUsage.STORAGE | BufferUsage.COPY_SRC });
  // the compacted variant: each instance's target, the chunk counter, the chunks (as many as a binding holds, at most
  // 256 MB), the indirect dispatch
  const instanceTarget = new Uint32Array(spheres.transforms.length / 16);
  spheres.targets.forEach((t, k) => instanceTarget.fill(k, t.instanceOffset, t.instanceOffset + t.instanceCount));
  const instanceTargetBuffer = storage(instanceTarget);
  const counter = device.createBuffer({ size: 16, usage: BufferUsage.STORAGE | BufferUsage.COPY_SRC });
  const chunkCapacity = Math.floor(Math.min(device.limits.maxStorageBufferBindingSize, 2 ** 28) / 16);
  const chunks = device.createBuffer({ size: chunkCapacity * 16, usage: BufferUsage.STORAGE });
  const args = device.createBuffer({ size: 16, usage: BufferUsage.STORAGE | BufferUsage.INDIRECT });
  const capacityUniform = storage(new Uint32Array([chunkCapacity, 0, 0, 0]), BufferUsage.UNIFORM);
  // the Hi-Z pyramid (`HizLevels` levels of tiles of 8 << l pixels), a uniform per level for its build
  const hizLevels = Array.from({ length: 5 }, (_, l) => {
    const tile = 8 << l;
    return { l, x: Math.ceil(width / tile), y: Math.ceil(height / tile) };
  });
  const hiz = device.createBuffer({
    size: hizLevels.reduce((s, d) => s + d.x * d.y, 0) * 4,
    usage: BufferUsage.STORAGE,
  });
  const hizLevelUniforms = hizLevels.map((d) => storage(new Uint32Array([d.l, 0, 0, 0]), BufferUsage.UNIFORM));
  // one-phase occlusion's pyramid: of the previous frame's visibility buffer, none at first (0: nothing occludes)
  const hizOne = device.createBuffer({
    size: hizLevels.reduce((s, d) => s + d.x * d.y, 0) * 4,
    usage: BufferUsage.STORAGE,
  });
  // each instance visible in the previous frame (two-phase occlusion), none at first
  const visibleBuffer = storage(new Uint32Array(instanceTarget.length));
  const color = device.createTexture({
    size: [width, height],
    format: 'rgba8unorm',
    usage: TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC,
  });
  const depth = device.createTexture({
    size: [width, height],
    format: 'depth32float',
    usage: TextureUsage.RENDER_ATTACHMENT,
  });

  /** the camera's matrices into `uniforms` (its first 32 floats) */
  const writeCamera = (uniforms: GPUBuffer) => {
    const data = new Float32Array(32);
    data.set(camera.view, 0);
    data.set(camera.projection, 16);
    device.queue.writeBuffer(uniforms, 0, data);
  };
  /** the camera's matrices as the previous frame's (`prevView`, `prevProj`, from float 48) */
  const writePreviousCamera = (uniforms: GPUBuffer) => {
    const data = new Float32Array(32);
    data.set(camera.view, 0);
    data.set(camera.projection, 16);
    device.queue.writeBuffer(uniforms, 48 * 4, data);
  };
  /** the camera one step on (`rotate`): the frame's previous camera and its own */
  const step = (uniforms: GPUBuffer) => {
    writePreviousCamera(uniforms);
    rotateCamera(canvas3d);
    camera.update();
    writeCamera(uniforms);
  };
  const params = (maxRadius: number) => {
    const data = new Float32Array(80);
    data.set(camera.view, 0);
    data.set(camera.projection, 16);
    data[32] = width;
    data[33] = height;
    data[34] = maxRadius;
    data[35] = camera.near;
    spheres.levels.forEach((l, i) => {
      data[36 + i] = l.minDistance;
      data[40 + i] = l.maxDistance;
    });
    new Uint32Array(data.buffer)[44] = spheres.levels.length;
    data.set(camera.view, 48);
    data.set(camera.projection, 64);
    return storage(data, BufferUsage.UNIFORM);
  };
  const pipeline = (code: string) =>
    device.createComputePipeline({
      layout: 'auto',
      compute: { module: device.createShaderModule({ code }), entryPoint: 'main' },
    });
  const dispatch = (pass: GPUComputePassEncoder, n: number) => {
    const x = Math.min(n, 65535);
    pass.dispatchWorkgroups(x, Math.ceil(n / x));
  };
  const clearPipeline = pipeline(ClearWgsl);
  const clears = [
    { buffer: vis, words: pixels * 2 },
    { buffer: depth32, words: pixels },
    { buffer: id32, words: pixels },
    { buffer: bins, words: 4 },
    { buffer: counter, words: 4 },
  ].map((c) => {
    const count = storage(new Uint32Array([c.words, 0, 0, 0]), BufferUsage.UNIFORM);
    const group = device.createBindGroup({
      layout: clearPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: c.buffer } },
        { binding: 1, resource: { buffer: count } },
      ],
    });
    return { ...c, count, group };
  });
  const clear = (encoder: GPUCommandEncoder, which: number[], timestampWrites?: GPUComputePassTimestampWrites) => {
    const pass = encoder.beginComputePass({ timestampWrites });
    pass.setPipeline(clearPipeline);
    for (const k of which) {
      pass.setBindGroup(0, clears[k].group);
      dispatch(pass, Math.ceil(clears[k].words / 64));
    }
    pass.end();
  };

  const kernels = Object.fromEntries(Object.entries(Kernels).map(([k, code]) => [k, pipeline(code)])) as {
    [k in keyof typeof Kernels]: GPUComputePipeline;
  };
  /** a bind group per target (its uniform at its offset) */
  const kernelGroups = (p: GPUComputePipeline, uniforms: GPUBuffer, extra: { binding: number; buffer: GPUBuffer }[]) =>
    spheres.targets.map((_, k) =>
      device.createBindGroup({
        layout: p.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: uniforms } },
          { binding: 1, resource: { buffer: atomBuffer } },
          { binding: 2, resource: { buffer: transformBuffer } },
          { binding: 3, resource: { buffer: targetUniformBuffer, offset: k * 256, size: 64 } },
          ...extra.map((e) => ({ binding: e.binding, resource: { buffer: e.buffer } })),
        ],
      }),
    );
  const kernelPass = (
    encoder: GPUCommandEncoder,
    p: GPUComputePipeline,
    groups: GPUBindGroup[],
    timestampWrites?: GPUComputePassTimestampWrites,
  ) => {
    const pass = encoder.beginComputePass({ timestampWrites });
    pass.setPipeline(p);
    spheres.targets.forEach((t, k) => {
      pass.setBindGroup(0, groups[k]);
      pass.dispatchWorkgroups(
        Math.ceil(t.atomCount / ChunkAtoms),
        Math.min(t.instanceCount, 65535),
        Math.ceil(t.instanceCount / 65535),
      );
    });
    pass.end();
  };
  const resolveModule = device.createShaderModule({ code: ResolveWgsl });
  const resolvePipeline = device.createRenderPipeline({
    layout: 'auto',
    vertex: { module: resolveModule, entryPoint: 'vs' },
    fragment: { module: resolveModule, entryPoint: 'fs', targets: [{ format: 'rgba8unorm' }] },
    depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'always' },
  });

  const group = (p: GPUComputePipeline, entries: [number, GPUBuffer][]) =>
    device.createBindGroup({
      layout: p.getBindGroupLayout(0),
      entries: entries.map(([binding, buffer]) => ({ binding, resource: { buffer } })),
    });
  const argsGroup = group(kernels.args, [
    [0, counter],
    [1, args],
    [2, capacityUniform],
  ]);
  const instanceGroups = Math.ceil(instanceTarget.length / 64);

  const notes: string[] = [];
  const f = (x: number) => x.toFixed(2);
  const n = (x: number) => x.toLocaleString('en');
  let binCounts: Uint32Array | undefined;
  let chunkCount = 0;
  for (const maxRadius of [2, 4, 8, 1e9]) {
    if (BenchmarkEnv.rasterProfile === 'raster' && maxRadius < 1e8) continue;
    const uniforms = params(maxRadius);
    const g64 = kernelGroups(kernels.raster64, uniforms, [{ binding: 4, buffer: vis }]);
    const gEarly = kernelGroups(kernels.raster64Early, uniforms, [
      { binding: 4, buffer: vis },
      { binding: 5, buffer: depth32 },
    ]);
    const gDepth = kernelGroups(kernels.depth32, uniforms, [{ binding: 5, buffer: depth32 }]);
    const gId = kernelGroups(kernels.id32, uniforms, [
      { binding: 5, buffer: depth32 },
      { binding: 6, buffer: id32 },
    ]);
    if (!binCounts) {
      const encoder = device.createCommandEncoder();
      clear(encoder, [3]);
      kernelPass(encoder, kernels.count, kernelGroups(kernels.count, uniforms, [{ binding: 7, buffer: bins }]));
      device.queue.submit([encoder.finish()]);
      binCounts = await readBuffer(device, bins, 16);
    }
    const resolveGroup = device.createBindGroup({
      layout: resolvePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: uniforms } },
        { binding: 1, resource: { buffer: atomBuffer } },
        { binding: 2, resource: { buffer: transformBuffer } },
        { binding: 3, resource: { buffer: targetBuffer } },
        { binding: 4, resource: { buffer: vis } },
      ],
    });
    const compact: [number, GPUBuffer][] = [
      [0, uniforms],
      [1, atomBuffer],
      [2, transformBuffer],
      [8, targetBuffer],
      [9, instanceTargetBuffer],
      [10, counter],
      [11, chunks],
    ];
    const cullGroup = group(
      kernels.cull,
      compact.filter(([b]) => b !== 1),
    );
    const compactGroup = group(kernels.raster64Compact, [...compact, [4, vis]]);
    const tightGroup = group(kernels.rasterTight, [...compact, [4, vis]]);
    const tileGroup = group(kernels.rasterTile, [...compact, [4, vis]]);
    const cullPhase1Group = group(kernels.cullPhase1, [...compact.filter(([b]) => b !== 1), [14, visibleBuffer]]);
    const cullPhase2Group = group(kernels.cullPhase2, [
      ...compact.filter(([b]) => b !== 1),
      [12, hiz],
      [14, visibleBuffer],
    ]);
    const hizGroups = hizLevelUniforms.map((u) =>
      group(kernels.hizBuild, [
        [0, uniforms],
        [4, vis],
        [12, hiz],
        [13, u],
      ]),
    );
    const hizOneGroups = hizLevelUniforms.map((u) =>
      group(kernels.hizBuild, [
        [0, uniforms],
        [4, vis],
        [12, hizOne],
        [13, u],
      ]),
    );
    const cullOneGroup = group(kernels.cullOne, [...compact.filter(([b]) => b !== 1), [12, hizOne]]);
    const buildPyramid = (
      encoder: GPUCommandEncoder,
      groups: GPUBindGroup[],
      timestampWrites?: GPUComputePassTimestampWrites,
    ) => {
      const build = encoder.beginComputePass({ timestampWrites });
      build.setPipeline(kernels.hizBuild);
      hizLevels.forEach((d, l) => {
        build.setBindGroup(0, groups[l]);
        build.dispatchWorkgroups(Math.ceil(d.x / 8), Math.ceil(d.y / 8));
      });
      build.end();
    };
    /**
     * One-phase occlusion (exact bounds), as the impostors': every instance tested against the previous frame's pyramid
     * with its camera, then the pyramid of this frame's for the next
     */
    const onePhase = (encoder: GPUCommandEncoder, ts: (GPUComputePassTimestampWrites | undefined)[] = []) => {
      clear(encoder, [0, 4]);
      cullPass(encoder, kernels.cullOne, cullOneGroup, ts[0]);
      compactPass(encoder, kernels.rasterTight, tightGroup, ts[1]);
      buildPyramid(encoder, hizOneGroups, ts[2]);
    };
    const compactPass = (
      encoder: GPUCommandEncoder,
      p: GPUComputePipeline,
      g: GPUBindGroup,
      timestampWrites?: GPUComputePassTimestampWrites,
    ) => {
      const pass = encoder.beginComputePass({ timestampWrites });
      pass.setPipeline(p);
      pass.setBindGroup(0, g);
      pass.dispatchWorkgroupsIndirect(args, 0);
      pass.end();
    };
    const cullPass = (
      encoder: GPUCommandEncoder,
      p: GPUComputePipeline,
      g: GPUBindGroup,
      timestampWrites?: GPUComputePassTimestampWrites,
    ) => {
      const pass = encoder.beginComputePass({ timestampWrites });
      pass.setPipeline(p);
      pass.setBindGroup(0, g);
      dispatch(pass, instanceGroups);
      pass.setPipeline(kernels.args);
      pass.setBindGroup(0, argsGroup);
      pass.dispatchWorkgroups(1);
      pass.end();
    };
    /**
     * Two-phase occlusion (exact bounds): the instances visible in the previous frame, the Hi-Z of their pixels, then
     * every instance tested against it (the newly visible drawn, the visible bits updated)
     */
    const twoPhase = (encoder: GPUCommandEncoder, ts: (GPUComputePassTimestampWrites | undefined)[] = []) => {
      clear(encoder, [0, 4]);
      cullPass(encoder, kernels.cullPhase1, cullPhase1Group, ts[0]);
      compactPass(encoder, kernels.rasterTight, tightGroup, ts[1]);
      const build = encoder.beginComputePass({ timestampWrites: ts[2] });
      build.setPipeline(kernels.hizBuild);
      hizLevels.forEach((d, l) => {
        build.setBindGroup(0, hizGroups[l]);
        build.dispatchWorkgroups(Math.ceil(d.x / 8), Math.ceil(d.y / 8));
      });
      build.end();
      clear(encoder, [4]);
      cullPass(encoder, kernels.cullPhase2, cullPhase2Group, ts[3]);
      compactPass(encoder, kernels.rasterTight, tightGroup, ts[4]);
    };
    if (BenchmarkEnv.rasterProfile === 'raster') {
      // a profiler's run: the one-phase raster and the resolve alone, frame after frame, without end
      for (let frame = 0; ; ++frame) {
        const encoder = device.createCommandEncoder();
        onePhase(encoder);
        const pass = encoder.beginRenderPass({
          colorAttachments: [{ view: color.createView(), loadOp: 'clear', storeOp: 'store', clearValue: [0, 0, 0, 0] }],
          depthStencilAttachment: {
            view: depth.createView(),
            depthLoadOp: 'clear',
            depthStoreOp: 'store',
            depthClearValue: 0,
          },
        });
        pass.setPipeline(resolvePipeline);
        pass.setBindGroup(0, resolveGroup);
        pass.draw(3);
        pass.end();
        device.queue.submit([encoder.finish()]);
        // frames queued in batches: the GPU busy throughout (a wait per frame left it idle 99 % of the time)
        if (frame % 20 === 19) await device.queue.onSubmittedWorkDone();
      }
    }
    const t = await timeFrames(
      device,
      15,
      [
        'clear',
        'raster64',
        'resolve',
        'raster64Early',
        'depth32',
        'id32',
        'cull',
        'raster64Compact',
        'rasterTight',
        'rasterTile',
        'cull1',
        'raster1',
        'hizBuild',
        'cull2',
        'raster2',
        'cullOne',
        'rasterOne',
        'hizOne',
      ],
      (encoder, ts) => {
        if (input.rotate) step(uniforms);
        clear(encoder, [0], ts(0));
        kernelPass(encoder, kernels.raster64, g64, ts(1));
        const pass = encoder.beginRenderPass({
          colorAttachments: [{ view: color.createView(), loadOp: 'clear', storeOp: 'store', clearValue: [0, 0, 0, 0] }],
          depthStencilAttachment: {
            view: depth.createView(),
            depthLoadOp: 'clear',
            depthStoreOp: 'store',
            depthClearValue: 0,
          },
          timestampWrites: ts(2),
        });
        pass.setPipeline(resolvePipeline);
        pass.setBindGroup(0, resolveGroup);
        pass.draw(3);
        pass.end();
        clear(encoder, [0, 1]);
        kernelPass(encoder, kernels.raster64Early, gEarly, ts(3));
        clear(encoder, [1, 2]);
        kernelPass(encoder, kernels.depth32, gDepth, ts(4));
        kernelPass(encoder, kernels.id32, gId, ts(5));
        // the compacted variant last: its buffer is the one the coverage reads
        clear(encoder, [0, 4]);
        cullPass(encoder, kernels.cull, cullGroup, ts(6));
        compactPass(encoder, kernels.raster64Compact, compactGroup, ts(7));
        // the same chunks: exact pixel bounds, then with the workgroup tile
        clear(encoder, [0]);
        compactPass(encoder, kernels.rasterTight, tightGroup, ts(8));
        clear(encoder, [0]);
        compactPass(encoder, kernels.rasterTile, tileGroup, ts(9));
        onePhase(encoder, [ts(15), ts(16), ts(17)]);
        twoPhase(encoder, [ts(10), ts(11), ts(12), ts(13), ts(14)]);
      },
    );
    // one-phase alone (its passes and the resolve, no other variant in the frame): against the times above
    const alone = await timeFrames(device, 15, ['cullOne', 'rasterOne', 'hizOne', 'resolve'], (encoder, ts) => {
      onePhase(encoder, [ts(0), ts(1), ts(2)]);
      const pass = encoder.beginRenderPass({
        colorAttachments: [{ view: color.createView(), loadOp: 'clear', storeOp: 'store', clearValue: [0, 0, 0, 0] }],
        depthStencilAttachment: {
          view: depth.createView(),
          depthLoadOp: 'clear',
          depthStoreOp: 'store',
          depthClearValue: 0,
        },
        timestampWrites: ts(3),
      });
      pass.setPipeline(resolvePipeline);
      pass.setBindGroup(0, resolveGroup);
      pass.draw(3);
      pass.end();
    });
    chunkCount = (await readBuffer(device, counter, 16))[0];
    const r = maxRadius > 1e8 ? 'all' : `R ${maxRadius}`;
    const occlusion = t.cull1 + t.raster1 + t.hizBuild + t.cull2 + t.raster2;
    notes.push(
      `${r}: clear ${f(t.clear)}, raster64 ${f(t.raster64)} + resolve ${f(t.resolve)} = ${f(t.raster64 + t.resolve)} ms; early-out ${f(t.raster64Early)}; 32-bit two passes ${f(t.depth32 + t.id32)}; compacted: cull ${f(t.cull)} + raster64 ${f(t.raster64Compact)} + resolve = ${f(t.cull + t.raster64Compact + t.resolve)} ms`,
      `${r}, one-phase alone: cull ${f(alone.cullOne)}, raster ${f(alone.rasterOne)}, build ${f(alone.hizOne)}, resolve ${f(alone.resolve)} = ${f(alone.cullOne + alone.rasterOne + alone.hizOne + alone.resolve)} ms`,
      `${r}, compacted with cull and resolve: exact bounds ${f(t.cull + t.rasterTight + t.resolve)}, with the tile ${f(t.cull + t.rasterTile + t.resolve)}, exact bounds with one-phase Hi-Z ${f(t.cullOne + t.rasterOne + t.hizOne + t.resolve)} (cull ${f(t.cullOne)}, raster ${f(t.rasterOne)}, build ${f(t.hizOne)}), with two-phase Hi-Z ${f(occlusion + t.resolve)} (phase 1: cull ${f(t.cull1)}, raster ${f(t.raster1)}; build ${f(t.hizBuild)}; phase 2: cull ${f(t.cull2)}, raster ${f(t.raster2)}) ms`,
    );
    if (maxRadius > 1e8) {
      // one more frame (a step on where rotating) without occlusion and with: the pixels whose sphere or depth differ
      if (input.rotate) step(uniforms);
      const frame = async (encode: (encoder: GPUCommandEncoder) => void) => {
        const encoder = device.createCommandEncoder();
        encode(encoder);
        device.queue.submit([encoder.finish()]);
        return readBuffer(device, vis, pixels * 8);
      };
      const reference = await frame((encoder) => {
        clear(encoder, [0, 4]);
        cullPass(encoder, kernels.cull, cullGroup);
        compactPass(encoder, kernels.rasterTight, tightGroup);
      });
      const differ = (culled: Uint32Array) => {
        let ids = 0;
        let depths = 0;
        for (let i = 0; i < pixels; ++i) {
          if (reference[2 * i] !== culled[2 * i]) ids += 1;
          if (reference[2 * i + 1] !== culled[2 * i + 1]) depths += 1;
        }
        return `${n(ids)} pixels with another sphere, ${n(depths)} with another depth`;
      };
      const one = differ(await frame((encoder) => onePhase(encoder)));
      const two = differ(await frame((encoder) => twoPhase(encoder)));
      notes.push(
        `against no occlusion${input.rotate ? ` (rotating, the camera turned ${Math.round((Vec3.angle(startDirection, direction()) * 180) / Math.PI)}° from this path's start)` : ''}: one phase ${one}; two-phase ${two}; phase 2's chunks of the last timed frame ${n(chunkCount)}`,
      );
    }
    uniforms.destroy();
  }

  // coverage at R unlimited (the compacted variant's buffer)
  const words = await readBuffer(device, vis, pixels * 8);
  const visible = new Set<number>();
  let covered = 0;
  for (let i = 0; i < pixels; ++i) {
    if (words[2 * i] === 0) continue;
    covered += 1;
    visible.add(words[2 * i]);
  }
  if (BenchmarkEnv.images) BenchmarkEnv.images[`${id}/resolve`] = await textureImage(device, color, width, height);

  const b = binCounts!;
  const passes = impostors.gpuPasses ?? {};
  const opaque = Object.keys(passes).filter((k) => k.endsWith('opaque'));
  const pass = (k: string) => `${k.split(' > ').pop()} ${f(passes[k])}`;
  notes.unshift(
    `${input.label}: ${n(count)} spheres (${spheres.targets.length} targets, ${n(spheres.transforms.length / 16)} instances) at ${width}x${height}, ${spheres.levels.length ? `${spheres.levels.length} levels by distance` : 'no levels'}`,
    `drawn in the frustum by radius in pixels (+1): ≤2 ${n(b[0])}, ≤4 ${n(b[1])}, ≤8 ${n(b[2])}, above ${n(b[3])}`,
    `covered ${f((100 * covered) / pixels)} % of the pixels, ${n(visible.size)} spheres visible; visibility buffer ${f((pixels * 8) / 1e6)} MB`,
    `impostors (GPU ms per frame): ${opaque.length ? opaque.map(pass).join(', ') : '(--gpu-timing)'}${
      passes['Canvas3D.render'] !== undefined ? `; render ${f(passes['Canvas3D.render'])}` : ''
    }`,
  );

  for (const buffer of [
    atomBuffer,
    transformBuffer,
    targetBuffer,
    targetUniformBuffer,
    vis,
    depth32,
    id32,
    bins,
    instanceTargetBuffer,
    counter,
    chunks,
    args,
    capacityUniform,
    hiz,
    hizOne,
    visibleBuffer,
    ...hizLevelUniforms,
  ]) {
    buffer.destroy();
  }
  for (const c of clears) c.count.destroy();
  color.destroy();
  depth.destroy();
  return { stages: [{ label: 'prototype', ms: 0 }], triangles: 0, notes };
}

function rasterCase(id: string, description: string, setup: () => Promise<RasterInput>): BenchmarkCase<RasterInput> {
  return { id, description, setup, paths: { webgpu: (graphics, input) => runPrototype(id, graphics, input) } };
}

const NoLod = { lodLevels: [], lodPixelSize: 0 };

/**
 * The prototype on the globule and PetWorld (`postsynstubby_2-no_bonds.bcif` in `.benchmark/`), whole and zoomed,
 * at 640x360, 1280x720 and 2560x1440 (on demand: `--filter 'sphere-raster/*' --gpu-timing --browser-arg
 * --enable-webgpu-developer-features`)
 */
export const SphereRasterCases: BenchmarkCase<RasterInput>[] = [
  { label: '360p', size: { width: 640, height: 360 } },
  { label: '720p', size: Size },
  { label: '1440p', size: LargeSize },
].flatMap(({ label, size }) => [
  rasterCase(`sphere-raster/globule-300k-${label}`, `300k atoms as spacefill, ${label}`, async () => {
    const s = await getGlobule();
    return {
      label: 'globule',
      size,
      spheres: () => structureSpheres(s),
      scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, NoLod)]),
    };
  }),
  ...[undefined, 0.3].map((zoom) =>
    rasterCase(
      `sphere-raster/postsyn${zoom ? '-zoom' : ''}-${label}`,
      `PetWorld particle targets as spacefill${zoom ? ', zoomed to 0.3 of the distance' : ''}, the \`quality\` levels by distance, the impostors culled on the GPU, ${label}`,
      async () => {
        const list = await getPostsyn();
        const levels = SpheresGeometry.LodLevelsPresets.quality;
        return {
          label: `PetWorld${zoom ? ' zoomed' : ''}`,
          size,
          zoom,
          spheres: () => particleSpheres(list, levels),
          scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, levels, undefined, 0)]),
        };
      },
    ),
  ),
]);

/** PetWorld zoomed at 720p: the impostors also culled by occlusion (Hi-Z), and both without levels (every sphere) */
/**
 * The particles' levels by size on screen (`lodPixelSize`, the chain of `Spheres.LodPixelPresets`) at `pixelSize`:
 * their strides, and as distances for a perspective camera of `fov` and `height` pixels (the impostors' ranges in Å a
 * pixel covers, `getLodLevelsValue`, at the nominal radius)
 */
function pixelLevels(pixelSize: number, fov: number, height: number): LodLevels {
  const chain = SpheresGeometry.LodPixelPresets.quality.lodLevels;
  const k = (2 * Math.tan(fov / 2)) / height;
  return chain.map((l, j) => ({
    ...l,
    minDistance: j === 0 ? 0 : getLodPixelBoundary(chain[j - 1], pixelSize) / k,
    maxDistance: j === chain.length - 1 ? 1e7 : getLodPixelBoundary(l, pixelSize) / k,
  }));
}

/**
 * PetWorld with the particles' levels by size on screen (their default, `lodPixelSize` 3), at the `ultra` preset's 1.5
 * and at 1, whole and zoomed, at 1280x720 and 2560x1440: the impostors with them, the raster with them as distances
 * (45° field of view)
 */
/**
 * 3j3q, the HIV-1 capsid (`3j3q.bcif` in `.benchmark/`, about 2.4 M atoms, 1356 chains), as spacefill with its
 * defaults (no levels of detail), whole and zoomed: for the integrated GPU, where PetWorld is too large (one structure,
 * one instance: the occlusion columns cull nothing)
 */
SphereRasterCases.push(
  ...[
    { label: '720p', size: Size },
    { label: '1440p', size: LargeSize },
  ].flatMap(({ label, size }) =>
    [undefined, 0.3].map((zoom) =>
      rasterCase(
        `sphere-raster/3j3q${zoom ? '-zoom' : ''}-${label}`,
        `3j3q as spacefill${zoom ? ', zoomed' : ''}, ${label}`,
        async () => {
          const s = await getBenchmarkEntry('3j3q.bcif');
          return {
            label: `3j3q${zoom ? ' zoomed' : ''}`,
            size,
            zoom,
            spheres: () => structureSpheres(s),
            scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]),
          };
        },
      ),
    ),
  ),
);

SphereRasterCases.push(
  ...[3, 1.5, 1].flatMap((pixelSize) =>
    [
      { label: '720p', size: Size },
      { label: '1440p', size: LargeSize },
    ].flatMap(({ label, size }) =>
      [undefined, 0.3].flatMap((zoom) =>
        // rotating (the popping of one phase) at the default and the `ultra` preset's
        (pixelSize === 1 ? [false] : [false, true]).map((rotate) =>
          rasterCase(
            `sphere-raster/postsyn${zoom ? '-zoom' : ''}-px${pixelSize}${rotate ? '-rotating' : ''}-${label}`,
            `PetWorld${zoom ? ' zoomed' : ''}, levels by size on screen (lodPixelSize ${pixelSize}), the impostors culled on the GPU, ${size.width}x${size.height}${rotate ? ', the camera rotating 2° a frame' : ''}`,
            async () => {
              const list = await getPostsyn();
              return {
                label: `PetWorld${zoom ? ' zoomed' : ''}, lodPixelSize ${pixelSize}${rotate ? ', rotating' : ''}`,
                size,
                zoom,
                rotate,
                spheres: () => particleSpheres(list, pixelLevels(pixelSize, Math.PI / 4, size.height)),
                scene: (g) =>
                  Promise.all([
                    postsynRepr(
                      g,
                      list,
                      true,
                      undefined,
                      SpheresGeometry.LodPixelPresets.quality.lodLevels,
                      undefined,
                      pixelSize,
                    ),
                  ]),
              };
            },
          ),
        ),
      ),
    ),
  ),
  rasterCase(
    'sphere-raster/postsyn-zoom-hiz-720p',
    'PetWorld zoomed, the `quality` levels, the impostors culled on the GPU with Hi-Z (the raster: frustum only)',
    async () => {
      const list = await getPostsyn();
      const levels = SpheresGeometry.LodLevelsPresets.quality;
      return {
        label: 'PetWorld zoomed (impostors with Hi-Z)',
        size: Size,
        zoom: 0.3,
        props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
        spheres: () => particleSpheres(list, levels),
        scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, levels, undefined, 0)]),
      };
    },
  ),
  ...[undefined, 0.3].map((zoom) =>
    rasterCase(
      `sphere-raster/postsyn${zoom ? '-zoom' : ''}-rotating-720p`,
      `PetWorld${zoom ? ' zoomed' : ''}, the \`quality\` levels, the camera rotating 2° a frame (both sides)`,
      async () => {
        const list = await getPostsyn();
        const levels = SpheresGeometry.LodLevelsPresets.quality;
        return {
          label: `PetWorld${zoom ? ' zoomed' : ''}, rotating`,
          size: Size,
          zoom,
          rotate: true,
          spheres: () => particleSpheres(list, levels),
          scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, levels, undefined, 0)]),
        };
      },
    ),
  ),
  rasterCase(
    'sphere-raster/postsyn-zoom-nolod-720p',
    'PetWorld zoomed without levels: every sphere, on both sides',
    async () => {
      const list = await getPostsyn();
      return {
        label: 'PetWorld zoomed, no levels',
        size: Size,
        zoom: 0.3,
        spheres: () => particleSpheres(list, []),
        scene: (g) => Promise.all([postsynRepr(g, list, false, undefined, [], undefined, 0)]),
      };
    },
  ),
);

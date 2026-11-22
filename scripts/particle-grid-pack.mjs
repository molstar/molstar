/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A particle pack's list repeated on an nx x ny x nz grid (side by side, a gap of a fraction of the extent between),
// written as a pack. Nested by default: the input list once, a template, and a particle a copy placing it (a pack of
// version 7, its size the input's). With --flat the copies' particles written out (the reference, read by any reader):
// tile-ordered with --tiles, coarse levels with --levels. Nested, its template's tiles an axis with --template-tiles
// (`WriteParticlePackOptions.templateTiles`, by default about a store page's particles a tile). A large scene for the particle tiles' scaling (flat, the
// yeast cell 2x2: 126 M particles, 2.8 GB; 3x3: 283 M, 6.1 GB). Runs on the built libraries (`npx tsc -b`); flat, a
// few GB of heap a hundred million particles (`node --max-old-space-size=40000`).
//
// usage: node scripts/particle-grid-pack.mjs <src.mspk> <out.mspk> [--nx 2] [--ny 2] [--nz 1] [--gap 0.05]
//                                            [--flat [--tiles 32] [--levels]] [--template-tiles 16]
import fs from 'node:fs';

/** a module of a package's built lib (`packages/<pkg>/lib/<path>.js`) */
const load = (pkg, path) => import(new URL(`../packages/${pkg}/lib/${path}.js`, import.meta.url).href);
const { createParticleListFromPack, readParticlePack, writeParticlePack } = await load(
  'model',
  'formats/particles/pack',
);
const { GridPositions } = await load('core', 'math/geometry/grid-positions');
const { Particle, ParticleValues } = await load('model', 'model/particles/particle-list');
const { CustomProperties } = await load('model', 'model/custom-property');

const args = process.argv.slice(2);
const opt = (name, fallback) => {
  const i = args.indexOf(`--${name}`);
  return i < 0 ? fallback : args[i + 1];
};
const [src, out] = args;
if (!src || !out) {
  console.error(
    'usage: node scripts/particle-grid-pack.mjs <src.mspk> <out.mspk> [--nx 2] [--ny 2] [--nz 1] [--gap 0.05] ' +
      '[--flat [--tiles 32] [--levels]] [--template-tiles 16]',
  );
  process.exit(1);
}
const nx = Number(opt('nx', 2));
const ny = Number(opt('ny', 2));
const nz = Number(opt('nz', 1));
const gap = Number(opt('gap', 0.05));
const tiles = opt('tiles', undefined) !== undefined ? Number(opt('tiles')) : undefined;
const levels = args.includes('--levels');
const flat = args.includes('--flat');
if (!flat && (tiles !== undefined || levels)) {
  console.error('--tiles and --levels with --flat only: a nested pack is not tile-ordered');
  process.exit(1);
}
const copies = nx * ny * nz;
const size = `${nx}x${ny}${nz > 1 ? `x${nz}` : ''}`;

const t0 = performance.now();
const log = (s) => console.log(`[${((performance.now() - t0) / 1000).toFixed(1)} s] ${s}`);

const list = createParticleListFromPack(readParticlePack(new Uint8Array(fs.readFileSync(src))));
const n = list.count;
log(`${src}: ${n} particles`);

// positions as float32 (the writer puts them on its grid again), their bounds
const xyz = new Float32Array(n * 3);
const p = [0, 0, 0];
const min = [Infinity, Infinity, Infinity];
const max = [-Infinity, -Infinity, -Infinity];
for (let i = 0; i < n; ++i) {
  GridPositions.get(list.coordinates, i, p);
  for (let k = 0; k < 3; ++k) {
    xyz[i * 3 + k] = p[k];
    if (p[k] < min[k]) min[k] = p[k];
    if (p[k] > max[k]) max[k] = p[k];
  }
}
const step = [0, 1, 2].map((k) => (max[k] - min[k]) * (1 + gap));
log(
  `bounds ${min.map((v) => v.toFixed(0))} to ${max.map((v) => v.toFixed(0))}, a copy every ${step.map((v) => v.toFixed(0))} Å`,
);

/** copy `c`'s offset: x fastest, then y, then z */
const offset = (c) => [(c % nx) * step[0], (Math.floor(c / nx) % ny) * step[1], Math.floor(c / (nx * ny)) * step[2]];

/** the bytes of `bytes` written to `out` in chunks: one write of over 2 GB fails */
function writeOut(bytes) {
  log(`written in memory: ${(bytes.byteLength / 1e6).toFixed(1)} MB`);
  const fd = fs.openSync(out, 'w');
  const Chunk = 1 << 28;
  for (let o = 0; o < bytes.byteLength; o += Chunk) {
    fs.writeSync(fd, bytes, o, Math.min(Chunk, bytes.byteLength - o));
  }
  fs.closeSync(fd);
  log(`${out}: ${(fs.statSync(out).size / 1e6).toFixed(1)} MB`);
}

if (!flat) {
  // a copy a particle of the top list placing the input whole (about the origin: its particles where the flat grid's
  // are), its radius the input's reach from the origin
  const { sphere } = Particle.getBoundary(list);
  const reach = Math.hypot(sphere.center[0], sphere.center[1], sphere.center[2]) + sphere.radius;
  const coordinates = new Float32Array(copies * 3);
  for (let c = 0; c < copies; ++c) coordinates.set(offset(c), c * 3);
  const top = {
    label: `${list.label ?? 'particles'} ${size}`,
    count: copies,
    coordinates,
    targets: ParticleValues.ofRuns(Int32Array.of(0, copies), Int32Array.of(0)),
    targetInfo: new Map([[0, { name: list.label ?? 'copy' }]]),
    targetMapping: new Map([[0, { kind: 'particles', particles: list, center: [0, 0, 0] }]]),
    radii: ParticleValues.ofRuns(Int32Array.of(0, copies), Float32Array.of(reach)),
    getParticleLabel: (i) => `Copy ${i + 1}`,
    sourceData: { kind: 'particle-grid', name: 'particle-grid', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
  log(`${copies} copies of ${n} particles (${n * copies} expanded) in a ${size} grid, nested; writing`);
  const templateTiles = opt('template-tiles', undefined);
  writeOut(writeParticlePack(top, templateTiles === undefined ? {} : { templateTiles: Number(templateTiles) }));
  process.exit(0);
}

const coordinates = new Float32Array(n * 3 * copies);
for (let c = 0; c < copies; ++c) {
  const [dx, dy, dz] = offset(c);
  const o = c * n * 3;
  for (let i = 0; i < n; ++i) {
    coordinates[o + i * 3] = xyz[i * 3] + dx;
    coordinates[o + i * 3 + 1] = xyz[i * 3 + 1] + dy;
    coordinates[o + i * 3 + 2] = xyz[i * 3 + 2] + dz;
  }
}

/** an array of a value or more a particle repeated `copies` times */
function repeat(a) {
  if (!a) return a;
  const r = new a.constructor(a.length * copies);
  for (let c = 0; c < copies; ++c) r.set(a, c * a.length);
  return r;
}

/** values as arrays or runs repeated: runs' offsets by the copy's first particle */
function repeatValues(v) {
  if (!v) return v;
  if (!ParticleValues.isRuns(v)) return repeat(v);
  const runs = v.values.length;
  const offsets = new Int32Array(runs * copies + 1);
  for (let c = 0; c < copies; ++c) for (let r = 0; r < runs; ++r) offsets[c * runs + r] = v.offsets[r] + c * n;
  offsets[runs * copies] = n * copies;
  return ParticleValues.ofRuns(offsets, repeat(v.values));
}

function repeatFibers(f) {
  if (!f) return f;
  const offsets = new Int32Array(f.count * copies + 1);
  const indices = new Int32Array(f.indices.length * copies);
  for (let c = 0; c < copies; ++c) {
    for (let k = 0; k < f.count; ++k) offsets[c * f.count + k] = f.offsets[k] + c * f.indices.length;
    for (let k = 0; k < f.indices.length; ++k) indices[c * f.indices.length + k] = f.indices[k] + c * n;
  }
  offsets[f.count * copies] = f.indices.length * copies;
  return { count: f.count * copies, offsets, indices };
}

/** attributes (columns) dropped: the yeast has none */
function repeatAttributes(attributes) {
  if (attributes?.size) log(`attributes dropped: ${[...attributes.keys()].join(', ')}`);
  return undefined;
}

const grid = {
  ...list,
  label: `${list.label ?? 'particles'} ${size}`,
  count: n * copies,
  cell: undefined,
  keys: undefined,
  coordinates,
  rotations: repeat(list.rotations),
  targets: repeatValues(list.targets),
  compartments: repeatValues(list.compartments),
  entities: repeatValues(list.entities),
  radii: repeatValues(list.radii),
  fibers: repeatFibers(list.fibers),
  attributes: repeatAttributes(list.attributes),
  getParticleLabel: (i) => list.getParticleLabel(i % n),
};
log(`${grid.count} particles in a ${size} grid, flat; writing (tiles ${tiles ?? 'none'}, levels ${levels})`);
writeOut(writeParticlePack(grid, { tiles, levels }));

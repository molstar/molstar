/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Headless bench of the text parsers on synthetic input (generated in memory, deterministic), the way the plugin reads
// a download or a file (`Data.Response`, its bytes):
// - `text`: the bytes decoded to a string (`utf8Read`, what `SO.Data.readText` does), then parsed: the path before the
//   text parsers read bytes, timed with its decode;
// - `bytes`: the bytes parsed (`SO.Data.readBytes`), the path since (`TextByteTokenizer`).
// Per format and input: the parse, then the parse and a read of every value of every column of the result (what a
// model or a shape built from it reads; the columns are lazy), each the median of the runs. On a commit whose parsers
// take strings only, run `--input text`: the before of a before/after (the file bundles on its own: copy it there).
//
// usage (through scripts/pipeline-bench.mjs, which bundles this file):
//   node scripts/pipeline-bench.mjs parsers [--only xyz,gro,...] [--input text,bytes] [--scale 1] [--runs 5]
//
// Results, macOS (arm64), Node 26, medians of 5 (parse, then parse and read, ms; `text` on the parent commit of the byte
// tokenizers, `bytes` after it; the checksums equal):
//   format       MB   text before      bytes after
//   xyz          43    336 / 370        143 / 208
//   gro          86    130 / 392         66 / 293
//   lammps-dump  75    215 / 392        126 / 273
//   lammps-data  47    138 / 253         77 / 170
//   psf          33     86 / 137         53 / 110
//   sdf          25    151 / 207        143 / 211
//   cube         46    223 / 226        146 / 148
//   dx           45    222 / 224        144 / 147
//   obj          72    369 / 379        220 / 228
//   ply          34    103 / 175         62 / 126
//   tbl          51    292 / 474        264 / 429
//   csv          33    104 / 185         64 / 244   (every column read as strings too: a unique
//                                                   value decodes, a `substring` did not)
//
// Formats (sizes at --scale 1): xyz, gro, lammps-dump (20 frames of 100k atoms each), lammps-data (1M atoms), psf (500k
// atoms), sdf (50k molecules), cube, dx (160³ values), obj, ply (500k vertices, 1M faces), tbl (500k particles), csv
// (1M rows).

import { Task } from '@molstar/core/task';
import { utf8Encode, utf8Read } from '@molstar/core/util/utf8';
import { parseCube } from '@molstar/io/reader/cube/parser';
import { parseCsv } from '@molstar/io/reader/csv/parser';
import { parseDx } from '@molstar/io/reader/dx/parser';
import { parseDynamoTbl } from '@molstar/io/reader/dynamo/tbl';
import { parseGRO } from '@molstar/io/reader/gro/parser';
import { parseLammpsData } from '@molstar/io/reader/lammps/data/parser';
import { parseLammpsTrajectory } from '@molstar/io/reader/lammps/traj/parser';
import { parseObj } from '@molstar/io/reader/obj/parser';
import { parsePly } from '@molstar/io/reader/ply/parser';
import { parsePsf } from '@molstar/io/reader/psf/parser';
import { parseSdf } from '@molstar/io/reader/sdf/parser';
import { parseXyz } from '@molstar/io/reader/xyz/parser';

const args = process.argv.slice(2);
const value = (name: string, fallback: string) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};
const only = value('--only', '').split(',').filter(Boolean);
const inputs = value('--input', 'text,bytes').split(',') as ('text' | 'bytes')[];
const scale = +value('--scale', '1');
const runs = +value('--runs', '5');

/** A deterministic pseudo-random number in [0, 1) */
let seed = 1;
const rand = () => ((seed = (Math.imul(seed, 48271) + 11) % 0x7fffffff) & 0xffffff) / 0x1000000;
const f = (v: number, d = 3) => v.toFixed(d);
const pad = (s: string | number, w: number, left = false) => (left ? `${s}`.padEnd(w) : `${s}`.padStart(w));
const n = (count: number) => Math.max(1, Math.round(count * scale));
const Elements = ['C', 'N', 'O', 'H', 'S'];

function xyz() {
  const atoms = n(100_000);
  const out: string[] = [];
  for (let fr = 0; fr < 20; fr++) {
    out.push(`${atoms}`, `frame ${fr}`);
    for (let i = 0; i < atoms; i++)
      out.push(`${Elements[i % 5]} ${f(100 * rand())} ${f(100 * rand())} ${f(100 * rand())}`);
  }
  return out.join('\n') + '\n';
}

function gro() {
  const atoms = n(100_000);
  const out: string[] = [];
  for (let fr = 0; fr < 20; fr++) {
    out.push(`synthetic t= ${f(fr, 5)}`, pad(atoms, 5));
    for (let i = 0; i < atoms; i++) {
      const res = (Math.floor(i / 10) % 99999) + 1;
      const atom = (i % 99999) + 1;
      out.push(
        `${pad(res, 5)}${pad('ALA', 5, true)}${pad(['N', 'CA', 'C', 'O', 'CB'][i % 5], 5)}${pad(atom, 5)}` +
          `${pad(f(10 * rand()), 8)}${pad(f(10 * rand()), 8)}${pad(f(10 * rand()), 8)}`,
      );
    }
    out.push('  10.00000  10.00000  10.00000');
  }
  return out.join('\n') + '\n';
}

function lammpsDump() {
  const atoms = n(100_000);
  const out: string[] = [];
  for (let fr = 0; fr < 20; fr++) {
    out.push('ITEM: TIMESTEP', `${fr * 100}`, 'ITEM: NUMBER OF ATOMS', `${atoms}`);
    out.push('ITEM: BOX BOUNDS pp pp pp', '0 100', '0 100', '0 100', 'ITEM: ATOMS id mol type x y z');
    for (let i = 0; i < atoms; i++) {
      out.push(
        `${i + 1} ${Math.floor(i / 10) + 1} ${(i % 3) + 1} ${f(100 * rand(), 5)} ${f(100 * rand(), 5)} ${f(100 * rand(), 5)}`,
      );
    }
  }
  return out.join('\n') + '\n';
}

function lammpsData() {
  const atoms = n(1_000_000);
  const out = [
    'LAMMPS data file',
    '',
    `${atoms} atoms`,
    '3 atom types',
    '',
    '0 100 xlo xhi',
    '0 100 ylo yhi',
    '0 100 zlo zhi',
  ];
  out.push('', 'Atoms # full', '');
  for (let i = 0; i < atoms; i++) {
    out.push(
      `${i + 1} ${Math.floor(i / 10) + 1} ${(i % 3) + 1} ${f(rand() - 0.5, 4)} ${f(100 * rand(), 5)} ${f(100 * rand(), 5)} ${f(100 * rand(), 5)}`,
    );
  }
  return out.join('\n') + '\n';
}

function psf() {
  const atoms = n(500_000);
  const out = ['PSF', '', '       1 !NTITLE', '* synthetic', '', `${pad(atoms, 8)} !NATOM`];
  for (let i = 0; i < atoms; i++) {
    out.push(
      `${pad(i + 1, 8)} PROT ${pad((Math.floor(i / 10) % 9999) + 1, 4, true)} ALA  ${pad(['N', 'CA', 'C', 'O', 'CB'][i % 5], 4, true)} ` +
        `${pad(22, 4)} ${pad(f(rand() - 0.5, 6), 10)}   ${pad('12.0110', 10)}           0`,
    );
  }
  out.push('', `${pad(0, 8)} !NBOND: bonds`, '');
  return out.join('\n') + '\n';
}

function sdf() {
  const molecules = n(50_000);
  const out: string[] = [];
  for (let m = 0; m < molecules; m++) {
    out.push(`mol${m}`, '  synthetic', '', '  5  4  0  0  0  0            999 V2000');
    for (let i = 0; i < 5; i++) {
      out.push(
        `${pad(f(10 * rand(), 4), 10)}${pad(f(10 * rand(), 4), 10)}${pad(f(10 * rand(), 4), 10)} ${pad(Elements[i], 3, true)} 0  0  0  0  0  0  0  0  0  0  0  0`,
      );
    }
    for (let i = 2; i <= 5; i++) out.push(`  1${pad(i, 3)}  1  0  0  0  0`);
    out.push('M  END', '> <ID>', `${m}`, '', '$$$$');
  }
  return out.join('\n') + '\n';
}

function cube() {
  const g = n(160);
  const out = ['synthetic cube', 'values', '    1    0.000000    0.000000    0.000000'];
  out.push(`${pad(g, 5)}    0.200000    0.000000    0.000000`, `${pad(g, 5)}    0.000000    0.200000    0.000000`);
  out.push(`${pad(g, 5)}    0.000000    0.000000    0.200000`, '    6    6.000000    1.000000    1.000000    1.000000');
  const line: string[] = [];
  for (let i = 0, il = g * g * g; i < il; i++) {
    line.push((rand() - 0.5).toExponential(5));
    if (line.length === 6 || i === il - 1) {
      out.push(' ' + line.join(' '));
      line.length = 0;
    }
  }
  return out.join('\n') + '\n';
}

function dx() {
  const g = n(160);
  const out = [
    `object 1 class gridpositions counts ${g} ${g} ${g}`,
    'origin 0 0 0',
    'delta 0.5 0 0',
    'delta 0 0.5 0',
    'delta 0 0 0.5',
  ];
  out.push(
    `object 2 class gridconnections counts ${g} ${g} ${g}`,
    `object 3 class array type double rank 0 items ${g * g * g} data follows`,
  );
  const line: string[] = [];
  for (let i = 0, il = g * g * g; i < il; i++) {
    line.push((rand() - 0.5).toExponential(5));
    if (line.length === 3 || i === il - 1) {
      out.push(line.join(' '));
      line.length = 0;
    }
  }
  out.push('attribute "dep" string "positions"', 'object "regular positions regular connections" class field');
  return out.join('\n') + '\n';
}

function obj() {
  const vertices = n(500_000);
  const out: string[] = [];
  for (let i = 0; i < vertices; i++) out.push(`v ${f(rand(), 6)} ${f(rand(), 6)} ${f(rand(), 6)}`);
  for (let i = 0; i < vertices; i++) out.push(`vn ${f(rand(), 6)} ${f(rand(), 6)} ${f(rand(), 6)}`);
  for (let i = 0; i < 2 * vertices; i++) {
    const a = 1 + Math.floor(rand() * vertices),
      b = 1 + Math.floor(rand() * vertices),
      c = 1 + Math.floor(rand() * vertices);
    out.push(`f ${a}//${a} ${b}//${b} ${c}//${c}`);
  }
  return out.join('\n') + '\n';
}

function ply() {
  const vertices = n(500_000);
  const out = [
    'ply',
    'format ascii 1.0',
    `element vertex ${vertices}`,
    'property float x',
    'property float y',
    'property float z',
  ];
  out.push(`element face ${2 * vertices}`, 'property list uchar int vertex_indices', 'end_header');
  for (let i = 0; i < vertices; i++) out.push(`${f(rand(), 6)} ${f(rand(), 6)} ${f(rand(), 6)}`);
  for (let i = 0; i < 2 * vertices; i++) {
    out.push(`3 ${Math.floor(rand() * vertices)} ${Math.floor(rand() * vertices)} ${Math.floor(rand() * vertices)}`);
  }
  return out.join('\n') + '\n';
}

function tbl() {
  const particles = n(500_000);
  const out: string[] = [];
  for (let i = 0; i < particles; i++) {
    const row = [i + 1, 1, 1, 0, 0, 0, f(360 * rand(), 2), f(180 * rand(), 2), f(360 * rand(), 2)];
    for (let c = 9; c < 19; c++) row.push(0);
    row.push(1 + (i % 7), 1, 1, 1, f(1000 * rand(), 2), f(1000 * rand(), 2), f(1000 * rand(), 2));
    for (let c = 26; c < 35; c++) row.push(0);
    row.push(6.5);
    out.push(row.join(' '));
  }
  return out.join('\n') + '\n';
}

function csv() {
  const rows = n(1_000_000);
  const out = ['id,name,x,y,z'];
  for (let i = 0; i < rows; i++)
    out.push(`${i},atom${i % 97},${f(100 * rand())},${f(100 * rand())},${f(100 * rand())}`);
  return out.join('\n') + '\n';
}

// the parsers take `StringLike | Uint8Array` since the byte tokenizers, a string before
type Parse = (data: any) => Task<any>;
const Formats: { name: string; make: () => string; parse: Parse }[] = [
  { name: 'xyz', make: xyz, parse: parseXyz },
  { name: 'gro', make: gro, parse: parseGRO },
  { name: 'lammps-dump', make: lammpsDump, parse: parseLammpsTrajectory },
  { name: 'lammps-data', make: lammpsData, parse: parseLammpsData },
  { name: 'psf', make: psf, parse: parsePsf },
  { name: 'sdf', make: sdf, parse: parseSdf },
  { name: 'cube', make: cube, parse: (d) => parseCube(d, 'bench') },
  { name: 'dx', make: dx, parse: (d) => (parseDx as any)(d, 'bench', { binary: false }) },
  { name: 'obj', make: obj, parse: parseObj },
  { name: 'ply', make: ply, parse: parsePly },
  { name: 'tbl', make: tbl, parse: parseDynamoTbl },
  { name: 'csv', make: csv, parse: parseCsv },
];

/** Reads every value of every column in the result (objects and arrays to a depth of 6), a checksum of the numbers */
function readAll(o: any, depth = 0): number {
  if (!o || typeof o !== 'object' || depth > 6) return 0;
  if (ArrayBuffer.isView(o)) {
    // eager values (grids, meshes): summed for the checksum
    const a = o as unknown as ArrayLike<number>;
    let sum = 0;
    for (let i = 0, il = a.length; i < il; i++) sum += +a[i];
    return sum;
  }
  if (typeof o.value === 'function' && typeof o.rowCount === 'number') {
    let sum = 0;
    for (let i = 0, il = o.rowCount; i < il; i++) {
      const v = o.value(i);
      sum += typeof v === 'number' ? v : typeof v === 'string' ? v.length : (v?.count ?? 0);
    }
    return sum;
  }
  if (typeof o.getColumn === 'function' && Array.isArray(o.columnNames)) {
    let sum = 0;
    for (const c of o.columnNames) {
      const col = o.getColumn(c);
      for (let i = 0; i < col.rowCount; i++) sum += col.float(i) + col.str(i).length;
    }
    return sum;
  }
  // PLY: elements and their properties behind getters
  if (typeof o.getElement === 'function') {
    return o.elementNames.reduce((s: number, e: string) => s + readAll(o.getElement(e), depth + 1), 0);
  }
  if (typeof o.getProperty === 'function') {
    return o.propertyNames.reduce((s: number, p: string) => s + readAll(o.getProperty(p), depth + 1), 0);
  }
  let sum = 0;
  for (const v of o instanceof Map ? o.values() : Object.values(o)) sum += readAll(v, depth + 1);
  return sum;
}

async function parseOnce(format: (typeof Formats)[number], bytes: Uint8Array, input: 'text' | 'bytes') {
  const data = input === 'text' ? utf8Read(bytes) : bytes;
  const result = await format.parse(data).run();
  if (result.isError) throw new Error(`${format.name}: ${result.message}`);
  return result.result;
}

const median = (xs: number[]) => xs.slice().sort((a, b) => a - b)[Math.floor(xs.length / 2)];

console.log(`format        MB     input   parse ms   parse+read ms   checksum`);
for (const format of Formats) {
  if (only.length && !only.includes(format.name)) continue;
  seed = 1;
  const bytes = utf8Encode(format.make());
  for (const input of inputs) {
    const parse: number[] = [];
    const total: number[] = [];
    let checksum = 0;
    for (let r = 0; r < runs + 1; r++) {
      (globalThis as any).gc?.();
      const t0 = performance.now();
      const result = await parseOnce(format, bytes, input);
      const t1 = performance.now();
      checksum = readAll(result);
      const t2 = performance.now();
      // the first run warms up
      if (r === 0) continue;
      parse.push(t1 - t0);
      total.push(t2 - t0);
    }
    console.log(
      `${format.name.padEnd(12)} ${pad(f(bytes.length / 2 ** 20, 0), 4)}   ${input.padEnd(5)} ${pad(f(median(parse), 0), 10)} ${pad(f(median(total), 0), 15)}   ${checksum.toPrecision(8)}`,
    );
  }
}

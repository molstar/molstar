/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the many-structures bench (many.ts): N separate files of small and normal-sized entries made from the
// local sources below, cycled, each moved to its own cell of a grid so that the scene is a plate of N proteins.
// Nothing is shared between two loaded structures in Mol*, so plain copies are as costly as N different entries.
//
// usage: node scripts/pipeline-bench.mjs many-fixtures [--out <dir>] [--count 500]

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import path from 'node:path';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { to_mmCIF } from '@molstar/model/model/structure/export/mmcif';
import { ROOT } from './common.js';

const PRISM = 'packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures';
/** 1.2k–11.5k atoms: protein–ligand, protein–DNA, multi-chain, an AlphaFold model; 5.3k atoms on average */
export const MANY_SOURCES = [
  'data/examples/1cbs_full.bcif',
  `${PRISM}/1tup.bcif`,
  `${PRISM}/2ocj.bcif`,
  `${PRISM}/2xwr.bcif`,
  `${PRISM}/AF-P04637-F1-model_v6.bcif`,
  'data/examples/7qpd.fw2.cif',
];

const args = process.argv.slice(2);
const arg = (n: string, d?: string) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 ? args[i + 1] : d;
};
const outDir = arg('out', path.join(ROOT, 'tests/build/pipeline/fixtures/many'))!;
const count = +arg('count', '500')!;

async function load(file: string) {
  const full = path.join(ROOT, file);
  const parsed = file.endsWith('.bcif')
    ? await CIF.parseBinary(new Uint8Array(readFileSync(full))).run()
    : await CIF.parse(readFileSync(full, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  return Structure.ofModel(model);
}

const sources = await Promise.all(MANY_SOURCES.map(load));
mkdirSync(outDir, { recursive: true });
const side = Math.ceil(Math.sqrt(count));
const spacing = 80;
let atoms = 0;
for (let i = 0; i < count; i++) {
  const src = sources[i % sources.length];
  const c = src.boundary.sphere.center;
  const offset = Vec3.create(spacing * (i % side) - c[0], spacing * Math.floor(i / side) - c[1], -c[2]);
  const moved = Structure.transform(src, Mat4.fromTranslation(Mat4(), offset));
  const name = `${String(i).padStart(3, '0')}`;
  writeFileSync(
    path.join(outDir, `${name}.bcif`),
    to_mmCIF(`MANY_${name}`, moved, true, { copyAllCategories: true }) as Uint8Array,
  );
  atoms += src.elementCount;
}
console.log(`${count} files, ${atoms} atoms -> ${outDir}`);
for (let i = 0; i < sources.length; i++)
  console.log(`  ${MANY_SOURCES[i]}: ${sources[i].elementCount} atoms, ${sources[i].units.length} units`);

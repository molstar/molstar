#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 *
 * Measures what a cached `StructureWrapper` retains (heap + array buffers, after a full GC) per byte of its source
 * file: the factor `structure-wrapper.ts` uses to turn source bytes into `approximateSize`. Run on the built server
 * (`npx tsc -b servers/model`) with `--expose-gc`:
 *
 *   node --expose-gc servers/model/scripts/measure-retained.mjs <file.bcif|file.cif> ...
 *
 * Files from `profile-memory.mjs --local <dir>` work. Each file is measured in a fresh process. Reports the wrapper
 * alone and after a `full` query, with the model's `Structure` (and the lookups the query built) kept as well.
 */

import { execFileSync } from 'node:child_process';
import * as fs from 'node:fs';
import * as path from 'node:path';
import { fileURLToPath } from 'node:url';

if (typeof globalThis.gc !== 'function') {
  console.error('Run with node --expose-gc.');
  process.exit(1);
}

const MB = 1024 * 1024;

async function fullGc() {
  for (let i = 0; i < 2; i++) {
    globalThis.gc();
    await new Promise((resolve) => setTimeout(resolve, 50));
  }
  globalThis.gc();
}

function used() {
  const m = process.memoryUsage();
  return m.heapUsed + m.arrayBuffers;
}

/** Measures one file in this (fresh) process and prints a JSON row. */
async function measure(file) {
  const { readStructureWrapper, resolveStructures } = await import('../lib/server/structure-wrapper.js');
  const { StructureQuery, StructureSelection, Queries } = await import('@molstar/model/model/structure');
  const log = console.log;
  console.log = () => {}; // the server logs every step

  await fullGc();
  const before = used();
  const wrapper = await readStructureWrapper(file, '_local_', path.resolve(file), undefined, undefined);
  await fullGc();
  const wrapperOnly = used() - before;

  // What a cached wrapper holds once it has served a query: the model's Structure and its lazily built lookups.
  const structures = await resolveStructures(wrapper);
  for (const s of structures) StructureSelection.unionStructure(StructureQuery.run(Queries.generators.all, s));
  await fullGc();
  const withStructure = used() - before;

  log(
    JSON.stringify({
      file: path.basename(file),
      sourceSize: fs.statSync(file).size,
      wrapperOnly,
      withStructure,
      approximateSize: wrapper.approximateSize,
      // reading them here keeps both referenced through the measurement
      models: wrapper.models.length,
      structures: structures.length,
    }),
  );
}

if (process.argv[2] === '--child') {
  await measure(process.argv[3]);
} else {
  const self = fileURLToPath(import.meta.url);
  const rows = [];
  for (const file of process.argv.slice(2)) {
    const out = execFileSync(process.execPath, ['--expose-gc', '--max-old-space-size=8192', self, '--child', file], {
      encoding: 'utf8',
      maxBuffer: 16 * MB,
    });
    const r = JSON.parse(out.trim().split('\n').pop());
    rows.push(r);
    const x = (v) => `${(v / MB).toFixed(1).padStart(7)} MB (x${(v / r.sourceSize).toFixed(2)})`;
    console.log(
      `${r.file.padEnd(14)} source ${(r.sourceSize / MB).toFixed(1).padStart(7)} MB  wrapper ${x(r.wrapperOnly)}  ` +
        `+structure ${x(r.withStructure)}  estimate ${x(r.approximateSize)}`,
    );
  }
  const sum = (k) => rows.reduce((n, r) => n + r[k], 0);
  const ratio = (k) => (sum(k) / sum('sourceSize')).toFixed(2);
  console.log(
    `\ntotal: source ${(sum('sourceSize') / MB).toFixed(1)} MB, wrapper x${ratio('wrapperOnly')}, ` +
      `+structure x${ratio('withStructure')}, estimate x${ratio('approximateSize')}`,
  );
}

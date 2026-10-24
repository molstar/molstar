/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Headless pipeline bench in Node (no browser, no GPU): bundles tests/src/pipeline into tests/build/pipeline (sources
// through the `molstar-src` condition, no package build needed) and runs it. Timings depend on the machine and its
// load: run on demand on an otherwise idle machine. See .v6/perf-audit/00-profile.md and 13-large-scale.md.
//
// usage: node scripts/pipeline-bench.mjs fixtures [--sizes 4x4x3,6x6x6] [--text]   # large fixtures, once
//        node scripts/pipeline-bench.mjs run --fixture <name> [bench options]       # stage table (bench.ts)
//        node scripts/pipeline-bench.mjs many-fixtures [--count 500]                # N small entries, once
//        node scripts/pipeline-bench.mjs many --count 300 [many options]            # N structures (many.ts)
//        node scripts/pipeline-bench.mjs playback [playback options]                # per-frame playback (playback.ts)
//        node scripts/pipeline-bench.mjs parsers [parsers options]                  # text parsers, text vs bytes (parsers.ts)
//        node scripts/pipeline-bench.mjs analyze <profile dir> [topN]               # hotspots of a --profile run
//                                                                                    (CPU, or allocations with --alloc)
//
// `run`, `many`, `playback` and `parsers` pass --heap <MB> as --max-old-space-size (default 14000) and everything else to bench.ts / many.ts / playback.ts.

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import { SourceMap } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import * as esbuild from 'esbuild';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outdir = path.join(root, 'tests/build/pipeline');
const [command, ...rest] = process.argv.slice(2).filter((a) => a !== '--');

async function bundle(name) {
  await esbuild.build({
    absWorkingDir: root,
    entryPoints: [path.join(root, `tests/src/pipeline/${name}.ts`)],
    tsconfig: path.join(root, 'tests/tsconfig.json'),
    outfile: path.join(outdir, `${name}.mjs`),
    bundle: true,
    platform: 'node',
    format: 'esm',
    target: 'node22',
    sourcemap: true,
    conditions: ['molstar-src', 'import', 'default'],
    // as the app builds (scripts/esbuild/app.mjs): no development checks (isProductionMode, mutative's draft check)
    define: { 'process.env.NODE_ENV': '"production"' },
    logLevel: 'warning',
    banner: { js: "import { createRequire as __cr } from 'node:module'; const require = __cr(import.meta.url);" },
  });
  return path.join(outdir, `${name}.mjs`);
}

function node(file, args, heapMB) {
  const flags = ['--expose-gc', `--max-old-space-size=${heapMB}`, '--enable-source-maps'];
  const r = spawnSync(process.execPath, [...flags, file, ...args], { stdio: 'inherit' });
  process.exit(r.status ?? 1);
}

/** Per stage profile: top self and total functions mapped back to the TS sources by the bundle's source map */
function analyze(dir, topN) {
  // the map of the bundle that wrote the profiles (copied by `run --profile`); a later run rebuilds tests/build
  const own = ['bench', 'many', 'playback'].map((b) => path.join(dir, `${b}.mjs.map`)).find((f) => fs.existsSync(f));
  const map = new SourceMap(JSON.parse(fs.readFileSync(own ?? path.join(outdir, 'bench.mjs.map'), 'utf8')));
  const cache = new Map();
  const loc = (cf) => {
    const k = `${cf.functionName}@${cf.url}:${cf.lineNumber}:${cf.columnNumber}`;
    if (cache.has(k)) return cache.get(k);
    let r;
    if (!cf.url) r = `(${cf.functionName || 'native'})`;
    else if (!/(bench|many|playback)\.mjs/.test(cf.url))
      r = `${cf.functionName || '(anon)'} ${cf.url.replace(/^.*node_modules\//, '')}:${cf.lineNumber + 1}`;
    else {
      const e = map.findEntry(cf.lineNumber, cf.columnNumber);
      const src = (e.originalSource || '?')
        .replace(/^file:\/\//, '')
        .replace(root + '/', '')
        .replace(/^(\.\.\/)+/, '');
      r = `${cf.functionName || '(anon)'} ${src}:${(e.originalLine ?? -1) + 1}`;
    }
    cache.set(k, r);
    return r;
  };
  // frames that only drive the stage (bench, task runner, representation update plumbing) are left out of "total"
  const plumbing =
    /^\((root|program|idle)\)|^\(anon\) |tests\/src\/pipeline|core\/src\/task|(createOrUpdate|update|create|runInContext|f|create2) packages\/graphics\/src\/repr\/(representation|structure\/(units|complex)-(representation|visual)|structure\/representation)/;
  const stages = JSON.parse(fs.readFileSync(path.join(dir, 'stages.json'), 'utf8'));
  const stageOf = (f, ext) => Object.entries(stages).find(([k]) => k.replace(/[:/]/g, '_') === f.replace(ext, ''))?.[1];
  // sampled allocations (--alloc): bytes allocated per function (self) and below it (total)
  for (const f of fs.readdirSync(dir).filter((f) => f.endsWith('.heapprofile'))) {
    const { head } = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8'));
    const self = new Map();
    const total = new Map();
    let all = 0;
    const walk = (n, stack) => {
      const key = loc(n.callFrame);
      self.set(key, (self.get(key) ?? 0) + n.selfSize);
      all += n.selfSize;
      const next = stack.includes(key) ? stack : [...stack, key];
      for (const k of next) total.set(k, (total.get(k) ?? 0) + n.selfSize);
      for (const c of n.children ?? []) walk(c, next);
    };
    walk(head, []);
    const fmt = (m, filter) =>
      [...m.entries()]
        .filter(([k]) => !filter || !plumbing.test(k))
        .sort((a, b) => b[1] - a[1])
        .slice(0, topN)
        .map(
          ([k, b]) =>
            `    ${(b / 2 ** 20).toFixed(0).padStart(8)} MB ${((100 * b) / all).toFixed(0).padStart(3)}%  ${k}`,
        )
        .join('\n');
    const st = stageOf(f, '.heapprofile');
    console.log(
      `## ${f.replace('.heapprofile', '')}  sampled=${(all / 2 ** 20).toFixed(0)}MB  wall=${st?.ms.toFixed(0)}ms`,
    );
    console.log(`  self:\n${fmt(self, false)}`);
    console.log(`  total:\n${fmt(total, true)}`);
  }
  for (const f of fs.readdirSync(dir).filter((f) => f.endsWith('.cpuprofile'))) {
    const p = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8'));
    const byId = new Map(p.nodes.map((n) => [n.id, n]));
    const parent = new Map();
    for (const n of p.nodes) for (const c of n.children ?? []) parent.set(c, n.id);
    const excluded = new Set();
    for (const n of p.nodes) {
      for (let cur = n.id; cur !== undefined; cur = parent.get(cur)) {
        const cf = byId.get(cur).callFrame;
        if (
          cf.url === 'node:inspector' ||
          cf.url === 'node:internal/util' ||
          cf.functionName === '(idle)' ||
          cf.functionName === '(program)'
        ) {
          excluded.add(n.id);
          break;
        }
      }
    }
    const selfUs = new Map();
    let totalUs = 0;
    for (let i = 0; i < p.samples.length; i++) {
      const dt = p.timeDeltas[i] ?? 0;
      if (excluded.has(p.samples[i])) continue;
      selfUs.set(p.samples[i], (selfUs.get(p.samples[i]) ?? 0) + dt);
      totalUs += dt;
    }
    const self = new Map();
    const total = new Map();
    for (const [id, us] of selfUs) {
      const key = loc(byId.get(id).callFrame);
      self.set(key, (self.get(key) ?? 0) + us);
      const seen = new Set();
      for (let cur = id; cur !== undefined; cur = parent.get(cur)) {
        const k = loc(byId.get(cur).callFrame);
        if (seen.has(k)) continue;
        seen.add(k);
        total.set(k, (total.get(k) ?? 0) + us);
      }
    }
    const fmt = (m, filter) =>
      [...m.entries()]
        .filter(([k]) => !filter || !plumbing.test(k))
        .sort((a, b) => b[1] - a[1])
        .slice(0, topN)
        .map(
          ([k, us]) =>
            `    ${(us / 1000).toFixed(1).padStart(8)} ms ${((100 * us) / totalUs).toFixed(0).padStart(3)}%  ${k}`,
        )
        .join('\n');
    const name = f.replace('.cpuprofile', '');
    const st = stageOf(f, '.cpuprofile');
    console.log(
      `## ${name}  sampled=${(totalUs / 1000).toFixed(0)}ms  wall=${st?.ms.toFixed(0)}ms alloc=${st?.allocMB.toFixed(0)}MB`,
    );
    console.log(`  self:\n${fmt(self, false)}`);
    console.log(`  total:\n${fmt(total, true)}`);
  }
}

if (command === 'fixtures') {
  node(await bundle('fixtures'), rest, 14000);
} else if (command === 'many-fixtures') {
  node(await bundle('many-fixtures'), rest, 14000);
} else if (command === 'run' || command === 'many' || command === 'playback' || command === 'parsers') {
  const name = command === 'run' ? 'bench' : command;
  const i = rest.indexOf('--heap');
  const heap = i >= 0 ? +rest.splice(i, 2)[1] : 14000;
  const file = await bundle(name);
  const p = rest.indexOf('--profile');
  if (p >= 0) {
    fs.mkdirSync(rest[p + 1], { recursive: true });
    fs.copyFileSync(`${file}.map`, path.join(rest[p + 1], `${name}.mjs.map`));
  }
  node(file, rest, heap);
} else if (command === 'analyze') {
  analyze(rest[0], +(rest[1] ?? 10));
} else {
  console.error(
    'usage: node scripts/pipeline-bench.mjs fixtures|run|many-fixtures|many|playback|parsers|analyze ... (see the header of this file)',
  );
  process.exit(1);
}

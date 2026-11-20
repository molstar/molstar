/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Bounds since the third slice's close (the author, 2026-11-20: "reduce render first, then move the size targets",
// "give the bounds a little wiggle room"): its measured sizes after the render reduction, with about 10 % room.
// The size targets (performance.md, Size): the entries of this folder bundled for the browser, minified, ESM, split
// (what loads before any dynamic import is measured: the entry chunk and its static imports),
// from the sources (the `molstar-src` condition), and the table of their sizes, minified and gzipped, against the
// targets, with each entry's largest contributors (by folder, from esbuild's metafile).
//
//   node packages/seq/_measure/sizes.mjs [--top 12] [--files] [--only spec-read.ts]
//
// `--files` lists each entry's largest input files too; `--only` measures the one entry (a check runs it so). An entry
// with `unlisted` modules is checked to hold none of them (bytes in any chunk, a lazy one too): what its composition does
// not list is not in its bundle (viewer-composition "Injected, opt in"); the script exits 1 when one is.
//
// The entries sit outside any package, so `@molstar/*` is resolved from the viewer's folder (a plugin below). An entry
// with `own` (an opt-in runtime) is measured alone: what its code imports from outside the folders `own` names is left
// out of its bundle (the viewer's already), so its row is what composing it adds. The opt-in runtimes (`runtimes/`)
// are in no viewer entry, which `unlisted` checks.

import { build } from 'esbuild';
import { gzipSync } from 'node:zlib';
import { dirname, join, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '../..');
const viewerDir = join(here, '../viewer');
const args = process.argv.slice(2);
const top = Number(args[args.indexOf('--top') + 1]) || 12;
const files = args.includes('--files');
const only = args.includes('--only') ? args[args.indexOf('--only') + 1] : undefined;

/**
 * What the minimal viewer does not list, by its modules' paths: the navigators, stories and transitions (the story
 * config's defaults are the viewer config's, kept), the history, the rows behavior, the records runtime, the built-in
 * tables. (The core's rows layer is still read by the tracks layout: not listed here.)
 */
const Unlisted = [
  /seq\/core\/src\/render\/navigator(-export)?\.ts$/,
  /seq\/viewer\/src\/(navigators?|mount-navigator|behaviors\/navigator)\.ts$/,
  /seq\/core\/src\/story\/(?!config\.ts$)/,
  /seq\/viewer\/src\/story\.ts$/,
  /ui\/src\/history\//,
  /seq\/viewer\/src\/ui-history\.ts$/,
  /seq\/viewer\/src\/behaviors\/rows/,
  /seq\/core\/src\/data\/records\.ts$/,
  /seq\/core\/src\/render\/built-ins\.ts$/,
  /seq\/core\/src\/context\/default\.ts$/,
  /seq\/core\/src\/render\/render-svg\.ts$/,
  /seq\/core\/src\/runtimes\//,
];

/** The opt-in runtimes: loaded only when a composition lists one. */
const Runtimes = [/seq\/core\/src\/runtimes\//];

/** Target and entry of each case, bytes (minified). */
const AllCases = [
  { name: 'reading a spec', entry: 'spec-read.ts', target: 195_000 },
  { name: 'headless rendering', entry: 'headless.ts', target: 780_000 },
  { name: 'the viewer, default', entry: 'viewer.ts', target: 1_250_000, unlisted: Runtimes },
  { name: 'viewer, minimal', entry: 'viewer-minimal.ts', target: 880_000, unlisted: Unlisted },
  { name: 'bigWig runtime', entry: 'bigwig.ts', target: 30_000, own: /seq\/core\/src\/runtimes\// },
];

const Cases = only ? AllCases.filter((c) => c.entry === only) : AllCases;
if (!Cases.length) throw new Error(`--only ${only}: no such entry (${AllCases.map((c) => c.entry).join(', ')})`);

const resolveWorkspace = {
  name: 'resolve-workspace',
  setup(b) {
    b.onResolve({ filter: /^@molstar\// }, async (a) => {
      if (a.pluginData?.skip) return;
      return b.resolve(a.path, { resolveDir: viewerDir, kind: a.kind, pluginData: { skip: true } });
    });
  },
};

/** Leaves out what an `own` entry's code imports from outside its folders: it is counted in the viewer's rows. */
const ownOnly = (own) => ({
  name: 'own-only',
  setup(b) {
    b.onResolve({ filter: /.*/ }, async (a) => {
      if (a.pluginData?.own || a.kind === 'entry-point' || !own.test(a.importer)) return;
      const r = await b.resolve(a.path, { resolveDir: a.resolveDir, kind: a.kind, pluginData: { own: true } });
      return r.errors.length || own.test(r.path) ? r : { path: a.path, external: true };
    });
  },
});

/** The chunks a split build loads before any dynamic import: the entry and what it imports statically, transitively. */
function loaded(meta) {
  const entry = Object.keys(meta.outputs).find((k) => meta.outputs[k].entryPoint !== undefined);
  const out = new Set([entry]);
  for (const k of out)
    for (const i of meta.outputs[k].imports) if (i.kind === 'import-statement' && !i.external) out.add(i.path);
  return [...out];
}

/** The inputs of the chunks loaded first, with their bytes in them. */
const loadedInputs = (meta) => loaded(meta).flatMap((k) => Object.entries(meta.outputs[k].inputs));

const kB = (n) => `${(n / 1000).toFixed(1)} kB`;

/** The bytes each input contributes to the output, grouped by its first `depth` path segments. */
function contributors(meta, depth) {
  const groups = new Map();
  for (const [file, v] of loadedInputs(meta)) {
    // `core/src/repr/tree/node.ts`, `seq/core/src/render/tracks.ts`: the package and its folders
    const parts = file.replace(/^packages\//, '').split('/');
    const key = parts.slice(0, depth).join('/');
    groups.set(key, (groups.get(key) ?? 0) + v.bytesInOutput);
  }
  return [...groups].sort((a, b) => b[1] - a[1]);
}

const rows = [];
for (const c of Cases) {
  const r = await build({
    entryPoints: [join(here, c.entry)],
    bundle: true,
    minify: true,
    format: 'esm',
    platform: 'browser',
    conditions: ['molstar-src'],
    write: false,
    metafile: true,
    logLevel: 'error',
    absWorkingDir: root,
    // split, as an application bundles: a dynamic import (MessagePack in the `seq` runtime) is a chunk of its own, not
    // paid by an entry that never loads it; what loads before any dynamic import is measured
    splitting: true,
    outdir: join(here, 'out'),
    plugins: c.own ? [ownOnly(c.own), resolveWorkspace] : [resolveWorkspace],
  });
  const chunks = loaded(r.metafile).map((k) => r.outputFiles.find((f) => f.path === join(root, k)).contents);
  const code = Buffer.concat(chunks);
  rows.push({ ...c, min: code.length, gzip: gzipSync(code, { level: 9 }).length, meta: r.metafile });
}

const pad = (s, n) => String(s).padEnd(n);
console.log(`${pad('entry', 22)}${pad('minified', 12)}${pad('gzip', 12)}${pad('target', 12)}result`);
for (const r of rows)
  console.log(
    `${pad(r.name, 22)}${pad(kB(r.min), 12)}${pad(kB(r.gzip), 12)}${pad(`< ${kB(r.target)}`, 12)}${r.min < r.target ? 'pass' : 'FAIL'}`,
  );

// each entry's unlisted modules: none may have a byte in any chunk (no lazy chunk either: the design has none)
let held = 0;
for (const r of rows) {
  if (!r.unlisted) continue;
  const bytes = new Map();
  for (const o of Object.values(r.meta.outputs))
    for (const [f, v] of Object.entries(o.inputs))
      if (v.bytesInOutput > 0) bytes.set(f, (bytes.get(f) ?? 0) + v.bytesInOutput);
  const found = [...bytes].filter(([f]) => r.unlisted.some((p) => p.test(f)));
  held += found.length;
  console.log(`\n${r.name}: unlisted modules ${found.length ? 'FOUND' : 'none'}`);
  for (const [f, n] of found) console.log(`  ${pad(kB(n), 10)}${f}`);
}

for (const r of rows) {
  console.log(`\n${r.name}: largest contributors (minified bytes in the output)`);
  for (const [k, v] of contributors(r.meta, 4).slice(0, top)) console.log(`  ${pad(kB(v), 10)}${k}`);
  if (files) {
    // the files of every chunk loaded first (a file can be in several chunks: its bytes add up)
    const bytes = new Map();
    for (const [f, v] of loadedInputs(r.meta)) bytes.set(f, (bytes.get(f) ?? 0) + v.bytesInOutput);
    console.log('  files:');
    for (const [f, n] of [...bytes].sort((a, b) => b[1] - a[1]).slice(0, top))
      console.log(`    ${pad(kB(n), 10)}${relative(root, join(root, f))}`);
  }
}

if (held) process.exitCode = 1;

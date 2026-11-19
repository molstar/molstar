/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The size targets (performance.md, Size): the three entries of this folder bundled for the browser, minified, ESM, split
// (what loads before any dynamic import is measured: the entry chunk and its static imports),
// from the sources (the `molstar-src` condition), and the table of their sizes, minified and gzipped, against the
// targets, with each entry's largest contributors (by folder, from esbuild's metafile).
//
//   node packages/seq/_measure/sizes.mjs [--top 12] [--files]
//
// The entries sit outside any package, so `@molstar/*` is resolved from the viewer's folder (a plugin below).

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

/** Target and entry of each case, bytes (minified). */
const Cases = [
  { name: 'reading a spec', entry: 'spec-read.ts', target: 160_000 },
  { name: 'headless rendering', entry: 'headless.ts', target: 300_000 },
  { name: 'the viewer, default', entry: 'viewer.ts', target: 600_000 },
];

const resolveWorkspace = {
  name: 'resolve-workspace',
  setup(b) {
    b.onResolve({ filter: /^@molstar\// }, async (a) => {
      if (a.pluginData?.skip) return;
      return b.resolve(a.path, { resolveDir: viewerDir, kind: a.kind, pluginData: { skip: true } });
    });
  },
};

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
    plugins: [resolveWorkspace],
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

for (const r of rows) {
  console.log(`\n${r.name}: largest contributors (minified bytes in the output)`);
  for (const [k, v] of contributors(r.meta, 4).slice(0, top)) console.log(`  ${pad(kB(v), 10)}${k}`);
  if (files) {
    const out = r.meta.outputs[entryOutput(r.meta)];
    console.log('  files:');
    for (const [f, v] of Object.entries(out.inputs)
      .sort((a, b) => b[1].bytesInOutput - a[1].bytesInOutput)
      .slice(0, top))
      console.log(`    ${pad(kB(v.bytesInOutput), 10)}${relative(root, join(root, f))}`);
  }
}

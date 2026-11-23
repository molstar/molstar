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
  // the structure translation is a host's library, in no viewer's composition
  /seq\/core\/src\/structure\//,
  // the track headers' behavior, its icons and its menu and parts (a host lists it: the default spec does)
  /seq\/viewer\/src\/(behaviors\/header[-a-z]*|header-[-a-z]+|mount-header-menu)\.ts$/,
];

/**
 * What loads only when used (a dynamic import): in no entry's first chunks, whatever it lists. The genome chunk and
 * its query forms, the first time a spec holds an assembly (or a query names a locus or a block).
 */
const Deferred = [
  /seq\/core\/src\/resolve\/(assemblies|assembly-names|genome-model|genome-chunk|genome-resolve|genome-issues)\.ts$/,
  /seq\/core\/src\/query\/genome-forms\.ts$/,
  /seq\/core\/src\/lang\/include-pattern\.ts$/,
  // the bridge chunk: a view with a bridge loads it (render's `bridge-load.ts`)
  /seq\/core\/src\/render\/(bridges|pairs-ribbons(-lay|-read|-draw)?|layout\/bridge)\.ts$/,
];

/** The bridge chunk's row (minified; the design's estimate was about 6 kB, measured 14.6 at LA, 20.8 with the per-side reads, merged blocks and the hatch; 24 until LA2's per-side reads over assemblies, raised to 28 by the author, 2026-11-22: "Raise to 28", measured 25.2): the default viewer's chunk holding `render/bridges.ts`. */
const BridgeChunk = { name: 'bridge chunk (lazy)', file: /seq\/core\/src\/render\/bridges\.ts$/, target: 28_000 };

/** The opt-in runtimes: loaded only when a composition lists one. */
const Runtimes = [/seq\/core\/src\/runtimes\//, /seq\/core\/src\/structure\//];

/** Target and entry of each case, bytes (minified). */
const AllCases = [
  { name: 'reading a spec', entry: 'spec-read.ts', target: 265_000 },
  { name: 'headless rendering', entry: 'headless.ts', target: 980_000 },
  { name: 'the viewer, default', entry: 'viewer.ts', target: 1_545_000, unlisted: Runtimes },
  { name: 'viewer, minimal', entry: 'viewer-minimal.ts', target: 1_090_000, unlisted: Unlisted },
  { name: 'bigWig runtime', entry: 'bigwig.ts', target: 30_000, own: /seq\/core\/src\/runtimes\// },
  // the genome ring extension (ring-extension.md, The package): its own code only, on the viewer's extension entry
  { name: 'genome ring (rondo)', entry: 'rondo.ts', target: 32_000, own: /seq\/ext-rondo\/src\// },
  // the BED family's runtime with io's BGZF (no viewer bundles it); its reads by region (io's tabix, the region lines)
  // are its lazy chunk, loaded when the first indexed source opens: listed apart (`lazy`)
  {
    name: 'BED runtime',
    entry: 'bed.ts',
    target: 30_000,
    own: /seq\/core\/src\/runtimes\/|(^|\/)packages\/io\/src\/common\//,
    lazy: true,
  },
  // GFF3 and GTF read whole (no viewer bundles it): its own code; io's BGZF (a bgzipped file read whole) a lazy chunk on a gzip file's magic; its reads
  // by region (io's tabix, the region lines) are its lazy chunk, loaded when the first indexed source opens: listed
  // apart (`lazy`)
  {
    name: 'GFF runtime',
    entry: 'gff.ts',
    target: 30_000,
    own: /seq\/core\/src\/runtimes\/|(^|\/)packages\/io\/src\/common\//,
    lazy: true,
  },
  // GenBank read whole (no viewer bundles it): its own code; io's BGZF (a bgzipped file read whole) a lazy chunk on a gzip file's magic
  {
    name: 'GenBank runtime',
    entry: 'genbank.ts',
    target: 30_000,
    own: /seq\/core\/src\/runtimes\/|(^|\/)packages\/io\/src\/common\//,

    lazy: true,
  },
  // VCF sites (no viewer bundles it): its own code; io's BGZF (a bgzipped file read whole) a lazy chunk on a gzip file's magic; its reads by region (io's
  // tabix, the region lines) are its lazy chunk, loaded when the first indexed source opens: listed apart (`lazy`)
  {
    name: 'VCF runtime',
    entry: 'vcf.ts',
    target: 30_000,
    own: /seq\/core\/src\/runtimes\/|(^|\/)packages\/io\/src\/common\//,
    lazy: true,
  },
  // the mmCIF runtime and Seq io's converter (in no viewer entry); `@molstar/io`'s parser is its lazy chunk, loaded
  // when the first mmCIF source opens: listed apart (`lazy`), measured with every import it brings that the default
  // viewer's bundle does not hold (its `@molstar/core` imports too): what a host pays on first open
  {
    name: 'mmCIF runtime',
    entry: 'mmcif.ts',
    target: 30_000,
    own: /seq\/core\/src\/(runtimes|io)\/|(^|\/)packages\/io\/src\//,
    lazy: true,
  },
];

const Cases = only ? AllCases.filter((c) => c.entry === only) : AllCases;
if (!Cases.length) throw new Error(`--only ${only}: no such entry (${AllCases.map((c) => c.entry).join(', ')})`);

const resolveWorkspace = {
  name: 'resolve-workspace',
  setup(b) {
    b.onResolve({ filter: /^@molstar\// }, async (a) => {
      if (a.pluginData?.skip) return;
      // a package the viewer does not depend on (`@molstar/io`, an opt-in runtime's) from its importer's folder
      const r = await b.resolve(a.path, { resolveDir: viewerDir, kind: a.kind, pluginData: { skip: true } });
      return r.errors.length
        ? b.resolve(a.path, { resolveDir: a.resolveDir, kind: a.kind, pluginData: { skip: true } })
        : r;
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
  // the measured entry (this folder's), never a dynamic import's chunk, which esbuild marks as an entry point too
  const entry = Object.keys(meta.outputs).find((k) => meta.outputs[k].entryPoint?.includes('_measure/'));
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

/** A split build of an entry, `own` leaving out what its folders import (none: every import bundled). */
const bundle = (entry, own) =>
  build({
    entryPoints: [join(here, entry)],
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
    plugins: own ? [ownOnly(own), resolveWorkspace] : [resolveWorkspace],
  });

/** The inputs of an entry's default viewer bundle (any chunk): what a host composing the viewer already pays. */
let viewerInputs;
const inViewer = async () => {
  if (!viewerInputs) {
    const v = await bundle('viewer.ts');
    viewerInputs = new Set(Object.values(v.metafile.outputs).flatMap((o) => Object.keys(o.inputs)));
  }
  return viewerInputs;
};

/** An entry's lazy chunks with every import bundled: the bytes past the viewer's, and the largest folders of them. */
async function lazyCost(entry) {
  const r = await bundle(entry);
  const first = new Set(loaded(r.metafile));
  const viewer = await inViewer();
  const groups = new Map();
  let bytes = 0;
  for (const [k, o] of Object.entries(r.metafile.outputs)) {
    if (first.has(k) || !k.endsWith('.js')) continue;
    for (const [f, v] of Object.entries(o.inputs)) {
      if (viewer.has(f)) continue;
      bytes += v.bytesInOutput;
      const key = f
        .replace(/^packages\//, '')
        .split('/')
        .slice(0, 4)
        .join('/');
      groups.set(key, (groups.get(key) ?? 0) + v.bytesInOutput);
    }
  }
  return { bytes, top: [...groups].sort((a, b) => b[1] - a[1]).slice(0, 5) };
}

const rows = [];
for (const c of Cases) {
  const r = await bundle(c.entry, c.own);
  const chunks = loaded(r.metafile).map((k) => r.outputFiles.find((f) => f.path === join(root, k)).contents);
  const code = Buffer.concat(chunks);
  // the chunks loaded after a dynamic import, with their imports: what an opt-in runtime loads on first use
  const lazy = c.lazy ? await lazyCost(c.entry) : undefined;
  rows.push({ ...c, min: code.length, gzip: gzipSync(code, { level: 9 }).length, lazyCost: lazy, meta: r.metafile });
}

const pad = (s, n) => String(s).padEnd(n);
console.log(`${pad('entry', 22)}${pad('minified', 12)}${pad('gzip', 12)}${pad('target', 12)}result`);
for (const r of rows)
  console.log(
    `${pad(r.name, 22)}${pad(kB(r.min), 12)}${pad(kB(r.gzip), 12)}${pad(`< ${kB(r.target)}`, 12)}${r.min < r.target ? 'pass' : 'FAIL'}`,
  );
// the bridge chunk: the default viewer's lazy chunk holding the bridges' entry (and what it alone imports)
const viewerRow = rows.find((r) => r.entry === 'viewer.ts');
if (viewerRow) {
  const first = new Set(loaded(viewerRow.meta));
  const chunk = Object.entries(viewerRow.meta.outputs).find(
    ([k, o]) => !first.has(k) && Object.keys(o.inputs).some((f) => BridgeChunk.file.test(f)),
  );
  const bytes = chunk ? chunk[1].bytes : 0;
  console.log(
    `${pad(BridgeChunk.name, 22)}${pad(kB(bytes), 12)}${pad('', 12)}${pad(`< ${kB(BridgeChunk.target)}`, 12)}${chunk && bytes < BridgeChunk.target ? 'pass' : 'FAIL'}`,
  );
  if (!chunk || bytes >= BridgeChunk.target) process.exitCode = 1;
}
for (const r of rows.filter((r) => r.lazyCost)) {
  console.log(
    `${r.name}: ${kB(r.lazyCost.bytes)} more on first open (its lazy chunks and their imports past the viewer's, minified):`,
  );
  for (const [k, v] of r.lazyCost.top) console.log(`  ${pad(kB(v), 10)}${k}`);
}

// each entry's lazy chunks (viewer-composition, Injected, opt in: a listed part used rarely or late may load on first
// use): each chunk's size and its largest inputs
for (const r of rows) {
  if (r.own) continue;
  const first = new Set(loaded(r.meta));
  const lazy = Object.entries(r.meta.outputs).filter(([k]) => !first.has(k) && k.endsWith('.js'));
  if (!lazy.length) continue;
  console.log(`\n${r.name}: lazy chunks (minified, loaded on first use)`);
  for (const [k, o] of lazy.sort((a, b) => b[1].bytes - a[1].bytes)) {
    const inputs = Object.entries(o.inputs)
      .sort((a, b) => b[1].bytesInOutput - a[1].bytesInOutput)
      .slice(0, 3)
      .map(([f, v]) => `${f.replace(/^packages\//, '')} ${kB(v.bytesInOutput)}`);
    console.log(`  ${pad(kB(o.bytes), 10)}${k.split('/').pop()}: ${inputs.join(', ')}`);
  }
}

// each entry's unlisted modules: none may have a byte in any chunk (a lazy one too)
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

// each entry's deferred modules: none may be in the chunks it loads first
for (const r of rows) {
  if (r.own) continue;
  const found = [...new Set(loadedInputs(r.meta).map(([f]) => f))].filter((f) => Deferred.some((p) => p.test(f)));
  held += found.length;
  if (found.length) console.log(`\n${r.name}: deferred modules loaded first: ${found.join(', ')}`);
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

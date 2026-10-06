/**
 * Bundle check of the slim-plugin acceptance target (.v6/designs/plugin-composition.md §12, plan step 4).
 *
 * Builds the example entry of `slim-exclusions.json` twice with esbuild, the way `scripts/esbuild/app.mjs` builds
 * examples (same `molstar-src` condition, production minification, source `.js` extension mapping):
 *   - split ESM (`splitting: true, format: 'esm'`), and
 *   - a single IIFE file,
 * and fails when the metafile of either build lists an excluded module. An excluded module is reported with the import
 * chain from the entry, taken from the metafile `imports`. The sizes of both builds and of the default-spec Viewer build
 * are printed in bytes for comparison (`--no-viewer` skips the Viewer build).
 *
 * Nothing is written to the repository: the builds run with `write: false`. Style and static asset imports resolve to
 * empty modules, so sizes are JavaScript bytes only.
 */
import * as esbuild from 'esbuild';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { exclusionErrors, findExcluded, loadExclusions } from './slim-exclusions.mjs';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');

const AssetLoaders = Object.fromEntries(
  [
    '.html',
    '.htm',
    '.ico',
    '.sdf',
    '.scss',
    '.css',
    '.png',
    '.jpg',
    '.jpeg',
    '.gif',
    '.svg',
    '.wasm',
    '.bin',
    '.dat',
  ].map((ext) => [ext, 'empty']),
);

/** Source files import `./x.js` for `./x.ts`; same as in `scripts/esbuild/app.mjs`. */
const sourceJsExtensionPlugin = {
  name: 'molstar-source-js-extension',
  setup(build) {
    build.onResolve({ filter: /^\./ }, (args) => {
      if (!args.path.endsWith('.js')) return null;
      const file = path.resolve(args.resolveDir, args.path);
      if (fs.existsSync(file)) return null;
      for (const ext of ['.ts', '.tsx']) {
        const source = file.slice(0, -'.js'.length) + ext;
        if (fs.existsSync(source)) return { path: source };
      }
      return null;
    });
  },
};

const BuildModes = {
  'split ESM': { format: 'esm', splitting: true },
  'single file': { format: 'iife', globalName: 'molstar' },
};

/** Builds `entry` (repository-relative) in memory and returns the esbuild result. */
export async function buildEntry({ root = repoRoot, entry, mode }) {
  const version = JSON.parse(fs.readFileSync(path.join(root, 'version.json'), 'utf8')).version;
  const entryDir = path.dirname(path.join(root, entry));
  let packageDir = entryDir;
  while (!fs.existsSync(path.join(packageDir, 'tsconfig.json')) && packageDir !== root)
    packageDir = path.dirname(packageDir);
  return esbuild.build({
    absWorkingDir: root,
    entryPoints: [entry],
    bundle: true,
    platform: 'browser',
    ...BuildModes[mode],
    write: false,
    outdir: path.join(os.tmpdir(), 'molstar-slim-bundle'),
    metafile: true,
    tsconfig: path.join(packageDir, 'tsconfig.json'),
    minify: true,
    minifyIdentifiers: false,
    conditions: ['molstar-src', 'import', 'default'],
    external: ['crypto', 'fs', 'path', 'stream'],
    loader: AssetLoaders,
    plugins: [sourceJsExtensionPlugin],
    define: {
      'process.env.NODE_ENV': JSON.stringify('production'),
      'process.env.DEBUG': JSON.stringify(false),
      __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(version),
    },
    logLevel: 'error',
  });
}

/** Total bytes of the JavaScript outputs of a metafile. */
export function outputBytes(metafile) {
  return Object.entries(metafile.outputs)
    .filter(([file]) => /\.m?js$/.test(file))
    .reduce((sum, [, output]) => sum + output.bytes, 0);
}

/** The shortest import chain from the entry points of the metafile to `target` (metafile input paths). */
export function importChain(metafile, target) {
  // With splitting, a dynamically imported module is an entry point of its own; the real entries are imported by nobody.
  const imported = new Set(Object.values(metafile.inputs).flatMap((input) => input.imports.map((i) => i.path)));
  const entries = Object.values(metafile.outputs)
    .map((output) => output.entryPoint)
    .filter((entryPoint) => entryPoint && metafile.inputs[entryPoint] && !imported.has(entryPoint));
  const previous = new Map(entries.map((entry) => [entry, undefined]));
  const queue = [...entries];
  for (let i = 0; i < queue.length; i++) {
    const current = queue[i];
    if (current === target) break;
    for (const imported of metafile.inputs[current]?.imports ?? []) {
      if (imported.external || previous.has(imported.path) || !metafile.inputs[imported.path]) continue;
      previous.set(imported.path, current);
      queue.push(imported.path);
    }
  }
  if (!previous.has(target)) return [target];
  const chain = [];
  for (let node = target; node !== undefined; node = previous.get(node)) chain.unshift(node);
  return chain;
}

/** Checks a metafile against the exclusions. Returns `{ errors, known }`; `known` are the known-leak modules it contains. */
export function checkMetafile(metafile, exclusions, label = 'bundle') {
  const files = Object.keys(metafile.inputs).map((file) => file.split(path.sep).join('/'));
  const { violations, known } = findExcluded(files, exclusions);
  const errors = [...violations].map(
    ([file, group]) =>
      `${label}: excluded module ${file} (${group}) is in the bundle, imported through: ${importChain(metafile, file).join(' -> ')}`,
  );
  return { errors, known };
}

const format = (n) => `${n.toLocaleString('en-US')} bytes`;

async function main() {
  const started = Date.now();
  const exclusions = loadExclusions(repoRoot);
  const errors = [...exclusionErrors(exclusions)];
  const sizes = [];
  for (const mode of Object.keys(BuildModes)) {
    const result = await buildEntry({ entry: exclusions.entry, mode });
    const bytes = outputBytes(result.metafile);
    sizes.push([`slim example, ${mode}`, bytes]);
    const checked = checkMetafile(result.metafile, exclusions, `slim example (${mode})`);
    errors.push(...checked.errors);
    for (const [file, leak] of checked.known) console.warn(`known leak (${mode}): ${file}: ${leak.reason}`);
    console.log(`slim example, ${mode}: ${Object.keys(result.metafile.inputs).length} modules, ${format(bytes)}`);
  }
  if (!process.argv.includes('--no-viewer')) {
    const result = await buildEntry({ entry: 'apps/viewer/src/index.ts', mode: 'single file' });
    const bytes = outputBytes(result.metafile);
    sizes.push(['default Viewer, single file', bytes]);
    console.log(`default Viewer, single file: ${Object.keys(result.metafile.inputs).length} modules, ${format(bytes)}`);
    const slim = sizes[1][1];
    console.log(`slim example is ${((slim / bytes) * 100).toFixed(1)}% of the default Viewer (single file).`);
  }
  if (errors.length) {
    console.error(errors.join('\n'));
    process.exitCode = 1;
  } else console.log(`Slim bundles contain no excluded module (${Date.now() - started} ms).`);
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) await main();

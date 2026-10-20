/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tests are not published: a package whose exports map reaches its sources through a wildcard (`./*` → `./src/*.ts`,
// `./lib/*.js`) needs `"./<dir>/_test/*": null` for every `_test/` folder, else the wildcard exports the tests and
// `check:workspace` fails on the build's missing `lib/**/_test/*.js` (it sees only built files, so only in CI). This is
// the same rule from the sources alone, without a build: fast enough for `pnpm lint`. `--fix` adds the entries.
//
// Usage: node scripts/workspace/test-exports.mjs [--fix]

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { resolveExport } from './exports.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const fix = process.argv.includes('--fix');
const inventory = JSON.parse(fs.readFileSync(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));

/** The `_test` folders under `src/`, relative to it, each with one of its sources. */
function testFolders(src) {
  const found = [];
  const walk = (dir) => {
    for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
      if (!ent.isDirectory() || ent.name === 'node_modules') continue;
      const full = path.join(dir, ent.name);
      if (ent.name === '_test') {
        const file = firstSource(full);
        if (file) found.push({ dir: path.relative(src, full).split(path.sep).join('/'), file });
        continue;
      }
      walk(full);
    }
  };
  walk(src);
  return found;
}

/** A `.ts`/`.tsx` file in the folder or below, relative to the folder, without its extension. */
function firstSource(dir) {
  for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
    if (ent.isFile() && /\.tsx?$/.test(ent.name) && !ent.name.endsWith('.d.ts')) return ent.name.replace(/\.tsx?$/, '');
  }
  for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
    if (!ent.isDirectory()) continue;
    const below = firstSource(path.join(dir, ent.name));
    if (below) return `${ent.name}/${below}`;
  }
  return undefined;
}

const problems = [];
let fixed = 0;
// The packages `check:workspace` checks the exports of: the published ones (a private package's map is not validated).
for (const pkg of inventory.packages ?? []) {
  if (pkg.private || pkg.kind === 'distribution') continue;
  const dir = path.join(root, pkg.path);
  const manifestPath = path.join(dir, 'package.json');
  const src = path.join(dir, 'src');
  if (!fs.existsSync(manifestPath) || !fs.existsSync(src)) continue;
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  const exports = manifest.exports;
  if (!exports || typeof exports !== 'object' || Array.isArray(exports)) continue;
  const missing = [];
  for (const { dir: testDir, file } of testFolders(src)) {
    if (resolveExport(exports, `./${testDir}/${file}`) == null) continue;
    missing.push(`./${testDir}/*`);
  }
  if (!missing.length) continue;
  if (!fix) {
    for (const key of missing) problems.push(`${pkg.name} (${pkg.path}/package.json): add "${key}": null to "exports"`);
    continue;
  }
  // The null entries kept together, sorted, after the entries they mask (Node picks the longest pattern either way).
  const entries = Object.entries(exports);
  const nulls = new Map(entries.filter(([, v]) => v === null));
  for (const key of missing) nulls.set(key, null);
  const others = entries.filter(([, v]) => v !== null);
  manifest.exports = Object.fromEntries([...others, ...[...nulls.entries()].sort(([a], [b]) => a.localeCompare(b))]);
  fs.writeFileSync(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`);
  fixed += missing.length;
  for (const key of missing) console.log(`${pkg.name}: added "${key}": null`);
}

if (problems.length) {
  console.error(
    `Test folders published through a wildcard export (each needs a null entry; \`node scripts/workspace/test-exports.mjs --fix\` adds them):\n${problems.map((p) => `  ${p}`).join('\n')}`,
  );
  process.exit(1);
}
if (fix && fixed) console.log(`Added ${fixed} null export${fixed === 1 ? '' : 's'}.`);

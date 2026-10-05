/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <malym@ibt.cas.cz>
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = new Set(process.argv.slice(2));
const known = new Set(['--build', '--lib', '--all']);
if (args.has('--help') || args.has('-h')) {
  console.log('Usage: node scripts/clean.js [--build] [--lib] [--all]');
  process.exit(0);
}
known.add('-h');
known.add('--help');
const unknown = [...args].filter((arg) => !known.has(arg));
if (unknown.length) throw new Error(`Unknown clean option(s): ${unknown.join(', ')}`);

const inventoryPath = path.join(root, 'scripts/workspace/inventory.json');
const inventory = JSON.parse(fs.readFileSync(inventoryPath, 'utf8'));
const cleanBuild = args.has('--build') || args.has('--all');
const cleanLib = args.has('--lib') || args.has('--all');
const targets = new Set();
if (cleanBuild) {
  targets.add(path.join(root, 'build'));
  targets.add(path.join(root, 'deploy/data'));
  for (const pkg of inventory.packages ?? []) targets.add(path.join(root, pkg.path, 'build'));
}
if (cleanLib) {
  targets.add(path.join(root, 'lib'));
  targets.add(path.join(root, 'tsconfig.tsbuildinfo'));
  for (const pkg of inventory.packages ?? []) {
    const packageRoot = path.join(root, pkg.path);
    targets.add(path.join(packageRoot, 'lib'));
    for (const name of fs.readdirSync(packageRoot)) {
      if (/^tsconfig.*\.tsbuildinfo$/u.test(name)) targets.add(path.join(packageRoot, name));
    }
  }
}

console.log('\n### cleaning', [...targets].map((target) => path.relative(root, target)).join(', '));
for (const target of targets) {
  try {
    fs.rmSync(target, { recursive: true, force: true });
  } catch (error) {
    console.warn(`Cleanup failed for ${path.relative(root, target)}: ${error}`);
  }
}

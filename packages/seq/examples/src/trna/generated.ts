/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { execFileSync } from 'node:child_process';
import { existsSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { derive, derivedModule } from './derive.js';

// Writing `derived-data.ts` (by the record script and the fixtures test's update path): the module as `derive.ts` prints
// it, then formatted by the repository's formatter, so a regeneration leaves `pnpm format:check` clean.

/** The repository's root: the nearest folder above `from` with a `biome.json`. */
function repositoryRoot(from: string): string {
  let dir = from;
  while (!existsSync(join(dir, 'biome.json'))) {
    const up = dirname(dir);
    if (up === dir) throw new Error('no biome.json above the package');
    dir = up;
  }
  return dir;
}

/** Writes the derived data of the seed and 1EHZ's letters to `file`, formatted. */
export function writeDerived(file: string, seed: string, letters: string): void {
  writeFileSync(file, derivedModule(derive(seed, letters)));
  execFileSync('pnpm', ['exec', 'biome', 'format', '--write', file], {
    cwd: repositoryRoot(dirname(file)),
    stdio: 'ignore',
  });
}

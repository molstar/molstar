/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';

// The agent module is opt-in (design/agent.md §6): it is its own package, `@molstar/prism-agent`, which an app gets
// only by depending on it. prism-ui knows it only through its ports (tabs, the command menu's assistant, the canvas's
// agent slot), so neither its manifest nor any of its files (stories and tests included) may name it.
const root = path.resolve(__dirname, '../..');
const Agent = '@molstar/prism-agent';

function* sources(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === 'node_modules' || entry.name === 'lib') continue;
      yield* sources(full);
    } else if (/\.(tsx?|mjs|mdx)$/.test(entry.name) && full !== __filename) yield full;
  }
}

describe('the agent module', () => {
  it('is not a dependency of prism-ui', () => {
    const manifest = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
    const named = ['dependencies', 'devDependencies', 'peerDependencies', 'optionalDependencies'].filter(
      (field) => manifest[field]?.[Agent],
    );
    expect(named).toEqual([]);
  });

  it('is imported by nothing in prism-ui', () => {
    const offenders: string[] = [];
    for (const dir of ['src', '.storybook', 'scripts']) {
      for (const file of sources(path.join(root, dir))) {
        if (fs.readFileSync(file, 'utf8').includes(`'${Agent}`)) offenders.push(path.relative(root, file));
      }
    }
    expect(offenders).toEqual([]);
  });
});

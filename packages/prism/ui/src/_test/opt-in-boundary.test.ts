/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';

// The opt-in modules are their own packages, which an app gets only by depending on them: the agent
// (`@molstar/prism-agent`, design/agent.md §6) and story mode (`@molstar/prism-stories`, design/modes.md). prism-ui
// knows them only through its ports (tabs, modes, the command menu's assistant, the canvas's agent slot), so neither
// its manifest nor any of its files (stories and tests included) may name them.
const root = path.resolve(__dirname, '../..');
const OptIn = ['@molstar/prism-agent', '@molstar/prism-stories'];

function* sources(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === 'node_modules' || entry.name === 'lib') continue;
      yield* sources(full);
    } else if (/\.(tsx?|mjs|mdx)$/.test(entry.name) && full !== __filename) yield full;
  }
}

describe.each(OptIn)('the opt-in module %s', (name) => {
  it('is not a dependency of prism-ui', () => {
    const manifest = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
    const named = ['dependencies', 'devDependencies', 'peerDependencies', 'optionalDependencies'].filter(
      (field) => manifest[field]?.[name],
    );
    expect(named).toEqual([]);
  });

  it('is imported by nothing in prism-ui', () => {
    const offenders: string[] = [];
    for (const dir of ['src', '.storybook', 'scripts']) {
      for (const file of sources(path.join(root, dir))) {
        if (fs.readFileSync(file, 'utf8').includes(`'${name}`)) offenders.push(path.relative(root, file));
      }
    }
    expect(offenders).toEqual([]);
  });
});

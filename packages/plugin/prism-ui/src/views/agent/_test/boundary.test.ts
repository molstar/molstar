/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';

// The agent module is opt-in (design/agent.md §6): an app gets it only by importing it, so nothing else in the package
// may import it; the core knows the agent only through its ports (the command menu's assistant, the canvas slot).
const src = path.resolve(__dirname, '../../..');
const agentModules = [/\/views\/agent\//, /\/model\/agent(\.js)?$/, /\/mocks\/agent(\.js)?$/];

function* sources(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === '_generated' || entry.name === 'node_modules') continue;
      yield* sources(full);
    } else if (/\.tsx?$/.test(entry.name)) yield full;
  }
}

/** Where the agent may be imported: the agent module itself, and stories (they show it). */
const allowed = (file: string) =>
  file.includes(`${path.sep}views${path.sep}agent${path.sep}`) ||
  file.endsWith(`${path.sep}model${path.sep}agent.ts`) ||
  file.endsWith(`${path.sep}mocks${path.sep}agent.ts`) ||
  file.endsWith('.stories.tsx');

describe('the agent module', () => {
  it('is imported by nothing outside it (it is opt-in, not in every build)', () => {
    const offenders: string[] = [];
    for (const file of sources(src)) {
      if (allowed(file)) continue;
      const text = fs.readFileSync(file, 'utf8');
      for (const m of text.matchAll(/(?:from|import)\s*\(?\s*['"]([^'"]+)['"]/g)) {
        const resolved = path.resolve(path.dirname(file), m[1]).split(path.sep).join('/');
        if (m[1].startsWith('.') && agentModules.some((re) => re.test(resolved)))
          offenders.push(`${path.relative(src, file)} imports ${m[1]}`);
      }
    }
    expect(offenders).toEqual([]);
  });
});

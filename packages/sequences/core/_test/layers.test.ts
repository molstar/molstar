/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The core's folder graph: what each folder of src/ may import of the others, as the packages it replaced declared
// it. Inside a folder imports are relative; across folders they name the entry (`@molstar/seq-core/model`), never a path
// into another folder. The core has no third-party dependencies: besides itself, Mol*'s packages (not in its leaves,
// query/ and links/) and Node's built-ins (the CLIs, the tile server). The model's own rule (the query language and the matchers from model/query only) is layers.grit's.

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const SRC = fileURLToPath(new URL('../src', import.meta.url));
/** `file`'s path in `SRC`, with `/` on every platform */
const inSrc = (file: string) => relative(SRC, file).split(sep).join('/');

/** Each folder and the folders it may import. */
const ALLOWED: Record<string, string[]> = {
  query: [],
  algorithms: [],
  model: ['query', 'algorithms'],
  translate: ['model'],
  format: ['model', 'algorithms', 'query'],
  builder: ['format'],
  links: ['format'],
  'stream/tiles': ['format'],
  'stream/gmod': ['format'],
};
/** The leaves import no Mol* package (platform-free, the language and the command links); the rest may. */
const NO_MOLSTAR = ['query', 'links'];
const FOLDERS = Object.keys(ALLOWED).sort((a, b) => b.length - a.length);
const folderOf = (path: string) => FOLDERS.find((f) => path === f || path.startsWith(`${f}/`));

function files(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    const f = join(dir, name);
    if (statSync(f).isDirectory()) files(f, out);
    else if (/\.tsx?$/.test(name)) out.push(f);
  }
  return out;
}

/** Static imports and re-exports (up to their `from`), bare imports, and dynamic imports. */
const SPECIFIER =
  /^(?:import|export)\s[^;]*?[\s}*]from\s*["']([^"']+)["']|^import\s*["']([^"']+)["']|\bimport\(\s*["']([^"']+)["']/gm;

describe("the core's folders", () => {
  const all = files(SRC);

  it('every file sits in a known folder', () => {
    expect(all.length).toBeGreaterThan(150);
    expect(all.map(inSrc).filter((f) => !folderOf(f))).toEqual([]);
  });

  it('import only what their folder may', () => {
    const wrong: string[] = [];
    for (const file of all) {
      const at = inSrc(file);
      const own = folderOf(at) as string;
      for (const m of readFileSync(file, 'utf8').matchAll(SPECIFIER)) {
        const spec = m[1] ?? m[2] ?? m[3];
        if (spec.startsWith('.')) {
          const target = folderOf(inSrc(resolve(dirname(file), spec)));
          if (target !== own) wrong.push(`${at}: ${spec} leaves ${own}/ (name the entry instead)`);
        } else if (spec.startsWith('@molstar/seq-core/')) {
          const target = folderOf(spec.slice('@molstar/seq-core/'.length));
          if (target === own) wrong.push(`${at}: ${spec} is its own folder (import it relatively)`);
          else if (!target || !ALLOWED[own].includes(target)) wrong.push(`${at}: ${spec} is not ${own}'s to import`);
        } else if (spec.startsWith('@molstar/') && !spec.startsWith('@molstar/seq-')) {
          if (NO_MOLSTAR.includes(own)) wrong.push(`${at}: ${spec} (${own}/ is a leaf: no Mol* package)`);
        } else if (!spec.startsWith('node:'))
          wrong.push(`${at}: ${spec} (the core depends on Mol*'s packages and Node's built-ins only)`);
      }
    }
    expect(wrong).toEqual([]);
  });

  it('each entry names what it exports (no `export *`): a name exported is API a host may use', () => {
    const pkg = JSON.parse(readFileSync(fileURLToPath(new URL('../package.json', import.meta.url)), 'utf8'));
    const entries = Object.values(pkg.exports as Record<string, { 'molstar-src': string } | null>)
      .map((t) => t?.['molstar-src'])
      .filter((t): t is string => typeof t === 'string' && t.endsWith('.ts'));
    expect(entries.length).toBeGreaterThan(10);
    const star = entries.filter((t) => /^export \* from/m.test(readFileSync(resolve(SRC, '..', t), 'utf8')));
    expect(star).toEqual([]);
  });
});

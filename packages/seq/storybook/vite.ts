/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync, realpathSync, statSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { join, normalize, relative, sep } from 'node:path';

// What the host Storybook (storybook/packages.ts) takes from this package: a Vite plugin that serves the examples'
// fixtures (`packages/seq/examples/fixtures/`) at `/fixtures/`, so a view's sources, named by relative location
// (`fixtures/p53/tp53-gene.seq`), are read by the page the way Node reads them from the package folder. In the
// development server a middleware; in a static build the files are emitted as assets under `fixtures/`.

const FixturesRoot = fileURLToPath(new URL('../examples/fixtures', import.meta.url));
const Prefix = '/fixtures/';

/** The file a URL path names under `root`, or undefined: not under the prefix, not a file, or out of the folder. */
export function fixtureFile(urlPath: string, root: string = FixturesRoot): string | undefined {
  // the query and the fragment are not the path: cut before decoding
  const path = urlPath.split(/[?#]/, 1)[0];
  if (!path.startsWith(Prefix)) return undefined;
  let rel: string;
  try {
    rel = decodeURIComponent(path.slice(Prefix.length));
  } catch {
    return undefined;
  }
  if (rel.includes('\0')) return undefined;
  const file = normalize(join(root, rel));
  if (file !== root && !file.startsWith(root + sep)) return undefined;
  try {
    if (!statSync(file).isFile()) return undefined;
    // the real path: a link in the folder may lead anywhere
    const real = realpathSync(file);
    const realRoot = realpathSync(root);
    return real.startsWith(realRoot + sep) ? file : undefined;
  } catch {
    return undefined;
  }
}

function listFiles(dir: string): string[] {
  const out: string[] = [];
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) out.push(...listFiles(p));
    else out.push(p);
  }
  return out;
}

interface Middlewares {
  use(fn: (req: { url?: string; method?: string }, res: ServerResponse, next: () => void) => void): void;
}

interface ServerResponse {
  setHeader(name: string, value: string): void;
  end(body: Uint8Array | string): void;
}

/** The plugin, to register in `storybook/packages.ts` as this package's `plugins`. */
export function seqExampleFixtures(root: string = FixturesRoot) {
  return {
    name: 'seqx-example-fixtures',
    configureServer(server: { middlewares: Middlewares }) {
      server.middlewares.use((req, res, next) => {
        // GET and HEAD only: any other method is the server's to answer
        if (req.method !== undefined && req.method !== 'GET' && req.method !== 'HEAD') return next();
        const file = req.url ? fixtureFile(req.url, root) : undefined;
        if (!file) return next();
        res.setHeader('Content-Type', 'text/plain; charset=utf-8');
        res.end(req.method === 'HEAD' ? '' : readFileSync(file));
      });
    },
    generateBundle(this: { emitFile(file: { type: 'asset'; fileName: string; source: Uint8Array }): void }) {
      for (const file of listFiles(root))
        this.emitFile({
          type: 'asset',
          fileName: `fixtures/${relative(root, file).split(sep).join('/')}`,
          source: readFileSync(file),
        });
    },
  };
}

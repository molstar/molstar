/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the host Storybook (storybook/packages.ts) takes from this package: two Vite plugins. The stories and docs pages are
// found by the host's globs (every `*.stories.tsx` and `*.mdx` under `src/`); its preview imports the decorator that
// mounts the React stories (`@molstar/seq-storybook/preview`) and the tokens (`@molstar/seq-storybook/tokens.css`).

import { fileURLToPath } from 'node:url';
import type { Plugin } from 'vite';

const src = fileURLToPath(new URL('./src', import.meta.url));

/**
 * The Seq* stories compile the format's acceptance cases on the format's test fixtures, which read the Seq* JSON
 * fixtures (packages/sequences/core/_test/fixtures) with `node:fs`: for those test modules only, `node:fs` is a small
 * module serving the bundled files (src/_stories/support/seqstar-fixture-fs.ts).
 */
export function seqstarFixtureFs() {
  const shim = `${src}/_stories/support/seqstar-fixture-fs.ts`;
  return {
    name: 'seqstar-fixture-fs',
    enforce: 'pre' as const,
    resolveId(id: string, importer?: string) {
      return id === 'node:fs' && importer && /\/packages\/sequences\/[^/]+\/_test\//.test(importer) ? shim : null;
    },
  };
}

/**
 * The Nightingale elements register themselves when their module runs (`@customElement`), and their packages declare
 * `"sideEffects": false`: a production build drops the bare `import '@nightingale-elements/…'` lines of
 * `@molstar/seq-nightingale`'s mount, and the elements never define (empty tracks). In a Storybook build the modules of
 * those packages count as having side effects.
 */
export function nightingaleSideEffects(): Plugin {
  return {
    name: 'seqstar-nightingale-side-effects',
    transform: {
      filter: { id: /[\\/]@nightingale-elements[\\/]/ },
      handler: (code) => ({ code, map: null, moduleSideEffects: true }),
    },
  };
}

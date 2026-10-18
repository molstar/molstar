/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A resolve hook for the hand-run generators (record-gmod, make-svsx-fixture, font-metrics, record-fixtures): the
// sources import each other with `.js` (NodeNext, as tsc and the bundlers read them), and Node's type stripping runs
// the `.ts` beside it. A relative `.js` specifier whose file does not exist resolves to the sibling `.ts` or `.tsx`.
// Register it with `--import`, add `--conditions=molstar-src` so workspace packages resolve to their sources, and
// `--experimental-transform-types` (Mol*'s sources use namespaces, which type stripping alone rejects):
//   node --conditions=molstar-src --experimental-transform-types --import ./packages/sequences/scripts/ts-resolve.ts <script>.ts
import { existsSync } from 'node:fs';
import { registerHooks } from 'node:module';
import { fileURLToPath } from 'node:url';

registerHooks({
  resolve(specifier, context, nextResolve) {
    if (!/^\.\.?\/.*\.js$/.test(specifier) || !context.parentURL?.startsWith('file:')) {
      return nextResolve(specifier, context);
    }
    const js = new URL(specifier, context.parentURL);
    if (!existsSync(fileURLToPath(js))) {
      for (const ext of ['.ts', '.tsx']) {
        const sibling = new URL(specifier.replace(/\.js$/, ext), context.parentURL);
        if (existsSync(fileURLToPath(sibling))) return nextResolve(sibling.href, context);
      }
    }
    return nextResolve(specifier, context);
  },
});

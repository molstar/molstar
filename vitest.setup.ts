/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, vi } from 'vitest';

// The root scopes of the DOM runtimes (Prism's `dom/scope.ts`, `@molstar/ui/dom`) a file leaves are disposed when it
// ends, and the document emptied: the worker's next files share the document, and a view left listening on it wakes up
// on their events (it renders, attaches listeners, sets state). After the file, not each test: a file keeps views across
// its tests (a viewer made in `beforeAll`), and makes root scopes inside them. `MOLSTAR_REPORT_ROOTS=1` lists where the
// roots left were made.
const TestRootScopes = Symbol.for('molstar.test.rootScopes');
const rootScopes = new Map<{ dispose(): void }, string | undefined>();
const reportRoots = !!process.env.MOLSTAR_REPORT_ROOTS;
(globalThis as Record<symbol, unknown>)[TestRootScopes] = {
  add: (scope: { dispose(): void }) => rootScopes.set(scope, reportRoots ? new Error().stack : undefined),
  delete: (scope: { dispose(): void }) => rootScopes.delete(scope),
};

// Files share their worker's module registry (vitest.config.ts): what a file changed of it is put back when the file
// ends, so the next file in the worker starts as it would alone. Mocked modules are not undone here: a file that mocks
// or resets modules runs in a registry of its own.
afterAll(() => {
  if (rootScopes.size) {
    if (reportRoots) {
      const sites = [...rootScopes.values()].map((s) =>
        s
          ?.split('\n')
          .find((l) => l.includes('.test.'))
          ?.trim(),
      );
      // to stderr: the console of a passing file is not shown
      process.stderr.write(`${rootScopes.size} root scopes left, disposed:\n${sites.join('\n')}\n`);
    }
    // newest first, each on its own: a failing cleanup does not keep the others from running
    for (const scope of [...rootScopes.keys()].reverse()) {
      try {
        scope.dispose();
      } catch (e) {
        console.error(e);
      }
    }
    rootScopes.clear();
  }
  if (typeof document !== 'undefined') document.body?.replaceChildren();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  vi.unstubAllEnvs();
  vi.useRealTimers();
});

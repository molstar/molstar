/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaultClientConditions, defineConfig } from 'vite';

// Workspace packages resolve to their sources through the `molstar-src` export condition (as in the tests and the
// other builds), so no package needs to be built first.
export default defineConfig({
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  build: { outDir: 'build', emptyOutDir: true, target: 'es2022', sourcemap: false },
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaultClientConditions, defineConfig } from 'vite';

// Workspace packages resolve to their sources through the `molstar-src` export condition (as in the tests and the
// other builds), so no package needs to be built first. `vite build --mode single` makes the single-file build (every
// lazy part inlined into one JS file, as the Prism viewer app's `prism-viewer.js`), into `build/single`.
export default defineConfig(({ mode }) => ({
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  build: {
    outDir: mode === 'single' ? 'build/single' : 'build',
    emptyOutDir: true,
    target: 'es2022',
    sourcemap: false,
    ...(mode === 'single' ? { rolldownOptions: { output: { codeSplitting: false } } } : {}),
  },
}));

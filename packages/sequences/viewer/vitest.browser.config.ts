/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fileURLToPath } from 'node:url';
import { playwright } from '@vitest/browser-playwright';
import { defineConfig } from 'vitest/config';
import { browserTests } from '../vitest.shared.js';

/** The Seq* fixtures read their JSON with `node:fs`; in the browser those test modules get the bundled copies. */
const shim = fileURLToPath(new URL('./_test/canvas/fixture-fs.ts', import.meta.url));
const fixtureFs = {
  name: 'seqstar-fixture-fs',
  enforce: 'pre' as const,
  resolveId(id: string, importer?: string) {
    return id === 'node:fs' && importer && /\/lib\/seqstar\/[^/]+\/_test\//.test(importer) ? shim : null;
  },
};

export default defineConfig({ ...browserTests(playwright), plugins: [fixtureFs] });

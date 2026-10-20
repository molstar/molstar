/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';
import { fileURLToPath } from 'node:url';

/** The repository root (the bundle runs from `tests/build/pipeline`) */
export const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../..');

export const SOURCE_1TUP = 'packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures/1tup.bcif';

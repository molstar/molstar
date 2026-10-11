/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { runHarnessConformance } from '@molstar/ai-conformance';
import { createSqliteStorage } from '../storage/sqlite.js';

// The whole suite against SQLite files, one fresh file per case.
const dir = mkdtempSync(join(tmpdir(), 'harness-sqlite-'));
let n = 0;
const opened: { close(): void }[] = [];
afterAll(() => {
  for (const storage of opened) storage.close();
  rmSync(dir, { recursive: true, force: true });
});
runHarnessConformance({
  name: 'sqlite',
  storage: () => {
    const storage = createSqliteStorage(join(dir, `case-${++n}.sqlite`));
    opened.push(storage);
    return storage;
  },
});

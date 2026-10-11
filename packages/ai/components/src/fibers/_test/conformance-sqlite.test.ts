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
import { runFiberConformance } from '@molstar/ai-conformance/fibers';
import { createFiberHost, type FiberTable } from '@molstar/ai-fibers';
import { createSqliteRows } from '../sqlite.js';

// The fibers conformance against SQLite files, one fresh file per case. A restart is a second connection to the file.
const dir = mkdtempSync(join(tmpdir(), 'fibers-sqlite-'));
let n = 0;
const opened: { close(): void }[] = [];
afterAll(() => {
  for (const rows of opened) rows.close();
  rmSync(dir, { recursive: true, force: true });
});
runFiberConformance({
  name: 'sqlite',
  durable: true,
  create: ({ clock }) => {
    const path = join(dir, `case-${++n}.sqlite`);
    const open = () => {
      const rows = createSqliteRows<FiberTable>(path);
      opened.push(rows);
      return rows;
    };
    return {
      fibers: createFiberHost({ rows: open(), clock }),
      reopen: () => createFiberHost({ rows: open(), clock }),
    };
  },
});

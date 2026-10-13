/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import 'fake-indexeddb/auto';
import { runFiberConformance } from '@molstar/ai-conformance/fibers';
import { createFiberHost, type FiberTable } from '@molstar/ai-fibers';
import { createIndexedDbRows } from '../indexeddb.js';

// The fibers conformance against IndexedDB (fake-indexeddb in Node), one fresh database per case. A restart is a
// second handle to the same database, as a second tab is.
jest.setTimeout(30_000);
let n = 0;
const opened: { close(): void }[] = [];
afterAll(() => {
  for (const rows of opened) rows.close();
});
runFiberConformance({
  name: 'indexeddb',
  durable: true,
  create: ({ clock }) => {
    const name = `fibers-case-${++n}`;
    const open = () => {
      const rows = createIndexedDbRows<FiberTable>(name);
      opened.push(rows);
      return rows;
    };
    return {
      fibers: createFiberHost({ rows: open(), clock }),
      reopen: () => createFiberHost({ rows: open(), clock }),
    };
  },
});

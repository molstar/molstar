/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, expect, it, vi } from 'vitest';
import 'fake-indexeddb/auto';
import { runHarnessConformance } from '@molstar/ai-conformance';
import { createIndexedDbStorage } from '../storage/indexeddb.js';

// The whole suite against IndexedDB (fake-indexeddb in Node), one fresh database per case.
// fake-indexeddb is slow: the long-history cases (2,000 appends) run close to Vitest's default of 5 s. Some cases wait a
// fixed number of rounds for work that includes WebCrypto hashing (resolved off the event loop), which this backend can
// outlast under load, so they are retried here (logged) until the suite waits on conditions (.v6/ai/plan.md).
vi.setConfig({ testTimeout: 30_000, retry: 2 });
let n = 0;
const opened: { close(): void }[] = [];
afterAll(() => {
  for (const storage of opened) storage.close();
});
runHarnessConformance({
  name: 'indexeddb',
  storage: () => {
    const storage = createIndexedDbStorage(`harness-case-${++n}`);
    opened.push(storage);
    return storage;
  },
});

it('refuses a database of the same name and version made by something else', async () => {
  const other = indexedDB.open('not-a-harness', 1);
  other.onupgradeneeded = () => other.result.createObjectStore('agent');
  await new Promise((resolve) => (other.onsuccess = resolve));
  other.result.close();
  const storage = createIndexedDbStorage('not-a-harness');
  opened.push(storage);
  await expect(storage.transact('h', async () => undefined)).rejects.toThrow(
    'IndexedDB database "not-a-harness" holds something else (no meta, records, events, blobs store)',
  );
});

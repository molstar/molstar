/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createThreadHost, createMemoryRows, type ThreadTable } from '../index.js';
import { runThreadConformance } from './index.js';

// The threads conformance against the memory rows: a restart is a second host over the same rows.
runThreadConformance({
  name: 'memory',
  durable: true,
  create: ({ clock }) => {
    const rows = createMemoryRows<ThreadTable>();
    return { threads: createThreadHost({ rows, clock }), reopen: () => createThreadHost({ rows, clock }) };
  },
});

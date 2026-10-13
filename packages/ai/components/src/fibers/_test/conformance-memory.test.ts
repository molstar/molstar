/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createFiberHost, createMemoryRows, type FiberTable } from '@molstar/ai-fibers';
import { runFiberConformance } from '@molstar/ai-conformance/fibers';

// The fibers conformance against the memory rows: a restart is a second host over the same rows.
runFiberConformance({
  name: 'memory',
  durable: true,
  create: ({ clock }) => {
    const rows = createMemoryRows<FiberTable>();
    return { fibers: createFiberHost({ rows, clock }), reopen: () => createFiberHost({ rows, clock }) };
  },
});

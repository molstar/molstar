/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { checkLinear, scaling } from '@molstar/test-harness/scaling';
import { describe, expect, it } from 'vitest';
import { LogModel } from '../../log/log-model.js';
import { StatusModel, type StatusTimers } from '../status-model.js';

// A burst of messages and of task progress grows about linearly (status §1): n and 4n posts into the status's log
// (one flush), n and 4n progress events over a few tasks (each a keyed change, no copy of the log).

const timers: StatusTimers = { setTimeout: () => 0, clearTimeout: () => {} };

describe('the status, by size', () => {
  it('a burst of n posts', () => {
    const s = scaling(
      {
        setup: (n) => {
          const log = new LogModel({ size: 4 * n, schedule: () => {} });
          return { status: new StatusModel({ log, timers }), n };
        },
        run: ({ status, n }) => {
          for (let i = 0; i < n; i++) status.post(`Message ${i}`, i % 9 ? 'info' : 'warning');
          status.log.flush();
        },
        teardown: ({ status, n }) => expect(status.log.entries.value).toHaveLength(n),
      },
      4000,
    );
    checkLinear('status: a burst of n posts', s);
  });

  it('n progress events over four tasks, one pill followed', () => {
    const s = scaling(
      {
        setup: (n) => {
          const status = new StatusModel({ timers });
          const stop = status.task('t0').listen!(() => {});
          return { status, n, stop };
        },
        run: ({ status, n }) => {
          for (let i = 0; i < n; i++) {
            status.setTasks([0, 1, 2, 3].map((t) => ({ id: `t${t}`, label: `Task ${t}`, progress: (i % 100) / 100 })));
          }
        },
        teardown: ({ status, stop }) => {
          stop();
          expect(status.taskIds.value).toHaveLength(4);
        },
      },
      4000,
    );
    checkLinear('status: n progress events', s);
  });
});

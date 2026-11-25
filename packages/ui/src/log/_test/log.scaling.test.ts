/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { checkLinear, quadraticControl, scaling } from '@molstar/test-harness/scaling';
import { describe, expect, it } from 'vitest';
import { LogModel, LogProblems } from '../log-model.js';

// A burst of entries is one change (log §2): writing n and 4n entries grows about linearly, in one flush (n entries
// held) and in many (a bounded log written a few at a time). The quadratic control is caught.

const noSchedule = () => {};

describe('the log, by size', () => {
  it('a burst of n entries in one flush', () => {
    const s = scaling(
      {
        setup: (n) => ({ log: new LogModel({ size: 4 * n, schedule: noSchedule }), n }),
        run: ({ log, n }) => {
          for (let i = 0; i < n; i++)
            log.write({ text: `Created ${i}`, kind: i % 7 ? 'info' : 'warning', source: 'S' });
          log.flush();
        },
        teardown: ({ log, n }) => expect(log.entries.value).toHaveLength(n),
      },
      5000,
    );
    checkLinear('log: a burst of n entries', s);
  });

  it('n entries ten a flush into a bounded log, its filter followed', () => {
    const s = scaling(
      {
        setup: (n) => {
          const log = new LogModel({ size: 300, schedule: noSchedule });
          log.setFilter(LogProblems);
          const stop = log.filtered.listen!(() => {});
          return { log, n, stop };
        },
        run: ({ log, n }) => {
          for (let i = 0; i < n; i++) {
            log.write({ text: `Updated ${i}`, kind: i % 5 ? 'info' : 'error', source: `S${i % 3}` });
            if (i % 10 === 9) log.flush();
          }
          log.flush();
        },
        teardown: ({ log, stop }) => {
          stop();
          expect(log.entries.value).toHaveLength(300);
        },
      },
      // Large enough that the base run is ~10 ms: at 4000 it was ~2 ms, and a GC in the larger run broke the ratio.
      16000,
    );
    checkLinear('log: n entries into a bounded log', s);
  });

  it('catches the quadratic control', () => {
    checkLinear('control', scaling(quadraticControl, 2000), false);
  });
});

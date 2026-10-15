/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { RunRecord } from '../../api/run.js';

/** Just what the record reads of the API (its finished runs), so tests can write hundreds of entries at once. */
export function fakeApi() {
  const finished = new Subject<RunRecord>();
  let id = 0;
  const api = { events: { finished }, describe: () => undefined } as unknown as ViewerApi;
  /** A person's change, finished. */
  const run = (did = `Change ${id + 1}`, origin: RunRecord['origin'] = { kind: 'person' }) =>
    finished.next({
      id: ++id,
      operation: 'styles.apply',
      title: 'Apply a style',
      values: {},
      origin,
      started: id,
      result: { ok: true, did } as RunRecord['result'],
      ms: 1,
      step: id,
    });
  return { api, run };
}

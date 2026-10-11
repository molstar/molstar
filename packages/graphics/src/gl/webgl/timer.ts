/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GLRenderingContext } from './compat.js';
import type { GraphicsStats } from '../stats.js';
import type { WebGLExtensions } from './extensions.js';
import { type BackendTimer, createGraphicsTimer, type TimerSource } from '../timer.js';

/** One `TIME_ELAPSED` query per pass (queries can't nest, passes are sequential) */
function createTimerSource(gl: GLRenderingContext, extensions: WebGLExtensions): TimerSource | undefined {
  const dtq = extensions.disjointTimerQuery;
  if (!dtq) return undefined;

  const queries = new Map<number, WebGLQuery>();
  let nextId = 0;
  let current = -1;

  return {
    begin: () => {
      const query = dtq.createQuery();
      if (!query) return -1;
      dtq.beginQuery(dtq.TIME_ELAPSED, query);
      const id = nextId++;
      queries.set(id, query);
      current = id;
      return id;
    },
    end: () => {
      if (current === -1) return;
      dtq.endQuery(dtq.TIME_ELAPSED);
      current = -1;
    },
    poll: () => {
      const finished: [number, number | undefined][] = [];
      queries.forEach((query, id) => {
        if (id === current) return;

        const available = dtq.getQueryParameter(query, dtq.QUERY_RESULT_AVAILABLE);
        const disjoint = gl.getParameter(dtq.GPU_DISJOINT);

        if (available && !disjoint) {
          finished.push([id, dtq.getQueryParameter(query, dtq.QUERY_RESULT) as number]);
        } else if (disjoint) {
          finished.push([id, undefined]);
        }

        if (available || disjoint) {
          dtq.deleteQuery(query);
          queries.delete(id);
        }
      });
      return finished;
    },
    clear: () => {
      if (current !== -1) dtq.endQuery(dtq.TIME_ELAPSED);
      current = -1;
      queries.forEach((query) => dtq.deleteQuery(query));
      queries.clear();
    },
  };
}

export function createTimer(
  gl: GLRenderingContext,
  extensions: WebGLExtensions,
  stats: GraphicsStats,
  options?: { avgCount: number },
): BackendTimer {
  return createGraphicsTimer(createTimerSource(gl, extensions), stats, options);
}

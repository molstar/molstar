/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, ApiResult } from './area.js';

// `camera` (operations.md §2): view state, never a history step (sequence.md §7.3). `reset` now; `zoom`, `orient` and
// `frame` with targets (plan phase 2).

export class CameraApi extends ApiArea {
  readonly id = 'camera';

  readonly reset = this.operation({
    id: 'camera.reset',
    title: 'Reset the camera',
    params: {
      animate: PD.Boolean(true, { description: 'A short transition instead of a cut' }),
    },
    history: false,
    modules: ['canvas'],
    async run({ animate }, ctx) {
      const reduced = ctx.modules.has('page') && prefersReducedMotion();
      ctx.plugin.managers.camera.reset(undefined, animate && !reduced ? 250 : 0);
      return ApiResult.ok('Camera reset.', { view: true });
    },
  });
}

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The quick allocation tier (`test:perf`) of a hover: a highlight moved from one sphere to the next and a Canvas3D tick
// (the marks resolved, the marking redrawn over the last image), as a pointer moving over objects causes. On WebGL (the
// GL shim) and WebGPU (the stub device); no GPU. Each with its positive control.

import { it } from 'vitest';
import { quickCase } from '@molstar/test-harness/quick';
import { createHoverCase } from './hover-cases.js';

for (const backend of ['webgl', 'webgpu'] as const) {
  it(`hover tick (${backend})`, async ({ task }) => {
    const { canvas3d, context, hover } = await createHoverCase(backend);
    quickCase(task, `hover: highlight moved and a tick (${backend})`, (i) => hover(i));
    canvas3d.dispose();
    context.dispose();
  });
}

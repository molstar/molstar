/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRequire } from 'node:module';
import { createContext } from '../webgl/context.js';

const require = createRequire(import.meta.url);

// The tests that need the native backend (`gl`, `pnpm native:install`) are medium (AGENTS.md, Tests): the medium CI
// job installs it and runs them under Xvfb; without it they skip.

export function getGLContext(width: number, height: number) {
  const gl = require('gl')(width, height, {
    alpha: true,
    depth: true,
    stencil: true,
    premultipliedAlpha: true,
    preserveDrawingBuffer: true,
    antialias: true,
  });
  return createContext(gl);
}

export function tryGetGLContext(width: number, height: number, required?: { fragDepth?: boolean }) {
  try {
    const ctx = getGLContext(width, height);
    if (required?.fragDepth && !ctx.capabilities.fragDepth) return;
    return ctx;
  } catch (e) {
    return;
  }
}

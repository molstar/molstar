/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { getDepthRangeCell } from '../depth-range.js';
import type { UniformBlockWriter } from '../uniform-layout.js';
import type { GlobalBindings } from '../render-item.js';
import type { RenderableSchema } from '../renderable/schema.js';
import type { Textures } from '../texture.js';
import type { UniformBuffer, UniformsList } from '../uniform.js';
import { getUniformBlockLayout } from '../uniform-layout.js';
import type { WebGLContext } from './context.js';
import type { Program } from './program.js';

export interface WebGLGlobalBindings extends GlobalBindings {
  /** Uses `program` and applies the bindings if it changed or they were invalidated */
  apply: (program: Program) => void;
}

/** Global bindings created by the WebGL backend, see `isWebGLGlobalBindings` */
const WebGLGlobalBindingsSet = new WeakSet<GlobalBindings>();

export function isWebGLGlobalBindings(globals: GlobalBindings): globals is WebGLGlobalBindings {
  return WebGLGlobalBindingsSet.has(globals);
}

/** Narrows to `WebGLGlobalBindings`, throws for global bindings of another backend */
export function asWebGLGlobalBindings(globals: GlobalBindings): WebGLGlobalBindings {
  if (!isWebGLGlobalBindings(globals)) {
    throw new Error('expected WebGL global bindings, got global bindings of another backend');
  }
  return globals;
}

export function createGlobalBindings(
  ctx: WebGLContext,
  schema: RenderableSchema,
  uniforms: UniformsList,
  textures: Textures,
): WebGLGlobalBindings {
  const { state, resources } = ctx;
  const globalTextureNames = Object.keys(schema).filter((k) => schema[k].type === 'texture');
  /** the global textures of `schema` not in `textures` (filled by passes as they run) */
  const getMissingTextures = () => {
    const missing: string[] = [];
    for (const k of globalTextureNames) if (!textures.some((t) => t[0] === k)) missing.push(k);
    return missing;
  };
  // the projections mapped to the context's clip depth
  uniforms = uniforms.map(([k, v]) => [k, getDepthRangeCell(schema, k, v, ctx.depthRangeConversion)]);

  // on WebGL2 most global uniforms are backed by a uniform buffer,
  // the rest (arrays, uDrawId) are still set as classic uniforms
  const layout = ctx.uniformBlocks ? getUniformBlockLayout('global', schema) : undefined;
  let block: { buffer: UniformBuffer; write: UniformBlockWriter } | undefined;
  if (layout) {
    const buffer = resources.uniformRingBuffer(layout);
    block = { buffer, write: buffer.createWriter(uniforms) };
  }
  const classicUniforms = block
    ? uniforms.filter(([k]) => {
        const spec = schema[k];
        return !(spec?.type === 'uniform' && spec.block === 'global');
      })
    : uniforms;

  let needsUpdate = true;
  const capPass = uniforms.find(([k]) => k === 'uSolidInteriorPass')?.[1] as ValueCell<number> | undefined;
  const capClip = uniforms.find(([k]) => k === 'uSolidInteriorClip')?.[1] as ValueCell<number> | undefined;

  const globals: WebGLGlobalBindings = {
    textures,
    invalidate: () => {
      needsUpdate = true;
      // render items rebind all of their values, the shared texture count may have changed
      state.currentRenderItemId = -1;
    },
    setSolidInteriorCap: (pass: number, clip: number) => {
      if (!capPass || !capClip || (capPass.ref.value === pass && capClip.ref.value === clip)) return;
      ValueCell.update(capPass, pass);
      ValueCell.update(capClip, clip);
      globals.invalidate();
    },
    apply: (program: Program) => {
      if (state.currentProgramId !== program.id) {
        needsUpdate = true;
        program.use();
      }
      if (!needsUpdate) return;

      if (block) {
        block.write();
        block.buffer.sync();
        block.buffer.bind();
      }
      program.setUniforms(classicUniforms);
      program.bindTextures(textures);
      // the program's global samplers this pass does not provide (the DPOIT textures before the first DPOIT pass):
      // unbound, else they sample what their unit last held, possibly a texture drawn into (a feedback loop)
      program.unbindTextures(getMissingTextures());
      needsUpdate = false;
    },
    destroy: () => {
      block?.write.dispose();
      block?.buffer.destroy();
    },
  };
  WebGLGlobalBindingsSet.add(globals);
  return globals;
}

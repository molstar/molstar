/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type RenderableSchema,
  type Values,
  type UnboxedValues,
  UniformSpec,
  TextureSpec,
  DefineSpec,
  type RenderableValues,
} from '../renderable/schema.js';
import { isWebGLContext, type WebGLContext } from '../webgl/context.js';
import { getRegularGrid3dDelta, type RegularGrid3d } from '@molstar/core/math/geometry/common';
import { grid3dTemplate_frag } from '../webgl/shader/util/grid3d-template.frag.js';
import { quad_vert } from '../webgl/shader/quad.vert.js';
import { registerGlslProgram } from '../webgl/programs.js';
import { ShaderCode } from '../shader-code.js';
import { UUID, ValueCell } from '@molstar/core/util';
import { objectForEach } from '@molstar/core/util/object';
import { isUniformValueScalar } from '../uniform.js';
import { QuadSchema, QuadValues } from './util.js';
import { createComputeRenderable } from '../renderable.js';
import { isLittleEndian } from '@molstar/core/util/is-little-endian';
import type { RuntimeContext } from '@molstar/core/task';
import { isTimingMode } from '@molstar/core/util/debug';
import type { GraphicsContext } from '../context.js';
import { RenderState } from '../render-state.js';

export function canComputeGrid3dOnGPU(graphics?: GraphicsContext): graphics is WebGLContext {
  return !!graphics && isWebGLContext(graphics);
}

export interface Grid3DComputeRenderableSpec<S extends RenderableSchema, P, CS> {
  schema: S;
  // indicate which params are loop bounds for WebGL1 compat
  loopBounds?: (keyof S)[];
  utilCode?: string;
  mainCode: string;
  returnCode: string;

  values(params: P, grid: RegularGrid3d): UnboxedValues<S>;

  cumulative?: {
    states(params: P): CS[];
    update(params: P, state: CS, values: Values<S>): void;
    // call gl.readPixes every 'yieldPeriod' states to split the computation
    // into multiple parts, if not set, the computation will be synchronous
    yieldPeriod?: number;
  };
}

const FrameBufferName = 'grid3d-computable' as const;
const Texture0Name = 'grid3d-computable-0' as const;
const Texture1Name = 'grid3d-computable-1' as const;

const SchemaBase = {
  ...QuadSchema,
  uDimensions: UniformSpec('v3'),
  uMin: UniformSpec('v3'),
  uDelta: UniformSpec('v3'),
  uWidth: UniformSpec('f'),
  uLittleEndian: UniformSpec('b'),
};

const CumulativeSumSchema = {
  tCumulativeSum: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
};

/** fragment shader code, declarations are generated from `getGrid3dComputeSchema` */
export function getGrid3dComputeFrag(spec: Grid3DComputeRenderableSpec<any, any, any>) {
  return grid3dTemplate_frag
    .replace('{UTILS}', spec.utilCode ?? '')
    .replace('{MAIN}', spec.mainCode)
    .replace('{RETURN}', spec.returnCode);
}

/** loop bounds are number defines on WebGL1 (constant loop bounds) */
export function getGrid3dComputeSchema(
  spec: Grid3DComputeRenderableSpec<any, any, any>,
  isWebGL2: boolean,
): RenderableSchema {
  const schema: { [k: string]: RenderableSchema[string] } = {
    ...SchemaBase,
    ...(spec.cumulative ? CumulativeSumSchema : {}),
    ...spec.schema,
  };
  if (!isWebGL2) {
    for (const b of spec.loopBounds ?? []) schema[b as string] = DefineSpec('number');
    schema.WEBGL1 = DefineSpec('boolean');
  }
  if (spec.cumulative) schema.CUMULATIVE = DefineSpec('boolean');
  return schema;
}

export function createGrid3dComputeRenderable<S extends RenderableSchema, P, CS>(
  spec: Grid3DComputeRenderableSpec<S, P, CS>,
) {
  const id = UUID.create22();

  objectForEach(spec.schema, (u) => {
    if (u.type === 'uniform' && u.kind.endsWith('[]')) throw new Error('array uniforms are not supported');
  });

  const shader = ShaderCode(id);
  registerGlslProgram(id, quad_vert, getGrid3dComputeFrag(spec));

  return async (ctx: RuntimeContext, webgl: WebGLContext, grid: RegularGrid3d, params: P) => {
    const schema = getGrid3dComputeSchema(spec, webgl.isWebGL2);

    if (!webgl.namedFramebuffers[FrameBufferName]) {
      webgl.namedFramebuffers[FrameBufferName] = webgl.resources.framebuffer();
    }
    const framebuffer = webgl.namedFramebuffers[FrameBufferName];

    if (!webgl.namedTextures[Texture0Name]) {
      webgl.namedTextures[Texture0Name] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    if (spec.cumulative && !webgl.namedTextures[Texture1Name]) {
      webgl.namedTextures[Texture1Name] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }

    const tex = [webgl.namedTextures[Texture0Name], webgl.namedTextures[Texture1Name]];

    const [nx, ny, nz] = grid.dimensions;
    const uWidth = Math.ceil(Math.sqrt(nx * ny * nz));

    const values: UnboxedValues<S & typeof SchemaBase> = {
      uDimensions: grid.dimensions,
      uMin: grid.box.min,
      uDelta: getRegularGrid3dDelta(grid),
      uWidth,
      uLittleEndian: isLittleEndian(),
      ...spec.values(params, grid),
    } as any;

    if (!webgl.isWebGL2) {
      (values as any).WEBGL1 = true;
    }
    if (spec.cumulative) {
      (values as any).tCumulativeSum = tex[0];
      (values as any).CUMULATIVE = true;
    }

    let renderable = webgl.namedComputeRenderables[id];
    let cells: RenderableValues;
    if (renderable) {
      cells = renderable.values as RenderableValues;
      objectForEach(values, (c, k) => {
        const s = schema[k];
        if (s?.type === 'value' || s?.type === 'attribute') return;

        if (!s || !isUniformValueScalar(s.kind as any)) {
          ValueCell.update(cells[k], c);
        } else {
          ValueCell.updateIfChanged(cells[k], c);
        }
      });
    } else {
      cells = {} as any;
      objectForEach(QuadValues, (v, k) => ((cells as any)[k] = v));
      objectForEach(values, (v, k) => ((cells as any)[k] = ValueCell.create(v)));
      renderable = createComputeRenderable(webgl.createComputeRenderItem('triangles', shader, schema, cells), cells);
    }

    const array = new Uint8Array(uWidth * uWidth * 4);
    if (spec.cumulative) {
      if (isTimingMode) webgl.timer.group('Grid3dCompute.renderCumulative');

      const states = spec.cumulative.states(params);

      tex[0].define(uWidth, uWidth);
      tex[1].define(uWidth, uWidth);

      resetGl(webgl, uWidth);

      tex[0].attachFramebuffer(framebuffer, 'color0');
      webgl.beginRenderPass(framebuffer, { label: 'grid3d clear', clear: { color: [0, 0, 0, 0] } });

      tex[1].attachFramebuffer(framebuffer, 'color0');
      webgl.beginRenderPass(framebuffer, { label: 'grid3d clear', clear: { color: [0, 0, 0, 0] } });

      if (spec.cumulative.yieldPeriod && !isTimingMode) {
        await ctx.update({ message: 'Computing...', isIndeterminate: false, current: 0, max: states.length });
      }

      const yieldPeriod = Math.max(1, spec.cumulative.yieldPeriod ?? 1 | 0);

      if (isTimingMode) webgl.timer.group('Grid3dCompute.renderBatch');
      for (let i = 0; i < states.length; i++) {
        ValueCell.update(cells.tCumulativeSum, tex[(i + 1) % 2]);
        tex[i % 2].attachFramebuffer(framebuffer, 'color0');
        webgl.beginRenderPass(framebuffer, { label: 'grid3d' });
        resetGl(webgl, uWidth);
        spec.cumulative.update(params, states[i], cells as any);
        renderable.update();
        renderable.render();

        if (spec.cumulative.yieldPeriod && i !== states.length - 1) {
          if (i % yieldPeriod === yieldPeriod - 1) {
            webgl.waitForGpuCommandsCompleteSync();
            if (isTimingMode) webgl.timer.groupEnd('Grid3dCompute.renderBatch');
            if (isTimingMode) webgl.timer.group('Grid3dCompute.renderBatch');
          }
          if (ctx.shouldUpdate && !isTimingMode) {
            await ctx.update({ current: i + 1 });
          }
        }
      }
      if (isTimingMode) webgl.timer.groupEnd('Grid3dCompute.renderBatch');
      if (isTimingMode) webgl.timer.groupEnd('Grid3dCompute.renderCumulative');
    } else {
      tex[0].define(uWidth, uWidth);
      tex[0].attachFramebuffer(framebuffer, 'color0');
      webgl.beginRenderPass(framebuffer, { label: 'grid3d' });
      resetGl(webgl, uWidth);
      renderable.update();
      renderable.render();
    }

    if (isTimingMode) webgl.timer.group('Grid3dCompute.readPixels');
    webgl.readPixels(0, 0, uWidth, uWidth, array);
    if (isTimingMode) webgl.timer.groupEnd('Grid3dCompute.readPixels');
    return new Float32Array(array.buffer, array.byteOffset, nx * ny * nz);
  };
}

function resetGl(graphics: GraphicsContext, w: number) {
  graphics.setViewport(0, 0, w, w);
  graphics.disableScissor();
  graphics.setRenderState(RenderState.Default);
}

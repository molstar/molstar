/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getRegularGrid3dDelta } from '@molstar/core/math/geometry/common';
import type { RuntimeContext } from '@molstar/core/task';
import { ValueCell } from '@molstar/core/util';
import { isTimingMode } from '@molstar/core/util/debug';
import { isLittleEndian } from '@molstar/core/util/is-little-endian';
import { objectForEach } from '@molstar/core/util/object';
import type { RegularGrid3d } from '@molstar/core/math/geometry/common';
import type { Grid3DComputeRenderableSpec, Grid3dFragmentImplementation } from '../../compute/grid3d.js';
import { QuadSchema, QuadValues } from '../../compute/util.js';
import type { GraphicsContext } from '../../context.js';
import { createComputeRenderable } from '../../renderable.js';
import {
  DefineSpec,
  type RenderableSchema,
  type RenderableValues,
  TextureSpec,
  type UnboxedValues,
  UniformSpec,
} from '../../renderable/schema.js';
import { RenderState } from '../../render-state.js';
import { ShaderCode } from '../../shader-code.js';
import { isUniformValueScalar } from '../../uniform.js';
import { isWebGLContext } from '../context.js';
import { registerGlslProgram } from '../programs.js';
import { quad_vert } from '../shader/quad.vert.js';
import { grid3dTemplate_frag } from '../shader/util/grid3d-template.frag.js';

/**
 * grid3d with fragment shaders (WebGL): a fragment per cell of a square texture, the values packed into rgba8 and read
 * back; registered with each WebGL context, see `createGrid3dComputeRenderable`
 */

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

/** The computation of `spec` (its GLSL program registered as `id`) */
function createWebGLGrid3dCompute<S extends RenderableSchema, P, CS>(
  id: string,
  spec: Grid3DComputeRenderableSpec<S, P, CS>,
) {
  const shader = ShaderCode(id);
  registerGlslProgram(id, quad_vert, getGrid3dComputeFrag(spec));

  return async (ctx: RuntimeContext, graphics: GraphicsContext, grid: RegularGrid3d, params: P) => {
    if (!isWebGLContext(graphics)) throw new Error('not a WebGL context');
    const webgl = graphics;
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

export const WebGLGrid3dFragment: Grid3dFragmentImplementation = {
  isSupported: (graphics) => isWebGLContext(graphics),
  create: createWebGLGrid3dCompute,
};

function resetGl(graphics: GraphicsContext, w: number) {
  graphics.setViewport(0, 0, w, w);
  graphics.disableScissor();
  graphics.setRenderState(RenderState.Default);
}

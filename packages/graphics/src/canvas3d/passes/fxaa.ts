/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Viewport } from '../camera/util.js';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';

export const FxaaParams = {
  edgeThresholdMin: PD.Numeric(
    0.0312,
    { min: 0.0312, max: 0.0833, step: 0.0001 },
    { description: 'Trims the algorithm from processing darks.' },
  ),
  edgeThresholdMax: PD.Numeric(
    0.063,
    { min: 0.063, max: 0.333, step: 0.001 },
    { description: 'The minimum amount of local contrast required to apply algorithm.' },
  ),
  iterations: PD.Numeric(12, { min: 0, max: 16, step: 1 }, { description: 'Number of edge exploration steps.' }),
  subpixelQuality: PD.Numeric(
    0.3,
    { min: 0.0, max: 1.0, step: 0.01 },
    { description: 'Choose the amount of sub-pixel aliasing removal.' },
  ),
};
export type FxaaProps = PD.Values<typeof FxaaParams>;

export class FxaaPass {
  private readonly renderable: FxaaRenderable;

  constructor(
    private graphics: GraphicsContext,
    input: Texture,
  ) {
    this.renderable = getFxaaRenderable(graphics, input);
  }

  private updateState(viewport: Viewport) {
    const { x, y, width, height } = viewport;
    this.graphics.setRenderState(RenderState.Default);
    this.graphics.setViewport(x, y, width, height);
    this.graphics.setScissor(x, y, width, height);
    this.graphics.clear({ color: [0, 0, 0, 1] });
  }

  setSize(width: number, height: number) {
    ValueCell.update(
      this.renderable.values.uTexSizeInv,
      Vec2.set(this.renderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
    );
  }

  update(input: Texture, props: FxaaProps) {
    const { values } = this.renderable;
    const { edgeThresholdMin, edgeThresholdMax, iterations, subpixelQuality } = props;

    let needsUpdate = false;

    if (values.tColor.ref.value !== input) {
      ValueCell.update(this.renderable.values.tColor, input);
      needsUpdate = true;
    }

    if (values.dEdgeThresholdMin.ref.value !== edgeThresholdMin) needsUpdate = true;
    ValueCell.updateIfChanged(values.dEdgeThresholdMin, edgeThresholdMin);

    if (values.dEdgeThresholdMax.ref.value !== edgeThresholdMax) needsUpdate = true;
    ValueCell.updateIfChanged(values.dEdgeThresholdMax, edgeThresholdMax);

    if (values.dIterations.ref.value !== iterations) needsUpdate = true;
    ValueCell.updateIfChanged(values.dIterations, iterations);

    if (values.dSubpixelQuality.ref.value !== subpixelQuality) needsUpdate = true;
    ValueCell.updateIfChanged(values.dSubpixelQuality, subpixelQuality);

    if (needsUpdate) {
      this.renderable.update();
    }
  }

  render(viewport: Viewport, target: RenderTarget | undefined) {
    this.graphics.beginRenderPass(target ?? 'drawing-buffer', { label: 'fxaa' });
    this.updateState(viewport);
    this.renderable.render();
  }
}

//

export const FxaaSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  uTexSizeInv: UniformSpec('v2'),

  dEdgeThresholdMin: DefineSpec('number'),
  dEdgeThresholdMax: DefineSpec('number'),
  dIterations: DefineSpec('number'),
  dSubpixelQuality: DefineSpec('number'),
};
export const FxaaShaderCode = ShaderCode('fxaa');
type FxaaRenderable = ComputeRenderable<Values<typeof FxaaSchema>>;

function getFxaaRenderable(ctx: GraphicsContext, colorTexture: Texture): FxaaRenderable {
  const width = colorTexture.getWidth();
  const height = colorTexture.getHeight();

  const values: Values<typeof FxaaSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

    dEdgeThresholdMin: ValueCell.create(0.0312),
    dEdgeThresholdMax: ValueCell.create(0.125),
    dIterations: ValueCell.create(12),
    dSubpixelQuality: ValueCell.create(0.3),
  };

  const schema = { ...FxaaSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', FxaaShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

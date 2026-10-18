/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Viewport } from '../camera/util.js';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';

export const CasParams = {
  sharpness: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }),
  denoise: PD.Boolean(true),
};
export type CasProps = PD.Values<typeof CasParams>;

export class CasPass {
  private readonly renderable: CasRenderable;

  constructor(
    private graphics: GraphicsContext,
    input: Texture,
  ) {
    this.renderable = getCasRenderable(graphics, input);
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

  update(input: Texture, props: CasProps) {
    const { values } = this.renderable;
    const { sharpness, denoise } = props;

    let needsUpdate = false;

    if (values.tColor.ref.value !== input) {
      ValueCell.update(this.renderable.values.tColor, input);
      needsUpdate = true;
    }

    ValueCell.updateIfChanged(values.uSharpness, 2 - 2 * Math.pow(sharpness, 0.25));

    if (values.dDenoise.ref.value !== denoise) needsUpdate = true;
    ValueCell.updateIfChanged(values.dDenoise, denoise);

    if (needsUpdate) {
      this.renderable.update();
    }
  }

  render(b: PassBuilder, viewport: Viewport, target: RenderTarget | undefined) {
    b.begin(target ?? 'drawing-buffer', { label: 'cas' });
    if (b.drawing) this.updateState(viewport);
    b.quad(this.renderable, RenderState.Default);
  }
}

//

export const CasSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  uTexSizeInv: UniformSpec('v2'),

  uSharpness: UniformSpec('f'),
  dDenoise: DefineSpec('boolean'),
};
export const CasShaderCode = ShaderCode('cas');
type CasRenderable = ComputeRenderable<Values<typeof CasSchema>>;

function getCasRenderable(ctx: GraphicsContext, colorTexture: Texture): CasRenderable {
  const width = colorTexture.getWidth();
  const height = colorTexture.getHeight();

  const values: Values<typeof CasSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

    uSharpness: ValueCell.create(0.5),
    dDenoise: ValueCell.create(true),
  };

  const schema = { ...CasSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', CasShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

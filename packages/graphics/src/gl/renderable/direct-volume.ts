/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { WebGLContext } from '../webgl/context.js';
import { createGraphicsRenderItem, type Transparency } from '../webgl/render-item.js';
import { AttributeSpec, type Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ElementsSpec, DefineSpec, type InternalValues, GlobalTextureSchema, BaseSchema, ValueSpec, type GlobalDefineValues, type GlobalDefines, GlobalDefineSchema } from './schema.js';
import { DirectVolumeShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const DirectVolumeSchema = {
    ...BaseSchema,

    aPosition: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),

    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uMaxSteps: UniformSpec('i'),
    uStepScale: UniformSpec('f'),
    uJumpLength: UniformSpec('f'),
    uTransform: UniformSpec('m4'),
    uGridDim: UniformSpec('v3'),
    tTransferTex: TextureSpec('image-uint8', 'alpha', 'ubyte', 'linear'),
    uTransferScale: UniformSpec('f', 'material'),

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    uGridTexDim: UniformSpec('v3'),
    tGridTex: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uGridStats: UniformSpec('v4'), // [min, max, mean, sigma]

    uCellDim: UniformSpec('v3'),
    uCartnToUnit: UniformSpec('m4'),
    uUnitToCartn: UniformSpec('m4'),
    dPackedGroup: DefineSpec('boolean'),
    dAxisOrder: DefineSpec('string', ['012', '021', '102', '120', '201', '210']),

    dIgnoreLight: DefineSpec('boolean'),
    dCelShaded: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),

    meta: ValueSpec('unknown')
};
export type DirectVolumeSchema = typeof DirectVolumeSchema
export type DirectVolumeValues = Values<DirectVolumeSchema>

export function DirectVolumeRenderable(ctx: WebGLContext, id: number, values: DirectVolumeValues, state: RenderableState, materialId: number, transparency: Transparency, globals: GlobalDefines): Renderable<DirectVolumeValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...GlobalDefineSchema, ...InternalSchema, ...DirectVolumeSchema };
    if (!ctx.isWebGL2) {
        // workaround for webgl1 limitation that loop counters need to be `const`
        (schema.uMaxSteps as any) = DefineSpec('number');
    }
const renderValues: DirectVolumeValues & InternalValues & GlobalDefineValues = {
        ...values,
        uObjectId: ValueCell.create(id),
        dLightCount: ValueCell.create(globals.dLightCount),
        dColorMarker: ValueCell.create(globals.dColorMarker),
    };
    const shaderCode = DirectVolumeShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, renderValues, materialId, transparency);
    return createRenderable(renderItem, renderValues, state);
}
/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, Transparency } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, ElementsSpec, DefineSpec, Values, InternalSchema, InternalValues, GlobalTextureSchema, ValueSpec, UniformSpec, GlobalDefineSchema, GlobalDefineValues, GlobalDefines } from './schema';
import { MeshShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const MeshSchema = {
    ...BaseSchema,
    aGroup: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    aNormal: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),
    dVaryingGroup: DefineSpec('boolean'),
    dFlatShaded: DefineSpec('boolean'),
    uDoubleSided: UniformSpec('b', 'material'),
    dFlipSided: DefineSpec('boolean'),
    dIgnoreLight: DefineSpec('boolean'),
    dCelShaded: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
    dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
    uBumpFrequency: UniformSpec('f', 'material'),
    uBumpAmplitude: UniformSpec('f', 'material'),
    uInteriorColor: UniformSpec('v4'),
    uInteriorSubstance: UniformSpec('v4'),
    meta: ValueSpec('unknown')
} as const;
export type MeshSchema = typeof MeshSchema
export type MeshValues = Values<MeshSchema>

export function MeshRenderable(ctx: WebGLContext, id: number, values: MeshValues, state: RenderableState, materialId: number, transparency: Transparency, globals: GlobalDefines): Renderable<MeshValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...GlobalDefineSchema, ...InternalSchema, ...MeshSchema };
    const renderValues: MeshValues & InternalValues & GlobalDefineValues = {
        ...values,
        uObjectId: ValueCell.create(id),
        dLightCount: ValueCell.create(globals.dLightCount),
        dColorMarker: ValueCell.create(globals.dColorMarker),
    };
    const shaderCode = MeshShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, renderValues, materialId, transparency);

    return createRenderable(renderItem, renderValues, state);
}
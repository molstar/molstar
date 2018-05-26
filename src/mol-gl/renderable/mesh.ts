/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, ElementsSpec, DefineSpec, Values } from '../renderable/schema';
import { MeshShaderCode } from '../shader-code';

export const MeshSchema = {
    ...BaseSchema,
    aNormal: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),
    dFlatShaded: DefineSpec('boolean'),
    dDoubleSided: DefineSpec('boolean'),
    dFlipSided: DefineSpec('boolean'),
}
export type MeshSchema = typeof MeshSchema
export type MeshValues = Values<MeshSchema>

export function MeshRenderable(ctx: Context, values: MeshValues, state: RenderableState): Renderable<MeshValues> {
    const schema = { ...GlobalUniformSchema, ...MeshSchema }
    const schaderCode = MeshShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', schaderCode, schema, values)

    return {
        draw: () => {
            renderItem.draw()
        },
        get values () { return values },
        get state () { return state },
        name: 'mesh',
        get program () { return renderItem.program },
        update: () => {
            renderItem.update()
        },
        dispose: () => {
            renderItem.destroy()
        }
    }
}
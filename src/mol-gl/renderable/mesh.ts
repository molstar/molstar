/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, ElementsSpec, DefineSpec, Values, InternalSchema } from '../renderable/schema';
import { MeshShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

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

export function MeshRenderable(ctx: Context, id: number, values: MeshValues, state: RenderableState): Renderable<MeshValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...MeshSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const schaderCode = MeshShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', schaderCode, schema, { ...values, ...internalValues })

    return {
        draw: () => {
            renderItem.draw()
        },
        pick: () => {
            renderItem.pick()
        },
        get values () { return values },
        get state () { return state },
        name: 'mesh',
        get drawProgram () { return renderItem.drawProgram },
        get pickProgram () { return renderItem.pickProgram },
        update: () => {
            renderItem.update()
        },
        dispose: () => {
            renderItem.destroy()
        }
    }
}
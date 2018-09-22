
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';

export type DefineKind = 'boolean' | 'string'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

const shaderCodeId = idFactory()

export interface ShaderCode {
    id: number
    vert: string
    frag: string
}

export function ShaderCode(vert: string, frag: string): ShaderCode {
    return { id: shaderCodeId(), vert, frag }
}

export const PointsShaderCode = ShaderCode(
    require('mol-gl/shader/points.vert'),
    require('mol-gl/shader/points.frag')
)

export const LinesShaderCode = ShaderCode(
    require('mol-gl/shader/lines.vert'),
    require('mol-gl/shader/lines.frag')
)

export const MeshShaderCode = ShaderCode(
    require('mol-gl/shader/mesh.vert'),
    require('mol-gl/shader/mesh.frag')
)

export type ShaderDefines = {
    [k: string]: ValueCell<DefineType>
}

function getDefinesCode (defines: ShaderDefines) {
    if (defines === undefined) return ''
    const lines = []
    for (const name in defines) {
        const define = defines[name]
        const v = define.ref.value
        if (v) {
            if (typeof v === 'string') {
                lines.push(`#define ${name}_${v}`)
            } else {
                lines.push(`#define ${name}`)
            }
        }
    }
    return lines.join('\n') + '\n'
}

export function addShaderDefines(defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const header = getDefinesCode(defines)
    return {
        id: shaderCodeId(),
        vert: `${header}${shaders.vert}`,
        frag: `${header}${shaders.frag}`
    }
}
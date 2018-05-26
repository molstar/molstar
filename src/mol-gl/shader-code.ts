
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';

export type DefineKind = 'boolean' | 'string'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

export interface ShaderCode {
    vert: string
    frag: string
}

export const PointShaderCode: ShaderCode = {
    vert: require('mol-gl/shader/point.vert'),
    frag: require('mol-gl/shader/point.frag')
}

export const MeshShaderCode: ShaderCode = {
    vert: require('mol-gl/shader/mesh.vert'),
    frag: require('mol-gl/shader/mesh.frag')
}

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

export function addShaderDefines(defines: ShaderDefines, shaders: ShaderCode) {
    const header = getDefinesCode(defines)
    return {
        vert: `${header}${shaders.vert}`,
        frag: `${header}${shaders.frag}`
    }
}

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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

type ShaderDefine = (
    'UNIFORM_COLOR' | 'ATTRIBUTE_COLOR' | 'INSTANCE_COLOR' | 'ELEMENT_COLOR' | 'ELEMENT_INSTANCE_COLOR' |
    'UNIFORM_SIZE' | 'ATTRIBUTE_SIZE' |
    'POINT_SIZE_ATTENUATION' |
    'FLAT_SHADED' | 'DOUBLE_SIDED' | 'FLIP_SIDED'
)
export type ShaderDefines = {
    [k in ShaderDefine]?: number|string
}

function getDefinesCode (defines: ShaderDefines) {
    if (defines === undefined) return ''
    const lines = []
    for (const name in defines) {
        const value = defines[ name as keyof ShaderDefines ]
        if (value) {
            lines.push(`#define ${name} ${value}`)
        } else {
            lines.push(`#define ${name}`)
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
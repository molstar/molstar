
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

interface Shaders {
    vert: string
    frag: string
}

export const PointShaders = {
    vert: require('mol-gl/shader/point.vert'),
    frag: require('mol-gl/shader/point.frag')
}

export const MeshShaders = {
    vert: require('mol-gl/shader/mesh.vert'),
    frag: require('mol-gl/shader/mesh.frag')
}

type ShaderDefine = (
    'UNIFORM_COLOR' | 'ATTRIBUTE_COLOR' | 'INSTANCE_COLOR' | 'ELEMENT_COLOR' | 'ELEMENT_INSTANCE_COLOR'
)
export type ShaderDefines = {
    [k in ShaderDefine]?: number|string
}

function getDefines (defines: ShaderDefines) {
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

export function addDefines(defines: ShaderDefines, shaders: Shaders) {
    const header = getDefines(defines)
    return {
        vert: `${header}${shaders.vert}`,
        frag: `${header}${shaders.frag}`
    }
}
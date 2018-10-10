/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';

export type DefineKind = 'boolean' | 'string' | 'number'
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

export const GaussianDensityShaderCode = ShaderCode(
    require('mol-gl/shader/gaussian-density.vert'),
    require('mol-gl/shader/gaussian-density.frag')
)

export const DirectVolumeShaderCode = ShaderCode(
    require('mol-gl/shader/direct-volume.vert'),
    require('mol-gl/shader/direct-volume.frag')
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
            } else if (typeof v === 'number') {
                lines.push(`#define ${name} ${v}`)
            } else if (typeof v === 'boolean') {
                lines.push(`#define ${name}`)
            } else {
                throw new Error('unknown define type')
            }
        }
    }
    return lines.join('\n') + '\n'
}

const glsl300VertPrefix = `#version 300 es
#define attribute in
#define varying out
#define texture2D texture
`

const glsl300FragPrefix = `#version 300 es
#define varying in
out highp vec4 out_FragColor;
#define gl_FragColor out_FragColor
#define gl_FragDepthEXT gl_FragDepth
#define texture2D texture
`

export function addShaderDefines(defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const isGlsl300es = defines.dGlslVersion && defines.dGlslVersion.ref.value === '300es'
    const header = getDefinesCode(defines)
    const vertPrefix = isGlsl300es ? glsl300VertPrefix : ''
    const fragPrefix = isGlsl300es ? glsl300FragPrefix : ''
    return {
        id: shaderCodeId(),
        vert: `${vertPrefix}${header}${shaders.vert}`,
        frag: `${fragPrefix}${header}${shaders.frag}`
    }
}
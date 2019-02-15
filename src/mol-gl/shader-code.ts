/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';
import { WebGLContext } from './webgl/context';

export type DefineKind = 'boolean' | 'string' | 'number'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

const shaderCodeId = idFactory()

export interface ShaderExtensions {
    readonly standardDerivatives: boolean
    readonly fragDepth: boolean
}

export interface ShaderCode {
    readonly id: number
    readonly vert: string
    readonly frag: string
    readonly extensions: ShaderExtensions
}

export function ShaderCode(vert: string, frag: string, extensions: ShaderExtensions): ShaderCode {
    return { id: shaderCodeId(), vert, frag, extensions }
}

export const PointsShaderCode = ShaderCode(
    require('mol-gl/shader/points.vert'),
    require('mol-gl/shader/points.frag'),
    { standardDerivatives: false, fragDepth: false }
)

export const SpheresShaderCode = ShaderCode(
    require('mol-gl/shader/spheres.vert'),
    require('mol-gl/shader/spheres.frag'),
    { standardDerivatives: false, fragDepth: true }
)

export const TextShaderCode = ShaderCode(
    require('mol-gl/shader/text.vert'),
    require('mol-gl/shader/text.frag'),
    { standardDerivatives: true, fragDepth: false }
)

export const LinesShaderCode = ShaderCode(
    require('mol-gl/shader/lines.vert'),
    require('mol-gl/shader/lines.frag'),
    { standardDerivatives: false, fragDepth: false }
)

export const MeshShaderCode = ShaderCode(
    require('mol-gl/shader/mesh.vert'),
    require('mol-gl/shader/mesh.frag'),
    { standardDerivatives: true, fragDepth: false }
)

export const GaussianDensityShaderCode = ShaderCode(
    require('mol-gl/shader/gaussian-density.vert'),
    require('mol-gl/shader/gaussian-density.frag'),
    { standardDerivatives: false, fragDepth: false }
)

export const DirectVolumeShaderCode = ShaderCode(
    require('mol-gl/shader/direct-volume.vert'),
    require('mol-gl/shader/direct-volume.frag'),
    { standardDerivatives: false, fragDepth: true }
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
        if (v !== undefined) {
            if (typeof v === 'string') {
                lines.push(`#define ${name}_${v}`)
            } else if (typeof v === 'number') {
                lines.push(`#define ${name} ${v}`)
            } else if (typeof v === 'boolean') {
                if (v) lines.push(`#define ${name}`)
            } else {
                throw new Error('unknown define type')
            }
        }
    }
    return lines.join('\n') + '\n'
}

function getGlsl100FragPrefix(ctx: WebGLContext, extensions: ShaderExtensions) {
    const prefix: string[] = []
    if (extensions.standardDerivatives) {
        prefix.push('#extension GL_OES_standard_derivatives : enable')
        prefix.push('#define enabledStandardDerivatives')
    }
    if (extensions.fragDepth) {
        if (ctx.extensions.fragDepth) {
            prefix.push('#extension GL_EXT_frag_depth : enable')
            prefix.push('#define enabledFragDepth')
        }
    }
    return prefix.join('\n') + '\n'
}

const glsl300VertPrefix = `#version 300 es
#define attribute in
#define varying out
#define texture2D texture
`

const glsl300FragPrefix = `#version 300 es
#define varying in
layout(location = 0) out highp vec4 out_FragColor;
#define gl_FragColor out_FragColor
#define gl_FragDepthEXT gl_FragDepth
#define texture2D texture

#define enabledStandardDerivatives
#define enabledFragDepth
`

export function addShaderDefines(ctx: WebGLContext, defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const { isWebGL2 } = ctx
    const header = getDefinesCode(defines)
    const vertPrefix = isWebGL2 ? glsl300VertPrefix : ''
    const fragPrefix = isWebGL2 ? glsl300FragPrefix : getGlsl100FragPrefix(ctx, shaders.extensions)
    return {
        id: shaderCodeId(),
        vert: `${vertPrefix}${header}${shaders.vert}`,
        frag: `${fragPrefix}${header}${shaders.frag}`,
        extensions: shaders.extensions
    }
}
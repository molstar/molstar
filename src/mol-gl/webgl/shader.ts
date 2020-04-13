/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { GLRenderingContext } from './compat';
import { isDebugMode } from '../../mol-util/debug';

const getNextShaderId = idFactory();

function addLineNumbers(source: string) {
    const lines = source.split('\n');
    for (let i = 0; i < lines.length; ++i) {
        lines[i] = (i + 1) + ': ' + lines[i];
    }
    return lines.join('\n');
}

export type ShaderType = 'vert' | 'frag'
export type ShaderProps = { type: ShaderType, source: string }
export interface Shader {
    readonly id: number
    attach: (program: WebGLProgram) => void
    reset: () => void
    destroy: () => void
}

function getShader(gl: GLRenderingContext, props: ShaderProps) {
    const { type, source } = props;
    const shader = gl.createShader(type === 'vert' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER);
    if (shader === null) {
        throw new Error(`Error creating ${type} shader`);
    }

    gl.shaderSource(shader, source);
    gl.compileShader(shader);

    if (isDebugMode && gl.getShaderParameter(shader, gl.COMPILE_STATUS) === false) {
        console.warn(`'${type}' shader info log '${gl.getShaderInfoLog(shader)}'\n${addLineNumbers(source)}`);
        throw new Error(`Error compiling ${type} shader`);
    }

    return shader;
}

export function createShader(gl: GLRenderingContext, props: ShaderProps): Shader {
    let shader = getShader(gl, props);

    return {
        id: getNextShaderId(),
        attach: (program: WebGLProgram) => {
            gl.attachShader(program, shader);
        },

        reset: () => {
            shader = getShader(gl, props);
        },
        destroy: () => {
            gl.deleteShader(shader);
        }
    };
}
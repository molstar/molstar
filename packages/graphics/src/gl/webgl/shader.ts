/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { GLRenderingContext } from './compat.js';
import { isDebugMode } from '@molstar/core/util/debug';

const getNextShaderId = idFactory();

function addLineNumbers(source: string) {
  const lines = source.split('\n');
  for (let i = 0; i < lines.length; ++i) {
    lines[i] = i + 1 + ': ' + lines[i];
  }
  return lines.join('\n');
}

export type ShaderType = 'vert' | 'frag';
export type ShaderProps = { type: ShaderType; source: string };
export interface Shader {
  readonly id: number;
  attach: (program: WebGLProgram) => void;
  /**
   * Throws where the compile failed, its info log and source warned: a synchronous query that waits for the compile,
   * so asked once a program using it linked, not right after `compileShader` (which loses the parallel compile)
   */
  check: () => void;
  destroy: () => void;
}

function compileShader(gl: GLRenderingContext, props: ShaderProps) {
  const { type, source } = props;
  const shader = gl.createShader(type === 'vert' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER);
  if (shader === null) {
    throw new Error(`Error creating ${type} shader`);
  }

  gl.shaderSource(shader, source);
  gl.compileShader(shader);
  return shader;
}

function checkShader(gl: GLRenderingContext, shader: WebGLShader, { type, source }: ShaderProps) {
  if (gl.getShaderParameter(shader, gl.COMPILE_STATUS) === false) {
    console.warn(`'${type}' shader info log '${gl.getShaderInfoLog(shader)}'\n${addLineNumbers(source)}`);
    throw new Error(`Error compiling ${type} shader`);
  }
}

/** A shader compiled, in debug mode checked at once: for a program used right away */
export function getShader(gl: GLRenderingContext, props: ShaderProps) {
  const shader = compileShader(gl, props);
  if (isDebugMode) checkShader(gl, shader, props);
  return shader;
}

export function createShader(gl: GLRenderingContext, props: ShaderProps): Shader {
  const shader = compileShader(gl, props);

  return {
    id: getNextShaderId(),
    attach: (program: WebGLProgram) => {
      gl.attachShader(program, shader);
    },
    check: () => checkShader(gl, shader, props),

    destroy: () => {
      gl.deleteShader(shader);
    },
  };
}

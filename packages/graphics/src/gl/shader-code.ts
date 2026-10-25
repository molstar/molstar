/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { ValueCell } from '@molstar/core/util';
import { idFactory } from '@molstar/core/util/id-factory';

export type DefineKind = 'boolean' | 'string' | 'number';
export type DefineType = boolean | string;
export type DefineValues = { [k: string]: ValueCell<DefineType> };

export type ShaderDefines = {
  [k: string]: ValueCell<DefineType>;
};

const shaderCodeId = idFactory();

type ShaderExtensionsValue = 'required' | 'optional';
export interface ShaderExtensions {
  readonly fragDepth?: ShaderExtensionsValue;
  readonly drawBuffers?: ShaderExtensionsValue;
  readonly shaderTextureLod?: ShaderExtensionsValue;
  /** Needed to enable the `gl_DrawID` built-in */
  readonly multiDraw?: ShaderExtensionsValue;
  readonly clipCullDistance?: ShaderExtensionsValue;
  readonly conservativeDepth?: ShaderExtensionsValue;
  /** Needed to enable the `gl_ViewID_OVR` built-in */
  readonly multiview2?: ShaderExtensionsValue;
}

/** Whether a define doesn't affect the program of a render variant (left out of program keys and sources) */
export type IgnoreDefine = (name: string, variant: string, defines: ShaderDefines) => boolean;

/**
 * A program, without sources: each backend looks its sources up by `name` (GLSL in `gl/webgl/programs.ts`, WGSL in
 * `gl/webgpu/shader/programs.ts`), so names are unique
 */
export interface ShaderCode {
  readonly id: number;
  readonly name: string;
  /** feature requirements, e.g. frag depth or draw buffers on WebGL */
  readonly extensions: ShaderExtensions;
  readonly ignoreDefine?: IgnoreDefine;
  /**
   * Whether the items are quads of six vertex ids without elements (`QuadIndicesPerItem`): drawn indexed by the
   * context's shared quad indices (`quad-indices.ts`)
   */
  readonly quads?: boolean;
}

export function ShaderCode(
  name: string,
  extensions: ShaderExtensions = {},
  ignoreDefine?: IgnoreDefine,
  quads?: boolean,
): ShaderCode {
  return { id: shaderCodeId(), name, extensions, ignoreDefine, quads };
}

// Note: `drawBuffers` need to be 'optional' for wboit

function ignoreDefine(name: string, variant: string, defines: ShaderDefines): boolean {
  // WebGL bounds the clip loops by the uniform `uClipObjectCount` (a define on WebGL1): one program for all counts
  if (name === 'dClipObjectCount') return true;
  if (variant.startsWith('color') || variant === 'tracing') {
    if (name === 'dLightCount') {
      return !!defines.dIgnoreLight?.ref.value;
    }
  } else {
    const ignore = ['dSubstanceType', 'dSubstance', 'dColorMarker', 'dCelShaded', 'dLightCount'];
    if (variant !== 'depth' && variant !== 'marking' && !variant.startsWith('pick')) {
      ignore.push('dXrayShaded');
    }
    if (variant !== 'emissive') {
      // emissive bloom seeds from the object's color, so keep the color/overpaint/emissive defines for it
      ignore.push('dColorType', 'dUsePalette', 'dOverpaintType', 'dOverpaint', 'dEmissiveType', 'dEmissive');
    }
    return ignore.includes(name);
  }
  return false;
}

function ignoreDefineUnlit(name: string, variant: string, defines: ShaderDefines): boolean {
  if (name === 'dLightCount') return true;
  return ignoreDefine(name, variant, defines);
}

export const PointsShaderCode = ShaderCode('points', { drawBuffers: 'optional' }, ignoreDefineUnlit, true);
export const SpheresShaderCode = ShaderCode(
  'spheres',
  { fragDepth: 'required', drawBuffers: 'optional' },
  ignoreDefine,
  true,
);
export const CylindersShaderCode = ShaderCode(
  'cylinders',
  { fragDepth: 'required', drawBuffers: 'optional' },
  ignoreDefine,
);
export const TextShaderCode = ShaderCode('text', { fragDepth: 'optional', drawBuffers: 'optional' }, ignoreDefineUnlit);
export const LinesShaderCode = ShaderCode('lines', { drawBuffers: 'optional' }, ignoreDefineUnlit);
export const MeshShaderCode = ShaderCode('mesh', { fragDepth: 'optional', drawBuffers: 'optional' }, ignoreDefine);
export const DirectVolumeShaderCode = ShaderCode(
  'direct-volume',
  { fragDepth: 'optional', drawBuffers: 'optional' },
  ignoreDefine,
);
export const ImageShaderCode = ShaderCode('image', { drawBuffers: 'optional' }, ignoreDefineUnlit);

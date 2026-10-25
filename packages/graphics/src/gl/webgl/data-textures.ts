/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  type DataFormat,
  type DataSpec,
  type RenderableSchema,
  TextureSpec,
  UniformSpec,
} from '../renderable/schema.js';
import { calculateTextureInfo } from '../renderable/util.js';

/**
 * WebGL's layout of the data arrays (`DataSpec`): ranges of shared data textures (`data-pages.ts`), the elements row by
 * row from an offset (no padding uploaded), the texture's size and the offset object uniforms of the backend
 * (`getDataTextureDimName`, `getDataTextureOffsetName`), read by generated accessors `read_<name>(i)`
 * (`getDataTextureAccessors`), as WebGPU reads its storage pages
 */

/** The size uniform of the data texture `name`, e.g. `uColorTexDim` for `tColor` */
export function getDataTextureDimName(name: string) {
  return `u${name.slice(1)}TexDim`;
}

/** The offset uniform of the data array `name` in its texture, in elements, e.g. `uColorOffset` for `tColor` */
export function getDataTextureOffsetName(name: string) {
  return `u${name.slice(1)}Offset`;
}

/** The size of a data texture of `count` elements: square-ish, at least 2 texels (1-texel textures misbehaved) */
export function getDataTextureSize(count: number) {
  const { width, height } = calculateTextureInfo(count, 1);
  return { width, height };
}

const dataSchemas = new WeakMap<RenderableSchema, RenderableSchema>();

/**
 * `schema` with its data arrays as textures (`TextureSpec.data`, nearest) and their size and offset uniforms; `schema`
 * itself without data arrays
 */
export function getDataTextureSchema(schema: RenderableSchema): RenderableSchema {
  let result = dataSchemas.get(schema);
  if (result) return result;
  let adapted: { [k: string]: RenderableSchema[string] } | undefined;
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'data') continue;
    adapted ??= { ...schema };
    adapted[k] = getDataTextureSpec(spec);
    adapted[getDataTextureDimName(k)] = UniformSpec('v2', undefined, 'object');
    adapted[getDataTextureOffsetName(k)] = UniformSpec('i', undefined, 'object');
  }
  result = adapted ?? schema;
  dataSchemas.set(schema, result);
  return result;
}

function getDataTextureSpec(spec: DataSpec) {
  const texture =
    spec.format === 'rgba32f'
      ? TextureSpec('image-float32', 'rgba', 'float', 'nearest', undefined, { usedIf: spec.usedIf })
      : TextureSpec('image-uint8', spec.format === 'r8' ? 'r' : 'rgba', 'ubyte', 'nearest', undefined, {
          usedIf: spec.usedIf,
        });
  texture.data = spec.format;
  return texture;
}

/** The names of the data textures of `schema` (see `getDataTextureSchema`) */
export function getDataTextureNames(schema: RenderableSchema) {
  return Object.keys(schema).filter((k) => {
    const spec = schema[k];
    return spec.type === 'texture' && spec.data !== undefined;
  });
}

/**
 * GLSL accessors `read_<name>(i)` of the data textures of `schema` (`vec4`, an `r8` element in `.r`), by `int` and by
 * `float` index, from the array's offset in its texture; after the declarations of the samplers and uniforms, so
 * self-contained. Integer math on WebGL2; float math on WebGL1, whose integers may be 16-bit (floats are exact up to
 * 2^24 elements)
 */
export function getDataTextureAccessors(schema: RenderableSchema, isWebGL2: boolean) {
  const lines: string[] = [];
  for (const name of getDataTextureNames(schema)) {
    const dim = getDataTextureDimName(name);
    const offset = getDataTextureOffsetName(name);
    const fn = `read_${name}`;
    if (isWebGL2) {
      lines.push(
        `vec4 ${fn}(const in int i) {`,
        `    int width = int(${dim}.x);`,
        `    int j = ${offset} + i;`,
        '    int y = j / width;',
        `    return texture2D(${name}, (vec2(float(j - y * width), float(y)) + 0.5) / ${dim});`,
        '}',
        `vec4 ${fn}(const in float i) { return ${fn}(int(i)); }`,
      );
    } else {
      lines.push(
        `vec4 ${fn}(const in float i) {`,
        `    float j = float(${offset}) + floor(i + 0.5);`,
        `    float y = floor((j + 0.5) / ${dim}.x);`,
        `    return texture2D(${name}, (vec2(j - y * ${dim}.x, y) + 0.5) / ${dim});`,
        '}',
        `vec4 ${fn}(const in int i) { return ${fn}(float(i)); }`,
      );
    }
  }
  return lines.join('\n');
}

/** Components per element of `format` as a texel (`r8` one, the others four) */
export function getDataTextureItemSize(format: DataFormat) {
  return format === 'r8' ? 1 : 4;
}

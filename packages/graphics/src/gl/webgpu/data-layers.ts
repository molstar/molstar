/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StorageBufferSpec } from '../compute.js';
import { type RenderableSchema, UniformSpec } from '../renderable/schema.js';

/**
 * The theme layers (`DataSpec`) WebGPU reads from shared storage pages: data by index
 * (group, instance or vertex), WebGL lays them out in data textures. `rgba8` one `u32` per element, `r8` four elements
 * per `u32`, `rgba32f` four `u32` (the floats' bits), the formats of their specs; the arrays are the layers'
 * `DataArray`s, element `i` at the same index. The positions and groups of spheres and points (`tPositionGroup`, one
 * element per item), the cylinders' data (`tCylinder`, three per cylinder), the lines' (`tLine`, two per line) and the
 * text's (`tText`, four per item) are layers too (`rgb32f-i32`: the groups integers, `readGroup_<layer>`), and so is
 * merged items' segment data (`tSegment`, `SegmentTexels` per segment). Spheres'
 * quantized centers (`tSite`, `rgba16i`: int16 x, y, z and a uint16, two `u32` an element) are one too, and their second
 * color (`tSiteColor`, a table by site type or an element a site). Sixteen layers: the page define's 32 bits
 */
export const DataLayers = {
  tColor: 'rgba8',
  tSize: 'rgba8',
  tMarker: 'r8',
  tOverpaint: 'rgba8',
  tTransparency: 'r8',
  tEmissive: 'r8',
  tSubstance: 'rgba8',
  tClipping: 'r8',
  tWiggle: 'r8',
  tPositionGroup: 'rgb32f-i32',
  tCylinder: 'rgb32f-i32',
  tSegment: 'rgba32f',
  tLine: 'rgb32f-i32',
  tText: 'rgb32f-i32',
  tSite: 'rgba16i',
  tSiteColor: 'rgba8',
} as const;
export type DataLayer = keyof typeof DataLayers;
/** Bytes per element of the encodings of `DataLayers` */
export const DataLayerElementBytes = { rgba8: 4, r8: 1, rgba32f: 16, 'rgb32f-i32': 16, rgba16i: 8 } as const;
export const DataLayerNames = Object.keys(DataLayers) as DataLayer[];

export function isDataLayer(name: string): name is DataLayer {
  return name in DataLayers;
}

/** The storage buffer of the pages, in the object's bind group */
export const DataPageName = 'sDataPage';

/**
 * The storage buffers an item's layers are read from: one page as a rule, more where its layers together exceed the
 * binding limit (`WebGPUContext.dataBindingSize`), each within it; at most 4 (of the 8 storage buffers a stage may bind)
 */
export const DataPageNames = [DataPageName, 'sDataPage1', 'sDataPage2', 'sDataPage3'] as const;
export const MaxDataPages = DataPageNames.length;

/**
 * The define of the page each layer is in (`DataPageNames`), two bits per layer in `DataLayerNames` order: 0 (all in
 * the first) for items within the limit, whose programs are those of before; not in the schema: in an item's defines
 * only where not 0
 */
export const DataLayerPagesName = 'dDataLayerPages';

/** The page of `layer` in the layout `pages` (`DataLayerPagesName`) */
export function getDataLayerPage(pages: number, layer: DataLayer) {
  return (pages >>> (2 * DataLayerNames.indexOf(layer))) & 3;
}

/** The layout of `pageOf` (a page per layer of `layers`), see `DataLayerPagesName` */
export function getDataLayerPages(layers: readonly DataLayer[], pageOf: readonly number[]) {
  let pages = 0;
  for (let i = 0; i < layers.length; ++i) pages |= pageOf[i] << (2 * DataLayerNames.indexOf(layers[i]));
  return pages;
}

/** The object uniform holding the word offset of `layer` in the page, e.g. `uColorOffset` */
export function getDataOffsetName(layer: DataLayer) {
  return `u${layer.slice(1)}Offset`;
}

const dataSchemas = new WeakMap<RenderableSchema, RenderableSchema>();

/**
 * `schema` with what the WebGPU backend reads its data layers through: per layer an object uniform with its offset
 * (`getDataOffsetName`) and the page (`DataPageName`); `schema` itself where it has no layers
 */
export function getDataLayerSchema(schema: RenderableSchema): RenderableSchema {
  let extended = dataSchemas.get(schema);
  if (extended) return extended;
  const layers = DataLayerNames.filter((n) => schema[n]?.type === 'data');
  if (layers.length === 0) {
    extended = schema;
  } else {
    const additions: { [k: string]: RenderableSchema[string] } = {};
    for (const layer of layers) additions[getDataOffsetName(layer)] = UniformSpec('i', undefined, 'object');
    for (const name of DataPageNames) additions[name] = StorageBufferSpec('u32', 'read');
    extended = { ...schema, ...additions };
  }
  dataSchemas.set(schema, extended);
  return extended;
}

/** The layers of `schema` read from the pages (its offsets, see `getDataLayerSchema`) */
function getSchemaLayers(schema: RenderableSchema) {
  return DataLayerNames.filter((n) => schema[getDataOffsetName(n)] !== undefined);
}

/**
 * WGSL accessors `read_<layer>(i: i32) -> vec4<f32>` of the layers of `schema` that `isReferenced` names, reading as
 * `textureLoad` reads the layer's texture, each from its page of the layout `pages` (`DataLayerPagesName`); programs of
 * other schemas may have functions of these names of their own
 */
export function getDataLayerAccessors(schema: RenderableSchema, isReferenced: (name: string) => boolean, pages = 0) {
  const lines: string[] = [];
  for (const layer of getSchemaLayers(schema)) {
    const fn = `read_${layer}`;
    const DataPageName = DataPageNames[getDataLayerPage(pages, layer)];
    const offset = `u32(uObject.${getDataOffsetName(layer)})`;
    if (DataLayers[layer] === 'rgb32f-i32' && isReferenced(`readGroup_${layer}`)) {
      // the group an integer: its word never read as a float (integer bits a denormal, which GPUs flush)
      lines.push(`fn readGroup_${layer}(i: i32) -> i32 {`);
      lines.push(`    return bitcast<i32>(${DataPageName}[${offset} + u32(i) * 4u + 3u]);`);
      lines.push('}');
    }
    if (!isReferenced(fn)) continue;
    if (DataLayers[layer] === 'rgba8') {
      lines.push(`fn ${fn}(i: i32) -> vec4<f32> {`);
      lines.push(`    return unpack4x8unorm(${DataPageName}[${offset} + u32(i)]);`);
    } else if (DataLayers[layer] === 'rgba32f' || DataLayers[layer] === 'rgb32f-i32') {
      // `rgb32f-i32`: the group's word too, as float bits (`readGroup_<layer>` the integer)
      lines.push(`fn ${fn}(i: i32) -> vec4<f32> {`);
      lines.push(`    let o = ${offset} + u32(i) * 4u;`);
      lines.push(
        `    return bitcast<vec4<f32>>(vec4<u32>(${DataPageName}[o], ${DataPageName}[o + 1u], ${DataPageName}[o + 2u], ${DataPageName}[o + 3u]));`,
      );
    } else if (DataLayers[layer] === 'rgba16i') {
      // little-endian int16 x, y in the first word, z and the uint16 in the second: the signed ones sign-extended
      lines.push(`fn ${fn}(i: i32) -> vec4<f32> {`);
      lines.push(`    let o = ${offset} + u32(i) * 2u;`);
      lines.push(`    let a = ${DataPageName}[o];`);
      lines.push(`    let b = ${DataPageName}[o + 1u];`);
      lines.push(
        '    return vec4<f32>(f32(bitcast<i32>(a << 16u) >> 16u), f32(bitcast<i32>(a) >> 16u), f32(bitcast<i32>(b << 16u) >> 16u), f32(b >> 16u));',
      );
    } else {
      lines.push(`fn ${fn}(i: i32) -> vec4<f32> {`);
      lines.push(`    let word = ${DataPageName}[${offset} + (u32(i) >> 2u)];`);
      lines.push(`    let value = (word >> ((u32(i) & 3u) * 8u)) & 0xffu;`);
      lines.push(`    return vec4<f32>(f32(value) / 255.0, 0.0, 0.0, 1.0);`);
    }
    lines.push('}');
  }
  return lines.join('\n');
}

/** Whether `isReferenced` names an accessor of a layer of `schema`: the page is bound then */
export function usesDataLayers(schema: RenderableSchema, isReferenced: (name: string) => boolean) {
  return getSchemaLayers(schema).some((n) => isReferenced(`read_${n}`) || isReferenced(`readGroup_${n}`));
}

/** The pages (`DataPageNames`) the accessors `isReferenced` names read in the layout `pages`: those bound */
export function getReadDataPages(schema: RenderableSchema, isReferenced: (name: string) => boolean, pages = 0) {
  const read = new Set<string>();
  for (const layer of getSchemaLayers(schema)) {
    if (isReferenced(`read_${layer}`) || isReferenced(`readGroup_${layer}`)) {
      read.add(DataPageNames[getDataLayerPage(pages, layer)]);
    }
  }
  return read;
}

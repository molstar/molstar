/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getAttributeLayout } from '../attribute-layout.js';
import { CulledData } from '../culled-data.js';
import { CopyLayerBlock, CopyLayerKinds } from '../instance-expansion.js';
import type { DefineLookup, RenderableSchema } from '../renderable/schema.js';
import { getWgslVertexAttributes, getWgslVertexInputs } from './wgsl-declarations.js';

/**
 * GPU-culled programs of merged render objects (`dCulled`, `CulledSchema` in `merged.ts`): drawn from indirect records
 * whose first instances index regions of instance indices (`sCulledIndices`), the instanced attributes read from storage
 * (`sCulledData`, `CulledData`) by the index at `instance_index`, which they can't follow as vertex inputs. The program
 * sources stay as they are: the stage input leaves the instanced attributes out, a generated reader fills `VertexInput`
 * and a generated `vertexMain` calls the program's own, renamed. Objects of poses drawn without GPU culling read theirs
 * the same way (`dInstanceStorage`, `instance-storage.ts`): from `sInstanceData` at `instance_index` itself
 */
export const CulledDefine = 'dCulled';
/**
 * a page store's instances (`paged-cull.ts`): `sCulledData` a page store's `CulledData` (its segments and its composed
 * poses, `CulledPagedFlag`), its pose read by the store index at `instance_index` of `sCulledIndices`, as a merged
 * culler's
 */
export const CulledPagedDefine = 'dCulledPaged';
const IndicesName = 'sCulledIndices';
const DataName = 'sCulledData';
const InstanceDataName = 'sInstanceData';

/**
 * What the generated reader reads: the culled storage buffers, or the instance data (referenced though the program's
 * source doesn't name them)
 */
export const CulledStorageNames: ReadonlySet<string> = new Set([IndicesName, DataName]);
export const InstanceStorageNames: ReadonlySet<string> = new Set([InstanceDataName]);

/** The instanced vertex inputs read from `sCulledData`, by name: their reader of instance `i` */
const CulledInputs: { readonly [name: string]: string } = {
  aTransform0: 'transformRows[0]',
  aTransform1: 'transformRows[1]',
  aTransform2: 'transformRows[2]',
  aSegment: 'culledSegment(i)',
  aInstance: 'culledInstance(i)',
};

/**
 * The vertex input declarations of a culled program of `schema`: `VertexStageInput` (the attributes used with `defines`,
 * per vertex), `VertexInput` (those and the instanced ones as the program reads them) and `readCulledVertexInput`;
 * `indexed`: through the culler's instance indices, else the instance data at `instance_index` (`dInstanceStorage`);
 * `paged`: a page store's (`dCulledPaged`): read alike, an instance's index its store index, its pose's slot
 */
export function getCulledVertexInputDeclaration(
  schema: RenderableSchema,
  defines: DefineLookup,
  indexed = true,
  paged = false,
) {
  const stageEntries = getWgslVertexAttributes(schema, defines);
  const instanced = getAttributeLayout(schema).entries.filter((e) => e.divisor !== 0 && !stageEntries.includes(e));
  const stage = [
    'struct VertexStageInput {',
    '    @builtin(vertex_index) vertexIndex: u32,',
    '    @builtin(instance_index) instanceIndex: u32,',
  ];
  const input = ['struct VertexInput {', '    vertexIndex: u32,', '    instanceIndex: u32,'];
  const reader = [
    'fn readCulledVertexInput(stage: VertexStageInput) -> VertexInput {',
    '    var in: VertexInput;',
    '    in.vertexIndex = stage.vertexIndex;',
    indexed ? `    let i = ${IndicesName}[stage.instanceIndex];` : '    let i = stage.instanceIndex;',
    '    in.instanceIndex = i;',
    // the transform's rows, built from the pose where the instance has one (`CulledData`; a page store's its store
    // index's)
    '    let transformRows = culledRows(i);',
  ];
  for (const entry of stageEntries) {
    if (entry.divisor !== 0)
      throw new Error(`instanced attribute '${entry.name}' of a culled program is a vertex input`);
    for (const { name, location, type } of getWgslVertexInputs(entry)) {
      stage.push(`    @location(${location}) ${name}: ${type},`);
      input.push(`    ${name}: ${type},`);
      reader.push(`    in.${name} = stage.${name};`);
    }
  }
  for (const entry of instanced) {
    for (const { name, type } of getWgslVertexInputs(entry)) {
      const read = CulledInputs[name];
      if (!read) throw new Error(`instanced attribute '${name}' has no culled data to be read from`);
      input.push(`    ${name}: ${type},`);
      reader.push(`    in.${name} = ${read};`);
    }
  }
  stage.push('};');
  input.push('};');
  reader.push('    return in;', '}');
  const data = CulledData.wgsl(indexed ? DataName : InstanceDataName);
  // added after the preprocessor: by the defines, not `#if`
  const copyLayers = defines('dCopyLayers') ? copy_layers_wgsl : '';
  return [stage.join('\n'), input.join('\n'), data, copyLayers, reader.join('\n')].join('\n\n');
}

/**
 * Layer `k` of a page store's copy's instance (`sCopyLayers`: two words a copy's kind, `CopyLayerKinds` a copy: a
 * value's flag and its first byte, or a block's and its first byte less its template's first store index times the
 * size), `tag` its copy plus one, `index` its template particle's store index; else `own`, its store index's
 */
const copy_layers_wgsl = /* wgsl */ `
fn copyLayerOr(k: u32, size: u32, tag: u32, index: u32, own: vec4<f32>) -> vec4<f32> {
    if (tag == 0u) { return own; }
    let o = 2u * (${CopyLayerKinds}u * (tag - 1u) + k);
    let flag = sCopyLayers[o];
    if (flag == 0u) { return own; }
    var b = sCopyLayers[o + 1u];
    if ((flag & ${CopyLayerBlock}u) != 0u) { b = u32(bitcast<i32>(b) + i32(index * size)); }
    let w = sCopyLayers[b / 4u];
    if (size == 4u) { return unpack4x8unorm(w); }
    return vec4<f32>(f32((w >> (8u * (b % 4u))) & 0xffu) / 255.0);
}
`;

const VertexMain = /@vertex\s+fn vertexMain\(in: VertexInput\) -> VertexOutput \{/;

/** The program's `vertexMain` renamed to `vertexBody`, called by a generated one reading the culled instance's data */
export function wrapCulledVertexMain(body: string) {
  if (!VertexMain.test(body)) throw new Error('culled programs need `fn vertexMain(in: VertexInput) -> VertexOutput`');
  const renamed = body.replace(VertexMain, 'fn vertexBody(in: VertexInput) -> VertexOutput {');
  return `${renamed}

@vertex
fn vertexMain(stage: VertexStageInput) -> VertexOutput {
    return vertexBody(readCulledVertexInput(stage));
}`;
}

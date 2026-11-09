/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getAttributeLayout } from '../attribute-layout.js';
import { CulledData } from '../culled-data.js';
import type { DefineLookup, RenderableSchema } from '../renderable/schema.js';
import { getWgslVertexAttributes, getWgslVertexInputs } from './wgsl-declarations.js';

/**
 * GPU-culled programs of merged render objects (`dCulled`, `CulledSchema` in `merged.ts`): drawn from indirect records
 * whose first instances index regions of instance indices (`sCulledIndices`), the instanced attributes read from storage
 * (`sCulledData`, `CulledData`) by the index at `instance_index`, which they can't follow as vertex inputs. The program
 * sources stay as they are: the stage input leaves the instanced attributes out, a generated reader fills `VertexInput`
 * and a generated `vertexMain` calls the program's own, renamed
 */
export const CulledDefine = 'dCulled';
const IndicesName = 'sCulledIndices';
const DataName = 'sCulledData';

/** What the generated reader reads: the culled storage buffers (referenced though the program's source doesn't name them) */
export const CulledStorageNames: ReadonlySet<string> = new Set([IndicesName, DataName]);

/** The instanced vertex inputs read from `sCulledData`, by name: their reader of instance `i` */
const CulledInputs: { readonly [name: string]: string } = {
  aTransform0: 'culledTransform0(i)',
  aTransform1: 'culledTransform1(i)',
  aTransform2: 'culledTransform2(i)',
  aSegment: 'culledSegment(i)',
  aInstance: 'culledInstance(i)',
};

/**
 * The vertex input declarations of a culled program of `schema`: `VertexStageInput` (the attributes used with `defines`,
 * per vertex), `VertexInput` (those and the instanced ones as the program reads them) and `readCulledVertexInput`
 */
export function getCulledVertexInputDeclaration(schema: RenderableSchema, defines: DefineLookup) {
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
    `    let i = ${IndicesName}[stage.instanceIndex];`,
    '    in.instanceIndex = i;',
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
  return [
    stage.join('\n'),
    input.join('\n'),
    CulledData.wgsl(DataName, 'u32(uObject.uCulledInstanceCount)'),
    reader.join('\n'),
  ].join('\n\n');
}

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

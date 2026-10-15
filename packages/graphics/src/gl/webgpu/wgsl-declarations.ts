/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type AttributeLayoutEntry, getAttributeLayout } from '../attribute-layout.js';
import { type DefineLookup, isUsed, type RenderableSchema, type UniformBlock } from '../renderable/schema.js';
import { type BindingLayout, type TextureBinding, UniformBlockInfo, UniformBlocks } from '../uniform-layout.js';
import { getWgslUniformDeclarations } from './uniform-declarations.js';

/** Shared immutable samplers of group 0, all clamp-to-edge */
export const WgslSamplers = ['sNearest', 'sLinear', 'sLinearMipmap'] as const;
export type WgslSampler = (typeof WgslSamplers)[number];

/** Binding of the first shared sampler in group 0, after the global block and textures */
export function getWgslSamplerBinding(layout: BindingLayout) {
  return layout.global.textures.length + 1;
}

export function getWgslTextureType(binding: Pick<TextureBinding, 'dimension' | 'sampleType'>) {
  const { dimension, sampleType } = binding;
  if (sampleType === 'depth') {
    if (dimension === '3d') throw new Error('3d depth textures are not supported');
    return dimension === 'cube' ? 'texture_depth_cube' : 'texture_depth_2d';
  }
  const type = sampleType === 'sint' ? 'i32' : sampleType === 'uint' ? 'u32' : 'f32';
  return `texture_${dimension === 'cube' ? 'cube' : dimension}<${type}>`;
}

/**
 * Texture and sampler declarations of `layout`; with `isReferenced`, only textures it accepts are
 * declared (the binding numbers stay those of the layout)
 */
export function getWgslTextureDeclarations(layout: BindingLayout, isReferenced?: (name: string) => boolean) {
  const lines: string[] = [];
  for (const block of UniformBlocks) {
    const group = UniformBlockInfo[block].bindingIndex;
    for (const t of layout[block].textures) {
      if (isReferenced && !isReferenced(t.name)) continue;
      lines.push(`@group(${group}) @binding(${t.binding}) var ${t.name}: ${getWgslTextureType(t)};`);
    }
  }
  const samplerBinding = getWgslSamplerBinding(layout);
  for (let i = 0; i < WgslSamplers.length; ++i) {
    const name = WgslSamplers[i];
    if (isReferenced && !isReferenced(name)) continue;
    lines.push(`@group(0) @binding(${samplerBinding + i}) var ${name}: sampler;`);
  }
  return lines.join('\n');
}

function getWgslAttributeType(itemSize: number) {
  switch (itemSize) {
    case 1:
      return 'f32';
    case 2:
    case 3:
    case 4:
      return `vec${itemSize}<f32>`;
    default:
      throw new Error(`unsupported attribute item size '${itemSize}'`);
  }
}

/** Vertex inputs of an attribute, `mat4` attributes are split into four `vec4<f32>` inputs (`aTransform0` to `aTransform3`) */
export function getWgslVertexInputs(entry: AttributeLayoutEntry): { name: string; location: number; type: string }[] {
  if (entry.locationCount === 4) {
    return [0, 1, 2, 3].map((i) => ({ name: `${entry.name}${i}`, location: entry.location + i, type: 'vec4<f32>' }));
  }
  return [{ name: entry.name, location: entry.location, type: getWgslAttributeType(entry.itemSize) }];
}

/** Attributes of `schema` that are vertex inputs, in location order */
export function getWgslVertexAttributes(schema: RenderableSchema, defines: DefineLookup) {
  return getAttributeLayout(schema).entries.filter((e) =>
    isUsed(schema[e.name] as { usedIf?: readonly string[] }, defines),
  );
}

/** `VertexInput` struct with the attribute locations of `getAttributeLayout` */
export function getWgslVertexInputDeclaration(schema: RenderableSchema, defines: DefineLookup) {
  const lines = ['struct VertexInput {'];
  lines.push('    @builtin(vertex_index) vertexIndex: u32,');
  lines.push('    @builtin(instance_index) instanceIndex: u32,');
  for (const entry of getWgslVertexAttributes(schema, defines)) {
    for (const { name, location, type } of getWgslVertexInputs(entry)) {
      lines.push(`    @location(${location}) ${name}: ${type},`);
    }
  }
  lines.push('};');
  return lines.join('\n');
}

export type WgslDeclarations = {
  readonly code: string;
  /** uniforms WebGPU can't declare (classic uniforms without a block) */
  readonly unsupported: readonly string[];
};

/**
 * All generated declarations of a program: uniform blocks, `Draw`, textures, samplers and vertex
 * inputs; with `isReferenced`, only referenced textures and samplers are declared
 */
export function getWgslDeclarations(
  schema: RenderableSchema,
  layout: BindingLayout,
  defines: DefineLookup,
  isReferenced?: (name: string) => boolean,
): WgslDeclarations {
  const uniforms = getWgslUniformDeclarations(schema);
  const code = [
    uniforms.code,
    getWgslTextureDeclarations(layout, isReferenced),
    getWgslVertexInputDeclaration(schema, defines),
  ].filter((c) => c.length > 0);
  return { code: code.join('\n\n'), unsupported: uniforms.unsupported };
}

/** Bind group index of a tier */
export function getWgslGroup(block: UniformBlock) {
  return UniformBlockInfo[block].bindingIndex;
}

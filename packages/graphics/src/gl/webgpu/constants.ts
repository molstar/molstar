/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// the flag namespaces (`GPUBufferUsage` etc.) exist at runtime but are not declared in lib.dom

export const BufferUsage = {
  MAP_READ: 0x1,
  MAP_WRITE: 0x2,
  COPY_SRC: 0x4,
  COPY_DST: 0x8,
  INDEX: 0x10,
  VERTEX: 0x20,
  UNIFORM: 0x40,
  STORAGE: 0x80,
  INDIRECT: 0x100,
  QUERY_RESOLVE: 0x200,
} as const;

export const TextureUsage = {
  COPY_SRC: 0x1,
  COPY_DST: 0x2,
  TEXTURE_BINDING: 0x4,
  STORAGE_BINDING: 0x8,
  RENDER_ATTACHMENT: 0x10,
} as const;

export const ShaderStage = {
  VERTEX: 0x1,
  FRAGMENT: 0x2,
  COMPUTE: 0x4,
} as const;

export const MapMode = {
  READ: 0x1,
  WRITE: 0x2,
} as const;

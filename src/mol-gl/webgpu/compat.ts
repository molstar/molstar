/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

// WebGPU bit masks from https://www.w3.org/TR/webgpu/.
// TypeScript 6's DOM library includes the interfaces but omits these namespaces.
export const GPUBufferUsage = { MAP_READ: 1, MAP_WRITE: 2, COPY_SRC: 4, COPY_DST: 8, INDEX: 16, VERTEX: 32, UNIFORM: 64, STORAGE: 128, INDIRECT: 256, QUERY_RESOLVE: 512 } as const;
export const GPUTextureUsage = { COPY_SRC: 1, COPY_DST: 2, TEXTURE_BINDING: 4, STORAGE_BINDING: 8, RENDER_ATTACHMENT: 16 } as const;
export const GPUShaderStage = { VERTEX: 1, FRAGMENT: 2, COMPUTE: 4 } as const;
export const GPUMapMode = { READ: 1, WRITE: 2 } as const;

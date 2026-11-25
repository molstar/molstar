/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

// Workgroup counters reset by the kernel itself: WGSL zero-initializes workgroup memory, but a device may skip it where
// the driver claims to (SwiftShader's WebGPU adapter kept the previous dispatch's values: the GPU isosurface's triangles
// written past their limit, most of the surface missing). A counter only added to must be stored first.

const shaderDir = path.join(path.dirname(fileURLToPath(import.meta.url)), '../webgpu/shader');

describe('WGSL workgroup counters', () => {
  it('are reset by their kernel (atomicStore), not left to the zero-initialization of workgroup memory', () => {
    const missing: string[] = [];
    for (const file of fs.readdirSync(shaderDir)) {
      if (!file.endsWith('.ts')) continue;
      const source = fs.readFileSync(path.join(shaderDir, file), 'utf8');
      for (const m of source.matchAll(/var<workgroup>\s+(\w+)\s*:\s*atomic</g)) {
        const name = m[1];
        if (!new RegExp(`atomicStore\\(\\s*&${name}\\s*,`).test(source)) missing.push(`${file}: ${name}`);
      }
    }
    expect(missing).toEqual([]);
  });
});

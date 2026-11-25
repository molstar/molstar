/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const GlSource = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

/** Relative imports kept at runtime (not `import type`, nor `export type`), as source paths */
const RuntimeImport = /^(?:import|export)\s+(?!type\b)[^'";]*?from\s+'(\.[^']+)';/gm;

/** The modules `file` reaches at runtime through relative imports, itself included, as paths below `gl/` */
function getRuntimeModules(file: string, seen = new Set<string>()) {
  const rel = path.relative(GlSource, file).split(path.sep).join('/');
  if (seen.has(rel)) return seen;
  seen.add(rel);
  const text = fs.readFileSync(file, 'utf8');
  for (const [, specifier] of text.matchAll(RuntimeImport)) {
    getRuntimeModules(path.resolve(path.dirname(file), specifier.replace(/\.js$/, '.ts')), seen);
  }
  return seen;
}

// An extension registers its programs and kernels without the package's own sources: a WebGL-only build registering an
// extension's WGSL (or a WebGPU-only one its GLSL) bundles nothing of the other backend but its registry
describe('the extension modules', () => {
  it('reach only the registries at runtime, no built-in sources nor contexts', () => {
    expect([...getRuntimeModules(path.join(GlSource, 'extend/webgl.ts'))].sort()).toEqual([
      'extend/webgl.ts',
      'webgl/program-registry.ts',
    ]);
    expect([...getRuntimeModules(path.join(GlSource, 'extend/webgpu.ts'))].sort()).toEqual([
      'extend/webgpu.ts',
      'webgpu/shader/kernel-registry.ts',
      'webgpu/shader/program-registry.ts',
    ]);
  });

  it('finds what a module with built-in sources reaches', () => {
    // the check itself: the built-in table imports its sources
    const reached = getRuntimeModules(path.join(GlSource, 'webgpu/shader/programs.ts'));
    expect(reached.has('webgpu/shader/mesh.wgsl.ts')).toBe(true);
    expect(reached.has('webgpu/shader/program-registry.ts')).toBe(true);
  });
});

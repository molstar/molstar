/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import * as fs from 'node:fs';
import * as path from 'node:path';

/** WebGL feature checks and parameters; use `capabilities` and `resources.limits` of the `GraphicsContext` instead */
const MemberPatterns = [
  /\.extensions\b/,
  /\bisWebGL2\b/,
  /(?<!limits)\.max(TextureSize|3dTextureSize|RenderbufferSize|DrawBuffers|TextureImageUnits)\b/,
];
/** destructured WebGL extensions, only checked in graphics where they are not false positives */
const DestructurePatterns = [/\bextensions\s*:\s*\{/, /\{[^}]*\bextensions\b[^}]*\}\s*=/];

/** Files still using WebGL specifics outside `gl/webgl/` */
const Allowlist: { [file: string]: string } = {};

function getSourceFiles(dir: string, files: string[] = []) {
  if (!fs.existsSync(dir)) return files;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name !== '_test' && entry.name !== 'webgl' && entry.name !== 'webgpu') getSourceFiles(p, files);
    } else if (/\.tsx?$/.test(entry.name)) {
      files.push(p.split(path.sep).join('/'));
    }
  }
  return files;
}

function getSourceRoots(parents: string[]) {
  const roots = ['packages/graphics/src'];
  for (const parent of parents) {
    for (const entry of fs.readdirSync(parent, { withFileTypes: true })) {
      if (entry.isDirectory()) roots.push(`${parent}/${entry.name}/src`);
    }
  }
  return roots;
}

describe('backend-neutral feature checks', () => {
  const matches = new Map<string, string[]>();
  for (const root of getSourceRoots(['packages/plugin', 'extensions'])) {
    const patterns = root === 'packages/graphics/src' ? [...MemberPatterns, ...DestructurePatterns] : MemberPatterns;
    for (const file of getSourceFiles(root)) {
      const lines = fs.readFileSync(file, 'utf8').split('\n');
      lines.forEach((line, i) => {
        if (!patterns.some((p) => p.test(line))) return;
        if (!matches.has(file)) matches.set(file, []);
        matches.get(file)!.push(`${file}:${i + 1}: ${line.trim()}`);
      });
    }
  }

  it('no WebGL extensions, isWebGL2 or parameters outside gl/webgl/', () => {
    const violations = [...matches].filter(([file]) => !Allowlist[file]).flatMap(([, lines]) => lines);
    expect(violations).toEqual([]);
  });

  it('allowlist has no stale entries', () => {
    expect(Object.keys(Allowlist).filter((file) => !matches.has(file))).toEqual([]);
  });
});

/** Import specifiers: `from '…'` (also of `export … from`), side-effect `import '…'` and dynamic `import('…')` */
const ImportPattern = /(?:\bfrom\s*|\bimport\s*\(?\s*)'([^']+)'/g;

const GraphicsSource = 'packages/graphics/src';

/**
 * The backend and module (`webgl/context`, also for nested paths, e.g. GLSL sources in `webgl/shader/`) an import of
 * `file` resolves to, by path: relative imports from the file's folder, `@molstar/graphics/…` from the graphics sources
 */
function getBackendModule(file: string, specifier: string): string | undefined {
  let resolved: string;
  if (specifier.startsWith('@molstar/graphics/')) {
    resolved = `${GraphicsSource}/${specifier.slice('@molstar/graphics/'.length)}`;
  } else if (specifier.startsWith('.')) {
    resolved = path.posix.join(path.posix.dirname(file), specifier);
  } else {
    return;
  }
  const m = resolved.match(/^packages\/graphics\/src\/gl\/(webgl|webgpu)\/([\w-]+)/);
  return m ? `${m[1]}/${m[2].replace(/\.js$/, '')}` : undefined;
}

/** Backend modules that neutral code may import */
const PendingModules: { [module: string]: string } = {
  'webgl/compat': 'isWebGL2, see the feature check allowlist',
  'webgl/backend': 'backend descriptors, consumers select the backends they ship',
  'webgpu/backend': 'backend descriptors, consumers select the backends they ship',
};

/** Other backend imports outside the backends: creating contexts */
const ImportAllowlist: { [file: string]: string } = {
  'packages/graphics/src/canvas3d/webgl.ts': 'createContext, getGLContext: `createWebGLCanvas3DContext` (WebGL)',
  'packages/plugin/headless/src/screenshot.ts': 'the contexts of the headless gl package and of Dawn (Node.js)',
};

describe('backend import resolution', () => {
  it('resolves relative and package imports by path', () => {
    const file = `${GraphicsSource}/canvas3d/webgl.ts`;
    expect(getBackendModule(file, '../gl/webgl/context.js')).toBe('webgl/context');
    expect(getBackendModule(file, '@molstar/graphics/gl/webgpu/backend')).toBe('webgpu/backend');
    expect(getBackendModule(`${GraphicsSource}/gl/backends/webgl.ts`, '../webgl/shader/mesh.vert.js')).toBe(
      'webgl/shader',
    );
    expect(getBackendModule(file, './canvas3d.js')).toBeUndefined();
    expect(getBackendModule(file, '../gl/context.js')).toBeUndefined();
    expect(getBackendModule(file, 'webgl-utils')).toBeUndefined();
    expect(
      [..."import { a } from '../gl/webgl/context.js'; import '../gl/webgpu/x.js';".matchAll(ImportPattern)].map(
        (m) => m[1],
      ),
    ).toEqual(['../gl/webgl/context.js', '../gl/webgpu/x.js']);
  });
});

describe('backend-neutral imports', () => {
  const imports = new Map<string, string[]>();
  for (const root of getSourceRoots(['packages/plugin', 'extensions', 'apps'])) {
    for (const file of getSourceFiles(root)) {
      const lines = fs.readFileSync(file, 'utf8').split('\n');
      lines.forEach((line, i) => {
        for (const [, specifier] of line.matchAll(ImportPattern)) {
          const module = getBackendModule(file, specifier);
          if (!module || PendingModules[module]) continue;
          if (!imports.has(file)) imports.set(file, []);
          imports.get(file)!.push(`${file}:${i + 1}: ${line.trim()}`);
        }
      });
    }
  }

  it('no gl/webgl/ or gl/webgpu/ imports outside the backends', () => {
    const violations = [...imports].filter(([file]) => !ImportAllowlist[file]).flatMap(([, lines]) => lines);
    expect(violations).toEqual([]);
  });

  it('import allowlist has no stale entries', () => {
    expect(Object.keys(ImportAllowlist).filter((file) => !imports.has(file))).toEqual([]);
  });
});

/**
 * Narrowing to a backend; decide by `capabilities` or services (e.g. `graphics.compute`, registered implementations),
 * or move the code into the backend
 */
const NarrowingPatterns = [/\bis(WebGL|WebGPU)Context\(/, /\bas (WebGL|WebGPU)Context\b/];

describe('backend narrowing', () => {
  it('no narrowing to a backend outside gl/webgl/ and gl/webgpu/', () => {
    const violations: string[] = [];
    for (const root of getSourceRoots(['packages/plugin', 'extensions', 'apps'])) {
      for (const file of getSourceFiles(root)) {
        fs.readFileSync(file, 'utf8')
          .split('\n')
          .forEach((line, i) => {
            if (NarrowingPatterns.some((p) => p.test(line))) violations.push(`${file}:${i + 1}: ${line.trim()}`);
          });
      }
    }
    expect(violations).toEqual([]);
  });
});

/** All `.ts` files below `dir`, including tests */
function getBackendFiles(dir: string, files: string[] = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) getBackendFiles(p, files);
    else if (entry.name.endsWith('.ts')) files.push(p.split(path.sep).join('/'));
  }
  return files;
}

describe('backend isolation', () => {
  const backends = [
    { dir: 'packages/graphics/src/gl/webgl', other: 'webgpu' },
    { dir: 'packages/graphics/src/gl/webgpu', other: 'webgl' },
  ];
  for (const { dir, other } of backends) {
    it(`no ${other} imports in ${dir}`, () => {
      const pattern = new RegExp(`from '[^']*/${other}/`);
      const violations: string[] = [];
      for (const file of getBackendFiles(dir)) {
        fs.readFileSync(file, 'utf8')
          .split('\n')
          .forEach((line, i) => {
            if (pattern.test(line)) violations.push(`${file}:${i + 1}: ${line.trim()}`);
          });
      }
      expect(violations).toEqual([]);
    });
  }
});

/** Fixed-function state, clears, draw buffers and copies; use `setRenderState`, `clear`, `setDrawBuffers` and `copyFromFramebuffer` */
const RenderStatePatterns = [
  /\bstate\.(enable|disable|blend\w*|colorMask|depthMask|depthFunc|clearColor|clearDepth|clearStencil|cullFace|frontFace|stencil\w*|viewport|scissor)\(/,
  /\bgl\.(clear\w*|copyTexSubImage2D|drawBuffers|enable|disable|blend\w*|colorMask|depthMask|viewport|scissor)\(/,
  /\.drawBuffers\(\[/,
  // binding caches are maintained by the backend (render items, textures, global bindings)
  /\bstate\.current(Program|Material|RenderItem)Id\b/,
  // downcasts to WebGL implementations, use the `asGLTexture`/`asWebGLGlobalBindings` guards
  /\bas (GLTexture|WebGLGlobalBindings)\b/,
];

/** Target binds outside of render passes; use `beginRenderPass` */
const TargetBindPatterns = [
  /\b\w*(Target|target|Framebuffer|framebuffer)(\[[^\]]*\])?\.bind\(\)/,
  /\bbindDrawingBuffer\(/,
];

/** Files binding targets for WebGL readback */
const TargetBindAllowlist: { [file: string]: string } = {
  'packages/graphics/src/gl/render-target.ts':
    '`RenderTarget.bind` implementation (binds its framebuffer for readPixels)',
};

describe('backend-neutral render state', () => {
  const stateViolations: string[] = [];
  const binds = new Map<string, string[]>();
  for (const root of getSourceRoots(['packages/plugin', 'extensions', 'apps'])) {
    for (const file of getSourceFiles(root)) {
      const lines = fs.readFileSync(file, 'utf8').split('\n');
      lines.forEach((line, i) => {
        const location = `${file}:${i + 1}: ${line.trim()}`;
        if (RenderStatePatterns.some((p) => p.test(line))) stateViolations.push(location);
        if (TargetBindPatterns.some((p) => p.test(line))) {
          if (!binds.has(file)) binds.set(file, []);
          binds.get(file)!.push(location);
        }
      });
    }
  }

  it('no WebGL state calls outside gl/webgl/', () => {
    expect(stateViolations).toEqual([]);
  });

  it('no target binds outside gl/webgl/, render passes instead', () => {
    const violations = [...binds].filter(([file]) => !TargetBindAllowlist[file]).flatMap(([, lines]) => lines);
    expect(violations).toEqual([]);
  });

  it('target bind allowlist has no stale entries', () => {
    expect(Object.keys(TargetBindAllowlist).filter((file) => !binds.has(file))).toEqual([]);
  });
});

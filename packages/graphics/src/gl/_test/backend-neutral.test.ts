/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
const Allowlist: { [file: string]: string } = {
    'packages/graphics/src/canvas3d/passes/hi-z.ts': 'fences, readback',
    'packages/graphics/src/canvas3d/passes/pick.ts': 'async readback',
    'packages/graphics/src/gl/compute/grid3d.ts': 'GLSL 100 schema',
    'packages/graphics/src/gl/renderable/direct-volume.ts': 'GLSL 100 loop bounds',
    'extensions/geo-export/src/mesh-exporter.ts': 'readback types',
};

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

/** `gl/webgl/<module>` imports, also of nested paths (e.g. GLSL sources in `shader/`), the first segment as module */
const WebGLImportPattern = /from '(?:@molstar\/graphics\/gl|\.{1,2}(?:\/\.\.)*)\/webgl\/([\w-]+)(?:\/[\w./-]+)?(?:\.js)?'/;

/** `gl/webgl/` modules that neutral code may still import */
const PendingModules: { [module: string]: string } = {
    compat: 'isWebGL2, see the feature check allowlist',
    backend: 'backend descriptors, consumers select the backends they ship',
};

/** Other `gl/webgl/` imports outside the backend */
const ImportAllowlist: { [file: string]: string } = {
    // readback, fences and sync
    'packages/graphics/src/canvas3d/passes/hi-z.ts': 'WebGLContext and PixelPackBuffer: fences, readback',
    'packages/graphics/src/canvas3d/passes/pick.ts': 'WebGLContext and PixelPackBuffer: fences, readback',
    'packages/graphics/src/gl/compute/util.ts': 'WebGLContext and asGLTexture: readback',
    'packages/graphics/src/gl/compute/histogram-pyramid/reduction.ts': 'WebGLContext: gl.finish, sum readback',
    'packages/graphics/src/gl/compute/histogram-pyramid/sum.ts': 'WebGLContext: gl.finish, readPixels',
    'packages/graphics/src/gl/compute/marching-cubes/active-voxels.ts': 'WebGLContext: gl.finish',
    'packages/graphics/src/gl/compute/marching-cubes/isosurface.ts': 'WebGLContext: gl.finish',
    'packages/graphics/src/gl/compute/grid3d.ts': 'WebGLContext: readPixels, sync; GLSL 100 schema',
    'packages/graphics/src/geo/gaussian-density/gpu.ts': 'WebGLContext: readPixels, sync',
    'packages/graphics/src/geo/geometry/texture-mesh/texture-mesh.ts': 'WebGLContext: readPixels',
    'extensions/geo-export/src/mesh-exporter.ts': 'WebGLContext: texture readback',
    'extensions/alpha-orbitals/src/gpu/compute.ts': 'WebGLContext: grid3d GPU compute',
    // pass the WebGLContext on to the WebGL-only functions above
    'packages/graphics/src/geo/gaussian-density.ts': 'WebGLContext for GPU gaussian density',
    'packages/graphics/src/repr/structure/visual/util/gaussian.ts': 'WebGLContext for GPU gaussian density',
    // WebGL-only GPU paths, narrowed with isWebGLContext
    'packages/graphics/src/repr/structure/visual/gaussian-density-volume.ts': 'isWebGLContext: GPU gaussian density',
    'packages/graphics/src/repr/structure/visual/gaussian-surface-mesh.ts': 'isWebGLContext: GPU gaussian surface',
    'packages/graphics/src/repr/volume/isosurface.ts': 'isWebGLContext: GPU isosurface',
    'packages/graphics/src/repr/volume/segment.ts': 'isWebGLContext: GPU isosurface',
    // context creation and other
    'packages/graphics/src/canvas3d/canvas3d.ts': 'createContext, getGLContext, WebGLContext; pixel pack buffers',
    'packages/plugin/headless/src/screenshot.ts': 'createContext for the headless gl package (WebGL-only by design)',
    'packages/graphics/src/gl/renderable/direct-volume.ts': 'isWebGLContext: GLSL 100 loop bounds',
};

describe('backend-neutral imports', () => {
    const imports = new Map<string, string[]>();
    for (const root of getSourceRoots(['packages/plugin', 'extensions', 'apps'])) {
        for (const file of getSourceFiles(root)) {
            const lines = fs.readFileSync(file, 'utf8').split('\n');
            lines.forEach((line, i) => {
                const m = line.match(WebGLImportPattern);
                if (!m || PendingModules[m[1]]) return;
                if (!imports.has(file)) imports.set(file, []);
                imports.get(file)!.push(`${file}:${i + 1}: ${line.trim()}`);
            });
        }
    }

    it('no gl/webgl/ imports outside the backend', () => {
        const violations = [...imports].filter(([file]) => !ImportAllowlist[file]).flatMap(([, lines]) => lines);
        expect(violations).toEqual([]);
    });

    it('import allowlist has no stale entries', () => {
        expect(Object.keys(ImportAllowlist).filter((file) => !imports.has(file))).toEqual([]);
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
const TargetBindPatterns = [/\b\w*(Target|target|Framebuffer|framebuffer)(\[[^\]]*\])?\.bind\(\)/, /\bbindDrawingBuffer\(/];

/** Files binding targets for WebGL readback */
const TargetBindAllowlist: { [file: string]: string } = {
    'packages/graphics/src/canvas3d/passes/pick.ts': 'readPixels, pixel pack buffers',
    'packages/graphics/src/geo/gaussian-density/gpu.ts': 'readPixels',
    'packages/graphics/src/geo/geometry/texture-mesh/texture-mesh.ts': 'readPixels',
    'packages/graphics/src/gl/compute/util.ts': 'readPixels',
    'packages/graphics/src/gl/render-target.ts': '`RenderTarget.bind` implementation (binds its framebuffer for readPixels)',
    'extensions/geo-export/src/mesh-exporter.ts': 'readPixels',
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

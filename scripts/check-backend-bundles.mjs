/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Bundles with a subset of the graphics backends contain no code of the others: builds small entries and the Viewer
// (full, and its variants with one backend, see `esbuild/backend-sets.mjs`) with esbuild (the resolution of the app
// builds) and checks the modules in the output by their metafile, `gl/webgl/` (including the GLSL
// sources) and `gl/webgpu/` (including the WGSL sources). Prints the bundle sizes; fails listing the modules of an
// excluded backend and the imports that pulled them in.
//
// usage: node scripts/check-backend-bundles.mjs [--list]
//
//   --list   lists the backend modules of each bundle with their size in the output

import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { getBackendSetAlias } from './esbuild/backend-sets.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
// an app's directory, so that the workspace packages resolve as in the app builds
const appDir = path.join(root, 'apps/viewer');
const esbuild = createRequire(path.join(root, 'package.json'))('esbuild');
const list = process.argv.includes('--list');

const Backends = {
  webgl: (file) => /[\\/]gl[\\/]webgl[\\/]/.test(file),
  webgpu: (file) => /[\\/]gl[\\/]webgpu[\\/]/.test(file),
};

const Cases = [
  {
    name: 'Canvas3D with the WebGPU backend',
    contents: `
      import { Canvas3D, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
      import { WebGPUBackend } from '@molstar/graphics/gl/webgpu/backend';
      console.log(Canvas3D, Canvas3DContext, WebGPUBackend);`,
    excludes: ['webgl'],
  },
  {
    name: 'Canvas3D with the WebGL backends',
    contents: `
      import { Canvas3D } from '@molstar/graphics/canvas3d/canvas3d';
      import { createWebGLCanvas3DContext } from '@molstar/graphics/canvas3d/webgl';
      import { WebGL2Backend, WebGLBackend } from '@molstar/graphics/gl/webgl/backend';
      console.log(Canvas3D, createWebGLCanvas3DContext, WebGL2Backend, WebGLBackend);`,
    excludes: ['webgpu'],
  },
  {
    name: 'plugin with the default backends (all)',
    contents: `
      import { PluginContext } from '@molstar/plugin/context';
      console.log(PluginContext);`,
    excludes: [],
  },
  {
    name: 'plugin with the default backends aliased to WebGL',
    contents: `
      import { PluginContext } from '@molstar/plugin/context';
      console.log(PluginContext);`,
    alias: getBackendSetAlias('webgl'),
    excludes: ['webgpu'],
  },
  {
    name: 'plugin with the default backends aliased to WebGPU',
    contents: `
      import { PluginContext } from '@molstar/plugin/context';
      console.log(PluginContext);`,
    alias: getBackendSetAlias('webgpu'),
    excludes: ['webgl'],
  },
  {
    name: 'Viewer',
    entry: 'apps/viewer/src/index.ts',
    excludes: [],
  },
  {
    name: 'Viewer with WebGL only (variant webgl)',
    entry: 'apps/viewer/src/index.ts',
    alias: getBackendSetAlias('webgl'),
    excludes: ['webgpu'],
  },
  {
    name: 'Viewer with WebGPU only (variant webgpu)',
    entry: 'apps/viewer/src/index.ts',
    alias: getBackendSetAlias('webgpu'),
    excludes: ['webgl'],
  },
];

const kib = (bytes) => `${(bytes / 1024).toFixed(0)} KiB`;
let failed = false;

for (const c of Cases) {
  const result = await esbuild.build({
    absWorkingDir: appDir,
    ...(c.entry
      ? { entryPoints: [path.join(root, c.entry)] }
      : { stdin: { contents: c.contents, resolveDir: appDir, loader: 'ts' } }),
    bundle: true,
    platform: 'browser',
    format: 'esm',
    write: false,
    outfile: 'bundle.js',
    minify: true,
    metafile: true,
    conditions: ['molstar-src', 'import', 'default'],
    external: ['crypto', 'fs', 'path', 'stream'],
    loader: {
      '.png': 'empty',
      '.jpg': 'empty',
      '.svg': 'empty',
      '.ico': 'empty',
      '.html': 'empty',
      '.css': 'empty',
      '.scss': 'empty',
    },
    alias: c.alias,
    logLevel: 'error',
  });
  const output = Object.values(result.metafile.outputs)[0];
  // modules with code in the output (imported modules whose code was all dropped don't count)
  const included = Object.entries(output.inputs).filter(([, input]) => input.bytesInOutput > 0);
  const parts = Object.entries(Backends).map(([backend, test]) => {
    const modules = included.filter(([file]) => test(file));
    return { backend, modules, bytes: modules.reduce((sum, [, input]) => sum + input.bytesInOutput, 0) };
  });
  const leaks = parts.filter((p) => c.excludes.includes(p.backend) && p.modules.length > 0);
  failed ||= leaks.length > 0;
  console.log(
    `${leaks.length ? 'FAIL' : 'ok  '} ${c.name}: ${kib(output.bytes)} (${parts
      .map((p) => `${p.backend} ${kib(p.bytes)}`)
      .join(', ')})`,
  );
  for (const p of list ? parts : leaks) {
    for (const [file, input] of p.modules.sort((a, b) => b[1].bytesInOutput - a[1].bytesInOutput)) {
      console.log(`       ${p.backend}: ${kib(input.bytesInOutput).padStart(7)} ${file}`);
    }
  }
  for (const leak of leaks) {
    // the imports from outside the backend that pulled it in
    const test = Backends[leak.backend];
    for (const [file, input] of Object.entries(result.metafile.inputs)) {
      if (test(file)) continue;
      for (const i of input.imports) if (test(i.path)) console.log(`       imported by ${file} -> ${i.path}`);
    }
  }
}
process.exit(failed ? 1 : 0);

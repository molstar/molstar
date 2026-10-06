import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';
import { pathToFileURL } from 'node:url';

const require = createRequire(import.meta.url);
function loadNativeModule(name, envName) {
  const modulePath = process.env[envName];
  if (!modulePath)
    throw new Error(
      `Missing native peer '${name}'. Set ${envName} to its installed package entry path, or make it resolvable from the repository node_modules.`,
    );
  try {
    return require(modulePath);
  } catch (error) {
    const message = error?.message ?? String(error);
    if (/Could not locate the bindings file/i.test(message) || error?.code === 'MODULE_NOT_FOUND') {
      throw new Error(`Native peer '${name}' is unavailable at ${modulePath}: ${message}`);
    }
    throw error;
  }
}

const gl = loadNativeModule('gl', 'MOLSTAR_SMOKE_GL');
const pngjs = loadNativeModule('pngjs', 'MOLSTAR_SMOKE_PNGJS');
const { HeadlessPluginContext } = await import('@molstar/plugin-headless/context');
const { DefaultPluginSpec } = await import('@molstar/plugin/default-spec');
const { setFSModule } = await import('@molstar/core/util/data-source');
setFSModule(fs);

const spec = DefaultPluginSpec();
const packageNames = JSON.parse(process.env.MOLSTAR_SMOKE_INTERNAL_PACKAGES ?? '[]');
assert(
  !packageNames.includes('@molstar/mp4-export-extension'),
  'Headless dependency closure must not include the MP4 export extension',
);
const defaultExtensionNames = [
  ...(spec.actions ?? []).map((entry) => entry.action?.id ?? entry.action?.name ?? ''),
  ...(spec.behaviors ?? []).map((entry) => entry.transformer?.id ?? entry.transformer?.name ?? ''),
  ...(spec.animations ?? []).map((animation) => animation?.id ?? animation?.name ?? animation?.constructor?.name ?? ''),
];
for (const name of defaultExtensionNames) {
  assert(!/mp4|h264/i.test(String(name)), `MP4/H264 extension was added to the base plugin spec: ${name}`);
}

let plugin;
try {
  plugin = new HeadlessPluginContext({ gl, pngjs }, spec, { width: 96, height: 80 });
  await plugin.init();
  const structurePath = path.resolve('tiny.pdb');
  const data = await plugin.builders.data.download({ url: pathToFileURL(structurePath).href, isBinary: false });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
  assert(plugin.managers.structure.hierarchy.current.structures.length > 0, 'Headless consumer loaded no structure');
  plugin.canvas3d?.commit(true);
  assert((plugin.canvas3d?.reprCount.value ?? 0) > 0, 'Headless consumer created no render representations');

  const output = path.resolve('headless-smoke.png');
  await plugin.saveImage(output, { width: 96, height: 80 }, undefined, 'png');
  const decoded = pngjs.PNG.sync.read(fs.readFileSync(output));
  assert.equal(decoded.width, 96);
  assert.equal(decoded.height, 80);
  let nonWhitePixels = 0;
  for (let i = 0; i < decoded.data.length; i += 4) {
    if (decoded.data[i] < 245 || decoded.data[i + 1] < 245 || decoded.data[i + 2] < 245) nonWhitePixels++;
  }
  assert(nonWhitePixels > 0, 'Captured PNG contains no rendered structure pixels');
  console.log(`Headless capture passed: ${decoded.width}x${decoded.height}, ${nonWhitePixels} non-white pixels`);
} catch (error) {
  const message = error?.message ?? String(error);
  if (/Could not locate the bindings file/i.test(message) || /Native peer 'gl' is unavailable/i.test(message)) {
    throw new Error(
      `Headless smoke unavailable: ${message}. Install a working node-gl build for this Node.js/runtime, then set MOLSTAR_SMOKE_GL.`,
    );
  }
  throw error;
} finally {
  plugin?.dispose();
}

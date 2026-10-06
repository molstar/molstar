#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Records what the default plugin spec and the Viewer register, as a baseline for the plugin-composition refactor
 * (.v6/plans/plugin-composition.md, step 0 and the acceptance comparison).
 *
 * Usage: node scripts/workspace/registry-dump.mjs [--out <dir>]
 *
 * Requires compiled library output (`pnpm build:lib`). Each target runs in its own process so the global transformer
 * registry only contains what that target imports.
 */

import { execFileSync, execSync } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { registerHooks } from 'node:module';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const targets = ['default', 'viewer'];

function arg(name) {
  const i = process.argv.indexOf(name);
  return i >= 0 ? process.argv[i + 1] : undefined;
}

function lib(path) {
  return pathToFileURL(join(root, path)).href;
}

async function createPlugin(target) {
  if (target === 'default') {
    const { DefaultPluginSpec } = await import(lib('packages/plugin/core/lib/default-spec.js'));
    const { PluginContext } = await import(lib('packages/plugin/core/lib/context.js'));
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    return plugin;
  }
  if (target === 'viewer') {
    // Mirrors Viewer.create with default options, without rendering the UI.
    const { createViewerSpec } = await import(lib('apps/viewer/lib/plugin-spec.js'));
    const { ViewerAutoPreset } = await import(lib('apps/viewer/lib/presets.js'));
    const { PluginUIContext } = await import(lib('packages/plugin/ui/lib/context.js'));
    const plugin = new PluginUIContext(createViewerSpec({}));
    await plugin.init();
    plugin.builders.structure.representation.registerPreset(ViewerAutoPreset);
    return plugin;
  }
  throw new Error(`Unknown target '${target}'`);
}

function names(list) {
  return list.map((e) => e.name);
}

function themes(ctx) {
  return { color: names(ctx.colorThemeRegistry.list), size: names(ctx.sizeThemeRegistry.list) };
}

function actions(manager, transformerByDisplay) {
  // Action ids are random UUIDs; identify transformer-derived actions by transformer id instead.
  const key = (a) => transformerByDisplay.get(a.definition.display) ?? `action:${a.definition.display.name}`;
  const byType = {};
  for (const [type, list] of manager.fromTypeIndex) byType[type.name] = list.map(key);
  return { all: [...manager.actions.values()].map((e) => key(e.action)), byFromType: byType };
}

async function dumpTarget(target) {
  globalThis.window ??= globalThis;
  registerHooks({
    load(url, context, nextLoad) {
      // Bundler-only asset imports (images, styles) resolve to their URL.
      if (/\.(jpe?g|png|gif|svg|webp|css|scss)$/.test(url)) {
        return { format: 'module', source: `export default ${JSON.stringify(url)};`, shortCircuit: true };
      }
      return nextLoad(url, context);
    },
  });

  const plugin = await createPlugin(target);
  const { StateTransformer } = await import(lib('packages/core/lib/state/transformer.js'));
  const transformers = StateTransformer.getAll();
  const transformerByDisplay = new Map(transformers.map((t) => [t.definition.display, t.id]));
  const r = plugin.representation;

  const registry = {
    representations: {
      structure: names(r.structure.registry.list),
      volume: names(r.volume.registry.list),
      particles: names(r.particles.registry.list),
    },
    themes: {
      structure: themes(r.structure.themes),
      volume: themes(r.volume.themes),
      particles: themes(r.particles.themes),
    },
    formats: names(plugin.dataFormats.list),
    presets: {
      hierarchy: plugin.builders.structure.hierarchy._providers.map((p) => p.id),
      representation: plugin.builders.structure.representation._providers.map((p) => p.id),
    },
    selectionQueries: plugin.query.structure.registry.list.map((q) => `${q.category || '(none)'}: ${q.label}`),
    lociLabelProviders: plugin.managers.lociLabels.providers.length,
    markdownExtensions: names(plugin.managers.markdownExtensions.extension),
    dragAndDrop: plugin.managers.dragAndDrop.handlers.map(([name]) => name),
    animations: names(plugin.managers.animation._animations),
    actions: {
      data: actions(plugin.state.data.actions, transformerByDisplay),
      behaviors: actions(plugin.state.behaviors.actions, transformerByDisplay),
    },
    behaviors: plugin.spec.behaviors.map((b) => b.transformer.id),
    customProperties: {
      model: [...plugin.customModelProperties.providers.keys()],
      structure: [...plugin.customStructureProperties.providers.keys()],
      volume: [...plugin.customVolumeProperties.providers.keys()],
    },
  };

  plugin.dispose();
  return { registry, transformerIds: transformers.map((t) => t.id).sort() };
}

function stable(value) {
  return `${JSON.stringify(value, null, 2)}\n`;
}

if (arg('--target')) {
  const result = await dumpTarget(arg('--target'));
  process.stdout.write(`\n@@registry-dump@@${JSON.stringify(result)}\n`);
  // Some optional dependencies (MP4 encoder) keep handles open.
  process.exit(0);
}

const out = resolve(root, arg('--out') ?? '.v6/baselines');
const commit = execSync('git rev-parse --short HEAD', { cwd: root }).toString().trim();
const registry = { commit, targets: {} };
const transformerIds = { commit, targets: {} };

for (const target of targets) {
  const stdout = execFileSync(process.execPath, [fileURLToPath(import.meta.url), '--target', target], {
    cwd: root,
    maxBuffer: 64 * 1024 * 1024,
    stdio: ['ignore', 'pipe', 'inherit'],
  }).toString();
  const marker = stdout.lastIndexOf('@@registry-dump@@');
  if (marker < 0) throw new Error(`No result from target '${target}'`);
  const result = JSON.parse(stdout.slice(marker + '@@registry-dump@@'.length));
  registry.targets[target] = result.registry;
  transformerIds.targets[target] = result.transformerIds;
}

mkdirSync(out, { recursive: true });
writeFileSync(join(out, 'registry.json'), stable(registry));
writeFileSync(join(out, 'transformer-ids.json'), stable(transformerIds));
console.log(`Wrote ${join(out, 'registry.json')} and ${join(out, 'transformer-ids.json')} (commit ${commit})`);

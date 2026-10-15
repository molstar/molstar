/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

/**
 * Per-kind plumbing of `PluginContext.register`: turns entries into registrations, checks them all for conflicts
 * before changing anything, applies them in the fixed order, and builds the undo.
 *
 * This module must not value-import catalogs; it only drives the registries and managers of the context.
 */
export { registerEntries };

/** One registration with the plugin: a provider and the registry or manager that counts it. */
interface Registration {
  /** What the conflict message calls the provider's registry. */
  readonly label: string;
  /** The registry or manager instance, to scope keys. */
  readonly scope: object;
  /** What the registry compares to tell the same provider from a different one under an existing key. */
  readonly identity: object;
  /** Keys the registry detects conflicts on, empty for identity-only registries. */
  readonly keys: readonly string[];
  findConflict(): string | undefined;
  add(): void;
  remove(): void;
}

interface Target<P> {
  readonly label: string;
  readonly scope: object;
  findConflict(provider: P): string | undefined;
  add(provider: P): void;
  remove(provider: P): void;
}

function registration<P>(target: Target<P>, provider: P, identity: object, keys: readonly string[] = []): Registration {
  return {
    label: target.label,
    scope: target.scope,
    identity,
    keys,
    findConflict: () => target.findConflict(provider),
    add: () => target.add(provider),
    remove: () => target.remove(provider),
  };
}

type ThemeRegistryLike<P> = {
  findConflict(p: P): string | undefined;
  add(p: P): void;
  remove(p: P): void;
};

function registryTarget<P>(label: string, registry: ThemeRegistryLike<P>): Target<P> {
  return {
    label,
    scope: registry,
    findConflict: (p) => registry.findConflict(p),
    add: (p) => registry.add(p),
    remove: (p) => registry.remove(p),
  };
}

function targets(plugin: PluginContext) {
  const { structure, volume, particles } = plugin.representation;
  const { hierarchy, representation } = plugin.builders.structure;
  const { dataFormats } = plugin;
  const { lociLabels, markdownExtensions, dragAndDrop, animation } = plugin.managers;
  const queries = plugin.query.structure.registry;
  const actions = plugin.state.data.actions;

  return {
    scopes: {
      structure: {
        colorThemes: registryTarget('Structure color theme', structure.themes.colorThemeRegistry),
        sizeThemes: registryTarget('Structure size theme', structure.themes.sizeThemeRegistry),
        representations: registryTarget('Structure representation', structure.registry),
      },
      volume: {
        colorThemes: registryTarget('Volume color theme', volume.themes.colorThemeRegistry),
        sizeThemes: registryTarget('Volume size theme', volume.themes.sizeThemeRegistry),
        representations: registryTarget('Volume representation', volume.registry),
      },
      particles: {
        colorThemes: registryTarget('Particles color theme', particles.themes.colorThemeRegistry),
        sizeThemes: registryTarget('Particles size theme', particles.themes.sizeThemeRegistry),
        representations: registryTarget('Particles representation', particles.registry),
      },
    },
    formats: registryTarget('Data format', dataFormats),
    hierarchyPresets: {
      label: 'Trajectory hierarchy preset',
      scope: hierarchy,
      findConflict: (p) => hierarchy.findConflict(p),
      add: (p) => hierarchy.registerPreset(p),
      remove: (p) => hierarchy.unregisterPreset(p),
    } satisfies Target<any>,
    representationPresets: {
      label: 'Structure representation preset',
      scope: representation,
      findConflict: (p) => representation.findConflict(p),
      add: (p) => representation.registerPreset(p),
      remove: (p) => representation.unregisterPreset(p),
    } satisfies Target<any>,
    selectionQueries: registryTarget('Structure selection query', queries),
    lociLabels: {
      label: 'Loci label provider',
      scope: lociLabels,
      findConflict: (p) => lociLabels.findConflict(p),
      add: (p) => lociLabels.addProvider(p),
      remove: (p) => lociLabels.removeProvider(p),
    } satisfies Target<any>,
    markdownExtensions: {
      label: 'Markdown extension',
      scope: markdownExtensions,
      findConflict: (p) => markdownExtensions.findConflict(p),
      add: (p) => markdownExtensions.registerExtension(p),
      remove: (p) => markdownExtensions.removeExtension(p),
    } satisfies Target<any>,
    dragAndDrop: {
      label: 'Drag and drop handler',
      scope: dragAndDrop,
      findConflict: (p) => dragAndDrop.findConflict(p),
      add: (p) => dragAndDrop.addEntry(p),
      remove: (p) => dragAndDrop.removeEntry(p),
    } satisfies Target<any>,
    actions: {
      label: 'Action',
      scope: actions,
      findConflict: (a) => actions.findConflict(a),
      add: (a) => actions.add(a),
      remove: (a) => actions.remove(a),
    } satisfies Target<any>,
    animations: {
      label: 'Animation',
      scope: animation,
      findConflict: (a) => animation.findConflict(a),
      add: (a) => animation.register(a),
      remove: (a) => animation.unregister(a),
    } satisfies Target<any>,
  };
}

type Targets = ReturnType<typeof targets>;

/** A `PluginSpec.Action` wrapper carries the action; bare actions and transformers are registered as they are. */
function unwrapAction(a: NonNullable<PluginRegistryEntry['actions']>[number]) {
  return 'definition' in a ? a : a.action;
}

/** The registrations of one entry in the fixed order of the registry contract. */
function entryRegistrations(t: Targets, entry: PluginRegistryEntry, out: Registration[]) {
  const scopes = [
    [entry.structure, t.scopes.structure],
    [entry.volume, t.scopes.volume],
    [entry.particles, t.scopes.particles],
  ] as const;
  for (const [scope, target] of scopes) {
    for (const p of scope?.themes?.color ?? []) out.push(registration(target.colorThemes, p, p, [p.name]));
    for (const p of scope?.themes?.size ?? []) out.push(registration(target.sizeThemes, p, p, [p.name]));
    for (const p of scope?.representations ?? []) out.push(registration(target.representations, p as any, p, [p.name]));
  }

  for (const p of entry.formats ?? []) out.push(registration(t.formats, p, p, [p.name]));

  for (const p of entry.structure?.presets?.hierarchy ?? []) {
    out.push(registration(t.hierarchyPresets, p, p, p.alias === undefined ? [p.id] : [p.id, p.alias]));
  }
  for (const p of entry.structure?.presets?.representation ?? []) {
    out.push(registration(t.representationPresets, p, p, p.alias === undefined ? [p.id] : [p.id, p.alias]));
  }
  for (const q of entry.structure?.selectionQueries ?? []) out.push(registration(t.selectionQueries, q, q));

  for (const p of entry.lociLabels ?? []) out.push(registration(t.lociLabels, p, p));
  for (const p of entry.markdownExtensions ?? []) out.push(registration(t.markdownExtensions, p, p, [p.name]));
  // A drag-and-drop provider's identity is its handler function; entries with the same name and handler are the same.
  for (const p of entry.dragAndDrop ?? []) out.push(registration(t.dragAndDrop, p, p.handle, [p.name]));
  for (const a of entry.actions ?? []) {
    const action = unwrapAction(a);
    // Actions are counted by id and never conflict. The action object (not the wrapper) is the identity.
    out.push(registration(t.actions, action, action));
  }
  for (const p of entry.animations ?? []) out.push(registration(t.animations, p, p, [p.name]));
}

/** Every conflict with existing registrations, and between registrations of the input. */
function findConflicts(registrations: readonly Registration[]) {
  const conflicts: string[] = [];
  const seen = new Map<object, Map<string, Registration>>();
  for (const r of registrations) {
    const existing = r.findConflict();
    if (existing) conflicts.push(existing);

    let keys = seen.get(r.scope);
    if (!keys) seen.set(r.scope, (keys = new Map()));
    for (const key of r.keys) {
      const other = keys.get(key);
      if (other && other.identity !== r.identity) {
        conflicts.push(`${r.label} '${key}' is listed more than once with different providers.`);
      } else if (!other) {
        keys.set(key, r);
      }
    }
  }
  return conflicts;
}

/**
 * Registers `entry` with `plugin` and returns an idempotent undo. Checks every provider for conflicts before changing
 * anything; on a conflict throws one error listing all of them and leaves the plugin unchanged.
 */
function registerEntries(plugin: PluginContext, entry: PluginRegistryEntry | readonly PluginRegistryEntry[]) {
  const t = targets(plugin);
  const list: readonly PluginRegistryEntry[] = Array.isArray(entry)
    ? (entry as readonly PluginRegistryEntry[])
    : [entry as PluginRegistryEntry];

  const registrations: Registration[] = [];
  for (const e of list) entryRegistrations(t, e, registrations);

  const conflicts = findConflicts(registrations);
  if (conflicts.length > 0) {
    throw new Error(
      `PluginContext.register: ${conflicts.length === 1 ? '1 conflict' : `${conflicts.length} conflicts`}, nothing was registered:\n${conflicts.map((c) => `  - ${c}`).join('\n')}`,
    );
  }

  const done: Registration[] = [];
  try {
    for (const r of registrations) {
      r.add();
      done.push(r);
    }
  } catch (e) {
    // A registry rejected something the checks did not predict; leave the plugin as it was.
    for (let i = done.length - 1; i >= 0; i--) done[i].remove();
    throw e;
  }

  let undone = false;
  return () => {
    if (undone) return;
    undone = true;
    for (let i = done.length - 1; i >= 0; i--) done[i].remove();
  };
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction, StateObject, StateTransformer } from '@molstar/core/state';
import { PluginContext } from '@molstar/plugin/context';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import type { PluginStateAnimation } from '@molstar/plugin/state/animation/model';

class Root extends StateObject.factory<{ name: string }>()({ name: 'Register test root' }) {}
let counter = 0;

/** Providers with unique names, so they never collide with anything a plugin registers by itself. */
function providers(tag: string) {
  const named = (name: string, extra?: object) => ({
    name: `${tag}-${name}`,
    label: `${tag}-${name}`,
    category: tag,
    ...extra,
  });
  return {
    color: named('color') as any,
    size: named('size') as any,
    representation: named('representation') as any,
    volumeColor: named('volume-color') as any,
    particlesSize: named('particles-size') as any,
    format: named('format') as any,
    hierarchyPreset: { id: `${tag}-hierarchy`, alias: `${tag}-hierarchy-alias`, display: { name: tag } } as any,
    representationPreset: { id: `${tag}-representation-preset`, display: { name: tag } } as any,
    query: { label: `${tag}-query`, category: tag } as any,
    lociLabel: { label: () => `${tag}-label` },
    markdown: { name: `${tag}-markdown` },
    dragAndDrop: { name: `${tag}-dnd`, handle: () => false },
    action: StateAction.create<Root, void, {}>({ from: [Root], display: { name: `${tag}-action` }, run() {} }),
    transformer: StateTransformer.create<Root, Root, {}>('register-test', {
      name: `${tag}-transformer-${counter++}`,
      from: [Root],
      to: [Root],
      display: { name: `${tag}-transformer` },
      apply: () => new Root({ name: tag }),
    }),
    animation: {
      name: `${tag}-animation`,
      display: { name: `${tag}-animation` },
      params: () => ({}),
      initialState: () => ({}),
      apply: async () => ({ kind: 'finished' }),
    } as unknown as PluginStateAnimation,
  };
}
type Providers = ReturnType<typeof providers>;

function fullEntry(p: Providers): PluginRegistryEntry {
  return {
    structure: {
      themes: { color: [p.color], size: [p.size] },
      representations: [p.representation],
      presets: { hierarchy: [p.hierarchyPreset], representation: [p.representationPreset] },
      selectionQueries: [p.query],
    },
    volume: { themes: { color: [p.volumeColor] } },
    particles: { themes: { size: [p.particlesSize] } },
    formats: [p.format],
    lociLabels: [p.lociLabel],
    markdownExtensions: [p.markdown],
    dragAndDrop: [p.dragAndDrop],
    actions: [p.action, p.transformer],
    animations: [p.animation],
  };
}

async function createPlugin(registry?: PluginRegistryEntry[]) {
  const plugin = new PluginContext({ behaviors: [], registry });
  await plugin.init();
  return plugin;
}

/** Everything `register` can change, in a comparable form. */
function snapshot(plugin: PluginContext) {
  const names = (list: { name: string }[]) => list.map((e) => e.name);
  const { structure, volume, particles } = plugin.representation;
  const scope = (s: typeof structure) => ({
    representations: names(s.registry.list),
    color: names(s.themes.colorThemeRegistry.list),
    size: names(s.themes.sizeThemeRegistry.list),
  });
  return {
    structure: scope(structure),
    volume: scope(volume as typeof structure),
    particles: scope(particles as typeof structure),
    formats: names(plugin.dataFormats.list),
    hierarchyPresets: plugin.builders.structure.hierarchy.providers.map((p) => p.id),
    representationPresets: plugin.builders.structure.representation.providers.map((p) => p.id),
    queries: plugin.query.structure.registry.list.map((q) => q.label),
    lociLabels: plugin.managers.lociLabels.providers.length,
    markdown: names(plugin.managers.markdownExtensions.extension),
    dragAndDrop: plugin.managers.dragAndDrop.list().map((e) => e.name),
    animations: names(plugin.managers.animation.animations),
    actions: (plugin.state.data.actions as any).actions.size as number,
  };
}

/** Records every registration call, in order, as `kind:provider`. */
function trace(plugin: PluginContext) {
  const log: string[] = [];
  const patch = (target: any, method: string, kind: string, nameOf: (p: any) => string) => {
    const original = target[method].bind(target);
    target[method] = (p: any, ...rest: any[]) => {
      log.push(`${kind}:${nameOf(p)}`);
      return original(p, ...rest);
    };
  };
  const byName = (p: any) => p.name;
  for (const scope of ['structure', 'volume', 'particles'] as const) {
    const s = plugin.representation[scope];
    patch(s.themes.colorThemeRegistry, 'add', `${scope}.color`, byName);
    patch(s.themes.sizeThemeRegistry, 'add', `${scope}.size`, byName);
    patch(s.registry, 'add', `${scope}.representation`, byName);
  }
  patch(plugin.dataFormats, 'add', 'format', byName);
  patch(plugin.builders.structure.hierarchy, 'registerPreset', 'hierarchy', (p) => p.id);
  patch(plugin.builders.structure.representation, 'registerPreset', 'representation-preset', (p) => p.id);
  patch(plugin.query.structure.registry, 'add', 'query', (p) => p.label);
  patch(plugin.managers.lociLabels, 'addProvider', 'loci-label', (p) => p.label());
  patch(plugin.managers.markdownExtensions, 'registerExtension', 'markdown', byName);
  patch(plugin.managers.dragAndDrop, 'addEntry', 'drag-and-drop', byName);
  patch(plugin.state.data.actions, 'add', 'action', (a) => a.definition.display.name);
  patch(plugin.managers.animation, 'register', 'animation', byName);
  return log;
}

describe('PluginContext.register', () => {
  it('throws before init()', () => {
    const plugin = new PluginContext({ behaviors: [] });
    expect(() => plugin.register({})).toThrow('PluginContext.register called before init()');
    expect(() => plugin.register([])).toThrow('PluginContext.register called before init()');
  });

  it('registers one entry, or a list of entries, and returns an undo', async () => {
    const plugin = await createPlugin();
    const before = snapshot(plugin);
    const a = providers('a');
    const b = providers('b');

    const undoOne = plugin.register(fullEntry(a));
    expect(snapshot(plugin).formats).toContain('a-format');
    const undoList = plugin.register([fullEntry(b)]);
    expect(snapshot(plugin).formats).toContain('b-format');

    undoList();
    undoOne();
    expect(snapshot(plugin)).toEqual(before);
  });

  it('registers in the fixed order within an entry and in input order across entries', async () => {
    const plugin = await createPlugin();
    const log = trace(plugin);
    const a = providers('a');
    const b = providers('b');

    plugin.register([fullEntry(a), { formats: [b.format], structure: { themes: { color: [b.color] } } }]);

    expect(log).toEqual([
      'structure.color:a-color',
      'structure.size:a-size',
      'structure.representation:a-representation',
      'volume.color:a-volume-color',
      'particles.size:a-particles-size',
      'format:a-format',
      'hierarchy:a-hierarchy',
      'representation-preset:a-representation-preset',
      'query:a-query',
      'loci-label:a-label',
      'markdown:a-markdown',
      'drag-and-drop:a-dnd',
      'action:a-action',
      'action:a-transformer',
      'animation:a-animation',
      // second entry: scopes (themes), then formats, even though it lists them the other way around
      'structure.color:b-color',
      'format:b-format',
    ]);
  });

  it('registers the action carried by a PluginSpec.Action wrapper', async () => {
    const plugin = await createPlugin();
    const log = trace(plugin);
    const a = providers('a');
    const undo = plugin.register({ actions: [{ action: a.action }, { action: a.transformer }] });
    expect(log).toEqual(['action:a-action', 'action:a-transformer']);

    const before = snapshot(plugin).actions;
    // The same action bare or wrapped counts under one id.
    plugin.register({ actions: [a.action] });
    expect(snapshot(plugin).actions).toBe(before);
    undo();
    expect(snapshot(plugin).actions).toBe(before - 1);
  });

  it('counts the same provider object registered by several calls', async () => {
    const plugin = await createPlugin();
    const before = snapshot(plugin);
    const a = providers('a');
    const undo1 = plugin.register(fullEntry(a));
    const undo2 = plugin.register(fullEntry(a));

    undo1();
    expect(snapshot(plugin).formats).toContain('a-format');
    expect(snapshot(plugin).queries).toContain('a-query');
    undo2();
    expect(snapshot(plugin)).toEqual(before);
  });

  describe('atomic conflict check', () => {
    it('throws one error listing every conflict and changes nothing', async () => {
      const plugin = await createPlugin();
      const a = providers('a');
      plugin.register(fullEntry(a));
      const before = snapshot(plugin);

      // Same keys as in `a`, but different objects, plus providers that would be fine on their own.
      const clash = providers('a');
      const fresh = providers('fresh');
      const entry: PluginRegistryEntry = {
        ...fullEntry(fresh),
        structure: {
          themes: { color: [fresh.color, clash.color], size: [clash.size] },
          representations: [clash.representation],
          presets: { hierarchy: [clash.hierarchyPreset], representation: [clash.representationPreset] },
          selectionQueries: [fresh.query],
        },
        volume: { themes: { color: [clash.volumeColor] } },
        formats: [fresh.format, clash.format],
        markdownExtensions: [clash.markdown],
        dragAndDrop: [clash.dragAndDrop],
        animations: [clash.animation],
      };

      let message = '';
      try {
        plugin.register(entry);
      } catch (e) {
        message = (e as Error).message;
      }

      expect(message).toMatch(/^PluginContext\.register: 10 conflicts, nothing was registered:/);
      for (const text of [
        /Theme 'a-color'/,
        /Theme 'a-size'/,
        /Representation 'a-representation'/,
        /'a-hierarchy'/,
        /'a-representation-preset'/,
        /Theme 'a-volume-color'/,
        /'a-format'/,
        /'a-markdown'/,
        /'a-dnd'/,
        /'a-animation'/,
      ]) {
        expect(message).toMatch(text);
      }
      expect(snapshot(plugin)).toEqual(before);
    });

    it('does not register anything before a conflict later in the input', async () => {
      const plugin = await createPlugin();
      const a = providers('a');
      plugin.register({ animations: [a.animation] });
      const before = snapshot(plugin);

      const fresh = providers('fresh');
      const clash = providers('a');
      expect(() => plugin.register([fullEntry(fresh), { animations: [clash.animation] }])).toThrow(/'a-animation'/);
      expect(snapshot(plugin)).toEqual(before);

      // The fresh providers are still free to register.
      expect(() => plugin.register(fullEntry(fresh))).not.toThrow();
    });

    it('reports conflicts between providers of the input', async () => {
      const plugin = await createPlugin();
      const before = snapshot(plugin);
      const a = providers('a');
      const same = providers('a');

      expect(() => plugin.register({ formats: [a.format, same.format] })).toThrow(
        /Data format 'a-format' is listed more than once with different providers/,
      );
      expect(() => plugin.register([{ formats: [a.format] }, { formats: [same.format] }])).toThrow(
        /Data format 'a-format' is listed more than once/,
      );
      expect(() => plugin.register({ animations: [a.animation, same.animation] })).toThrow(/'a-animation'/);
      expect(() => plugin.register({ dragAndDrop: [a.dragAndDrop, same.dragAndDrop] })).toThrow(/'a-dnd'/);
      expect(snapshot(plugin)).toEqual(before);
    });

    it('detects an alias of one preset clashing with the id or alias of another in the input', async () => {
      const plugin = await createPlugin();
      const before = snapshot(plugin);
      const preset = (id: string, alias?: string) => ({ id, alias, display: { name: id } }) as any;

      expect(() =>
        plugin.register({ structure: { presets: { hierarchy: [preset('one', 'x'), preset('two', 'x')] } } }),
      ).toThrow(/Trajectory hierarchy preset 'x' is listed more than once/);
      expect(() =>
        plugin.register({ structure: { presets: { representation: [preset('one'), preset('two', 'one')] } } }),
      ).toThrow(/Structure representation preset 'one' is listed more than once/);
      expect(snapshot(plugin)).toEqual(before);
    });

    it('accepts the same object listed twice in the input', async () => {
      const plugin = await createPlugin();
      const before = snapshot(plugin);
      const a = providers('a');
      const undo = plugin.register({
        formats: [a.format, a.format],
        structure: { selectionQueries: [a.query, a.query] },
      });
      expect(snapshot(plugin).formats).toContain('a-format');
      undo();
      expect(snapshot(plugin)).toEqual(before);
    });

    it('treats drag-and-drop entries with the same name and handler as the same provider', async () => {
      const plugin = await createPlugin();
      const a = providers('a');
      plugin.register({ dragAndDrop: [a.dragAndDrop] });
      expect(() => plugin.register({ dragAndDrop: [{ ...a.dragAndDrop }] })).not.toThrow();
    });
  });

  describe('undo', () => {
    it('decrements exactly the registrations it made and is idempotent', async () => {
      const plugin = await createPlugin();
      const a = providers('a');
      // Registered by someone else first: the undo must leave this one alone.
      const other = plugin.register({ formats: [a.format], animations: [a.animation] });

      const undo = plugin.register(fullEntry(a));
      undo();
      const after = snapshot(plugin);
      expect(after.formats).toContain('a-format');
      expect(after.animations).toContain('a-animation');
      expect(after.queries).not.toContain('a-query');
      expect(after.structure.color).not.toContain('a-color');

      // A second call is a no-op: it must not take the other registration's count.
      undo();
      expect(snapshot(plugin)).toEqual(after);

      other();
      expect(snapshot(plugin).formats).not.toContain('a-format');
      expect(snapshot(plugin).animations).not.toContain('a-animation');
    });

    it('is a no-op for an empty entry', async () => {
      const plugin = await createPlugin();
      const before = snapshot(plugin);
      plugin.register({})();
      plugin.register([])();
      expect(snapshot(plugin)).toEqual(before);
    });
  });

  describe('init()', () => {
    it('registers spec.registry before behaviors are initialized', async () => {
      const a = providers('a');
      const b = providers('b');
      const plugin = new PluginContext({ behaviors: [], registry: [fullEntry(a), { formats: [b.format] }] });
      await plugin.init();

      const s = snapshot(plugin);
      expect(s.formats).toEqual(expect.arrayContaining(['a-format', 'b-format']));
      expect(s.structure.color).toContain('a-color');
      expect(s.queries).toContain('a-query');
      expect(s.animations).toContain('a-animation');
      expect(plugin.managers.dragAndDrop.list().map((e) => e.name)).toContain('a-dnd');
    });

    it('rejects init() on a conflict between spec.registry entries', async () => {
      const a = providers('a');
      const same = providers('a');
      const plugin = new PluginContext({
        behaviors: [],
        registry: [{ formats: [a.format] }, { formats: [same.format] }],
      });
      plugin.initialized.catch(() => {});
      await expect(plugin.init()).rejects.toThrow(/'a-format'/);
      expect(plugin.isInitialized).toBe(false);
    });

    it.each(['actions', 'animations', 'customFormats'])('rejects the removed spec field %s', (key) => {
      const spec = { behaviors: [], [key]: [] } as any;
      expect(() => new PluginContext(spec)).toThrow(
        `PluginSpec.${key} was removed in 6.0; use registry entries (see the migration guide)`,
      );
    });

    it('ignores removed spec keys whose value is undefined', async () => {
      const spec = { behaviors: [], actions: undefined, animations: undefined, customFormats: undefined } as any;
      const plugin = new PluginContext(spec);
      await plugin.init();
      expect(plugin.managers.animation.animations).toEqual([]);
      expect(plugin.managers.animation.current).toBeUndefined();
    });

    it('does not let a behavior-style removal remove a provider the spec also registered', async () => {
      const a = providers('a');
      const plugin = await createPlugin([fullEntry(a)]);

      // A behavior registers what the spec already lists, then unregisters it again.
      plugin.managers.lociLabels.addProvider(a.lociLabel);
      plugin.managers.lociLabels.removeProvider(a.lociLabel);
      plugin.managers.animation.register(a.animation);
      plugin.managers.animation.unregister(a.animation);
      plugin.dataFormats.add(a.format);
      plugin.dataFormats.remove(a.format);
      plugin.managers.markdownExtensions.registerExtension(a.markdown);
      plugin.managers.markdownExtensions.removeExtension(a.markdown);
      plugin.managers.dragAndDrop.addEntry(a.dragAndDrop);
      plugin.managers.dragAndDrop.removeEntry(a.dragAndDrop);
      plugin.state.data.actions.add(a.action);
      plugin.state.data.actions.remove(a.action);

      const s = snapshot(plugin);
      expect(s.lociLabels).toBeGreaterThan(0);
      expect(plugin.managers.lociLabels.providers).toContain(a.lociLabel);
      expect(s.animations).toContain('a-animation');
      expect(s.formats).toContain('a-format');
      expect(s.markdown).toContain('a-markdown');
      expect(s.dragAndDrop).toContain('a-dnd');
      expect(
        plugin.state.data.actions.fromCell({ obj: new Root({ name: 'x' }), transform: {} } as any, plugin),
      ).toContain(a.action);
    });
  });
});

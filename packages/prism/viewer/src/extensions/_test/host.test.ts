/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginSpec } from '@molstar/plugin/spec';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { ApiArea, type ApiAreaContext, ApiResult } from '../../api/area.js';
import type { LazyViewerExtension, ViewerExtension } from '../../spec.js';

// The extension API (viewer-spec.md §3–§4): eager and lazy extensions, stubs, replay, overrides.

/** A behavior that does nothing, for the plugin part of a lazy extension. */
const Quiet = PluginBehavior.create({
  name: 'prism-test-quiet',
  category: 'misc',
  display: { name: 'Quiet' },
  ctor: class extends PluginBehavior.Handler {
    register() {}
  },
});

class ShoutApi extends ApiArea {
  readonly id = 'shout';
  readonly say = this.operation({
    id: 'shout.say',
    title: 'Say something',
    params: { text: PD.Text('hello') },
    async run({ text }, ctx) {
      // A real step: it loads a file through another operation, part of this one.
      const r = await ctx.api.load.files({ files: [cbs()] });
      return r.ok ? ApiResult.ok(`Said ${text}.`, { nodes: [] }) : r;
    },
  });
}

const shout: ViewerExtension = {
  id: 'shout',
  title: 'Shout',
  plugin: { behaviors: [PluginSpec.Behavior(Quiet)] },
  api: { shout: (ctx: ApiAreaContext) => new ShoutApi(ctx) },
  commands: (api) => [
    {
      id: 'shout.say',
      title: 'Say something',
      category: 'Shout',
      run: async () => {
        const r = await api.run('shout.say', {});
        return { message: r.ok ? r.did : 'failed' };
      },
    },
  ],
};

/** A lazy `shout` with its own load counter (`count()`), so tests do not share state. */
const lazyShout = (o: Partial<LazyViewerExtension> = {}) => {
  let loads = 0;
  const extension: LazyViewerExtension = {
    id: 'shout',
    title: 'Shout',
    stub: {
      parts: [],
      operations: [{ id: 'shout.say', title: 'Say something' }],
      commands: [{ id: 'shout.say', title: 'Say something' }],
      owns: { fileExtensions: ['shout'] },
    },
    load: async () => {
      loads++;
      return shout;
    },
    ...o,
  };
  return { extension, count: () => loads };
};

describe('extensions', () => {
  it('loads a lazy extension at its first call, and replays the call as one step', async () => {
    const lazy = lazyShout();
    const v = await headlessViewer({ extensions: [lazy.extension] });
    expect(v.extensions.status('shout')?.value).toBe('stub');
    expect(v.api.describe('shout.say')).toBeDefined();
    const r = await v.api.run('shout.say', { text: 'hi' });
    expectOk(r);
    expect(r.did).toBe('Said hi.');
    expect(lazy.count()).toBe(1);
    expect(v.extensions.status('shout')?.value).toBe('loaded');
    expect(v.extensions.loaded.value).toContain('shout');
    // Its plugin part is on the running plugin.
    expect(v.plugin.state.behaviors.transforms.has(Quiet.id)).toBe(true);
    // The replayed call and the file it loaded are one step.
    expect(v.managers.history.steps.value).toHaveLength(1);
    // A second call runs the real area directly.
    expectOk(await v.api.run('shout.say', { text: 'again' }));
    expect(lazy.count()).toBe(1);
  });

  it("makes the replayed call's step with the real label and the caller's origin, and Undo removes what it loaded", async () => {
    const lazy = lazyShout();
    const v = await headlessViewer({ extensions: [lazy.extension] });
    expectOk(await v.api.as({ kind: 'person' }).run('shout.say', { text: 'hi' }));
    const steps = v.managers.history.steps.value;
    expect(steps).toHaveLength(1);
    expect(steps[0].label).toBe('Said hi');
    expect(steps[0].origin).toEqual({ kind: 'person' });
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(1);

    expectOk(await v.api.history.undo());
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(0);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('loads once for two first calls at once, and makes two steps', async () => {
    const lazy = lazyShout();
    const v = await headlessViewer({ extensions: [lazy.extension] });
    const [a, b] = await Promise.all([v.api.run('shout.say', { text: 'a' }), v.api.run('shout.say', { text: 'b' })]);
    expectOk(a);
    expectOk(b);
    expect(lazy.count()).toBe(1);
    expect(v.extensions.status('shout')?.value).toBe('loaded');
    expect(v.managers.history.steps.value).toHaveLength(2);
    expect(v.managers.history.steps.value.map((s) => s.label).sort()).toEqual(['Said a', 'Said b']);
  });

  it('leaves no step and no structure when the replayed call fails', async () => {
    // The real operation fails after it loaded a file: the step restores the checkpoint.
    class FailingApi extends ApiArea {
      readonly id = 'shout';
      readonly say = this.operation({
        id: 'shout.say',
        title: 'Say something',
        params: { text: PD.Text('hello') },
        async run(_, ctx) {
          expectOk(await ctx.api.load.files({ files: [cbs()] }));
          return ApiResult.fail('failed', 'It broke.');
        },
      });
    }
    const failing: ViewerExtension = { ...shout, api: { shout: (ctx: ApiAreaContext) => new FailingApi(ctx) } };
    const lazy = lazyShout({ load: async () => failing });
    const v = await headlessViewer({ extensions: [lazy.extension] });
    const r = await v.api.run('shout.say', { text: 'hi' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(v.managers.history.steps.value).toHaveLength(0);
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(0);
  });

  it("lists a lazy extension's commands before it loads; running one loads it", async () => {
    const lazy = lazyShout();
    const v = await headlessViewer({ extensions: [lazy.extension] });
    expect(v.managers.commands.get('shout.say')?.title).toBe('Say something');
    const opened = v.managers.commands.open('shout.say');
    if (opened.kind === 'ran') await opened.result;
    expect(lazy.count()).toBe(1);
    expect(v.extensions.status('shout')?.value).toBe('loaded');
  });

  it('loads what a file needs first, and preloads when asked', async () => {
    const v = await headlessViewer({ extensions: [lazyShout().extension] });
    await v.extensions.ensureFor({ fileExtensions: ['SHOUT'] });
    expect(v.extensions.status('shout')?.value).toBe('loaded');
    const preloaded = lazyShout({ preload: true });
    const p = await headlessViewer({ extensions: [preloaded.extension] });
    expect(p.extensions.status('shout')?.value).toBe('loaded');
    expect(preloaded.count()).toBe(1);
  });

  it('fails a load that brings what its stub did not declare, and says why', async () => {
    const v = await headlessViewer({
      extensions: [lazyShout({ stub: { parts: [], operations: [{ id: 'other.op', title: 'X' }] } }).extension],
    });
    const r = await v.api.run('other.op', {});
    expect(r).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.extensions.status('shout')?.value).toBe('failed');
  });

  it('puts nothing of a lazy extension on a viewer disposed while its chunk loads', async () => {
    let release!: () => void;
    const gate = new Promise<void>((r) => {
      release = r;
    });
    let areas = 0;
    const counted: ViewerExtension = {
      ...shout,
      api: {
        shout: (ctx: ApiAreaContext) => {
          areas++;
          return new ShoutApi(ctx);
        },
      },
    };
    const lazy = lazyShout({
      load: async () => {
        await gate;
        return counted;
      },
    });
    const v = await headlessViewer({ extensions: [lazy.extension] });
    const loading = v.extensions.load('shout');
    await v.dispose();
    release();
    expect(await loading).toBe(false);
    expect(areas).toBe(0);
    expect(v.extensions.loaded.value).not.toContain('shout');
    expect(await v.extensions.load('shout')).toBe(false);
  });

  it('removes or replaces a command by id', async () => {
    const v = await headlessViewer({ extensions: [shout], spec: { overrides: { 'shout.say': null } } });
    expect(v.managers.commands.get('shout.say')).toBeUndefined();
  });
});

describe('extensions mounted at runtime (viewer-spec.md §3)', () => {
  it('mounts an extension the spec did not list, with its plugin part, area and commands; unmounting removes them', async () => {
    const v = await headlessViewer();
    expect(v.api.area('shout')).toBeUndefined();
    const unmount = await v.extensions.mount(shout);
    expect(v.extensions.loaded.value).toContain('shout');
    expect(v.managers.commands.get('shout.say')).toBeDefined();
    expectOk(await v.api.run('shout.say', {}));
    await unmount();
    expect(v.api.area('shout')).toBeUndefined();
    expect(v.managers.commands.get('shout.say')).toBeUndefined();
    expect(v.extensions.loaded.value).not.toContain('shout');
  });

  it('mounts a lazy one as its stub, loading on the first call', async () => {
    const v = await headlessViewer();
    const { extension, count } = lazyShout();
    await v.extensions.mount(extension);
    expect(count()).toBe(0);
    expectOk(await v.api.run('shout.say', {}));
    expect(count()).toBe(1);
  });

  it('refuses an id the viewer has, and an area another extension has, naming it', async () => {
    const v = await headlessViewer({ extensions: [shout] });
    await expect(v.extensions.mount(shout)).rejects.toThrow(/in the viewer already/);
    const echo: ViewerExtension = { id: 'echo', title: 'Echo', api: { shout: (ctx) => new ShoutApi(ctx) } };
    await expect(v.extensions.mount(echo)).rejects.toThrow(/'shout', which extension 'shout' has already/);
    const core: ViewerExtension = { id: 'core-clash', title: 'Clash', api: { canvas: (ctx) => new ShoutApi(ctx) } };
    await expect(v.extensions.mount(core)).rejects.toThrow(/'canvas', which the viewer has already/);
    // Nothing of a refused one stays.
    expect(v.extensions.loaded.value).not.toContain('echo');
  });
});

describe("extensions' command-menu prefixes", () => {
  const prefixed = (id: string, char: string): ViewerExtension => ({
    id,
    title: id,
    commandPrefixes: ({ api }) => [
      { char, title: `${id} words`, items: (text) => [{ id: `${id}:${text}`, label: text, run: () => void api }] },
    ],
  });

  it('lists an extension’s prefixes while it is mounted, and none after', async () => {
    const v = await headlessViewer();
    const unmount = await v.extensions.mount(prefixed('echo', 'echo '));
    expect(v.extensions.prefixes.value.map((p) => p.char)).toEqual(['echo ']);
    expect(v.extensions.prefixes.value[0].items('hi').map((i) => i.label)).toEqual(['hi']);
    await unmount();
    expect(v.extensions.prefixes.value).toEqual([]);
  });

  it('refuses a second extension with a prefix another has, and keeps nothing of it', async () => {
    const v = await headlessViewer({ extensions: [prefixed('echo', 'say ')] });
    await expect(v.extensions.mount(prefixed('other', 'say '))).rejects.toThrow(/prefix 'say ', which another has/);
    expect(v.extensions.prefixes.value.map((p) => p.char)).toEqual(['say ']);
    expect(v.extensions.loaded.value).not.toContain('other');
  });
});

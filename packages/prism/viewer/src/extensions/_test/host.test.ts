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

// The extension API (viewer-spec.md §3–§4, plan 3.1): eager and lazy extensions, stubs, replay, overrides.

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

let loads = 0;
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

const lazyShout = (o: Partial<LazyViewerExtension> = {}): LazyViewerExtension => ({
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
});

describe('extensions', () => {
  it('loads a lazy extension at its first call, and replays the call as one step', async () => {
    loads = 0;
    const v = await headlessViewer({ extensions: [lazyShout()] });
    expect(v.extensions.status('shout')?.value).toBe('stub');
    expect(v.api.describe('shout.say')).toBeDefined();
    const r = await v.api.run('shout.say', { text: 'hi' });
    expectOk(r);
    expect(r.did).toBe('Said hi.');
    expect(loads).toBe(1);
    expect(v.extensions.status('shout')?.value).toBe('loaded');
    expect(v.extensions.loaded.value).toContain('shout');
    // Its plugin part is on the running plugin.
    expect(v.plugin.state.behaviors.transforms.has(Quiet.id)).toBe(true);
    // The replayed call and the file it loaded are one step.
    expect(v.managers.history.steps.value).toHaveLength(1);
    // A second call runs the real area directly.
    expectOk(await v.api.run('shout.say', { text: 'again' }));
    expect(loads).toBe(1);
  });

  it("lists a lazy extension's commands before it loads; running one loads it", async () => {
    loads = 0;
    const v = await headlessViewer({ extensions: [lazyShout()] });
    expect(v.managers.commands.get('shout.say')?.title).toBe('Say something');
    const opened = v.managers.commands.open('shout.say');
    if (opened.kind === 'ran') await opened.result;
    expect(loads).toBe(1);
    expect(v.extensions.status('shout')?.value).toBe('loaded');
  });

  it('loads what a file needs first, and preloads when asked', async () => {
    loads = 0;
    const v = await headlessViewer({ extensions: [lazyShout()] });
    await v.extensions.ensureFor({ fileExtensions: ['SHOUT'] });
    expect(v.extensions.status('shout')?.value).toBe('loaded');
    loads = 0;
    const p = await headlessViewer({ extensions: [lazyShout({ preload: true })] });
    expect(p.extensions.status('shout')?.value).toBe('loaded');
    expect(loads).toBe(1);
  });

  it('fails a load that brings what its stub did not declare, and says why', async () => {
    const v = await headlessViewer({
      extensions: [lazyShout({ stub: { parts: [], operations: [{ id: 'other.op', title: 'X' }] } })],
    });
    const r = await v.api.run('other.op', {});
    expect(r).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.extensions.status('shout')?.value).toBe('failed');
  });

  it('removes or replaces a command by id', async () => {
    const { DefaultPrismViewerSpec } = await import('../../default-spec.js');
    const { createPrismViewerHeadless } = await import('../../headless.js');
    const v = await createPrismViewerHeadless({
      spec: { ...DefaultPrismViewerSpec(), extensions: [shout], overrides: { 'shout.say': null } },
    });
    expect(v.managers.commands.get('shout.say')).toBeUndefined();
    await v.dispose();
  });
});

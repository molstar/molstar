/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { StructureFeature } from '../../features/structure/index.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { ViewerExtension } from '../../spec.js';
import type { PrismViewer, ViewerPointerEvent } from '../../viewer.js';
import { ApiArea, type ApiAreaContext, ApiResult } from '../area.js';

// The run path (operations.md §1, api.md §5), through the viewer's API on a headless plugin.

const dataDir = path.resolve(import.meta.dirname, '../../../../../../data/examples');
const crambin = () => new File([fs.readFileSync(path.join(dataDir, '1crn.cif'))], '1crn.cif');
const cbs = () => new File([fs.readFileSync(path.join(dataDir, '1cbs_updated.cif'))], '1cbs.cif');

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

/** A gate the test opens. */
function gate() {
  let open!: () => void;
  const p = new Promise<void>((r) => (open = r));
  return { p, open };
}

const gates: { current?: ReturnType<typeof gate> } = {};

/** An extension's area whose operations call others, wait, and fail on purpose. */
class DemoApi extends ApiArea {
  readonly id = 'demo';
  readonly loadBoth = this.operation({
    id: 'demo.load-both',
    title: 'Load both',
    params: {},
    async run(_, ctx) {
      const a = await ctx.api.load.files({ files: [crambin()] });
      const b = await ctx.api.load.files({ files: [cbs()] });
      return a.ok && b.ok ? ApiResult.ok('Loaded both.', { nodes: [] }) : ApiResult.fail('failed', 'One failed.');
    },
  });
  readonly waitThenFail = this.operation({
    id: 'demo.wait-then-fail',
    title: 'Wait, then fail',
    params: {},
    async run(_, ctx) {
      await ctx.api.load.files({ files: [crambin()] });
      await gates.current?.p;
      return ApiResult.fail('failed', 'Failed on purpose.');
    },
  });
  readonly slow = this.operation({
    id: 'demo.slow',
    title: 'Slow',
    params: {},
    async run(_, ctx) {
      await ctx.api.load.files({ files: [crambin()] });
      await gates.current?.p;
      return ApiResult.ok('Done.', { nodes: [] });
    },
  });
}

const demo: ViewerExtension = { id: 'demo', title: 'Demo', api: { demo: (ctx: ApiAreaContext) => new DemoApi(ctx) } };

// The scene's API (the \`structure\` feature) and the test's extensions; none of the other features.
async function create(extensions: ViewerExtension[] = [demo]) {
  viewer = await createPrismViewerHeadless({
    spec: { ...DefaultPrismViewerSpec(), extensions: [StructureFeature, ...extensions] },
  });
  return viewer;
}

const demoApi = (v: PrismViewer) => v.api.area('demo') as DemoApi;

describe('the run path', () => {
  it('runs an operation that calls others as one step, without waiting on itself', async () => {
    const v = await create();
    const r = await demoApi(v).loadBoth();
    expect(r.ok).toBe(true);
    expect(v.managers.history.steps.value).toHaveLength(1);
    expect(v.managers.history.steps.value[0].runs).toHaveLength(3);
    expect(v.api.structure.list.value).toHaveLength(2);
    await v.api.history.undo();
    expect(v.state.empty.value).toBe(true);
  });

  it('restores everything a failed operation and its calls did', async () => {
    const v = await create();
    gates.current = gate();
    gates.current.open();
    const r = await demoApi(v).waitThenFail();
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(v.state.empty.value).toBe(true);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it("keeps another origin's call out of a running operation: it waits, and stands as its own step", async () => {
    const v = await create();
    gates.current = gate();
    const running = demoApi(v).waitThenFail();
    await new Promise((r) => setTimeout(r, 20));
    const click = v.api.as({ kind: 'person' }).load.files({ files: [cbs()] });
    await new Promise((r) => setTimeout(r, 20));
    gates.current.open();
    expect((await running).ok).toBe(false);
    expect((await click).ok).toBe(true);
    expect(v.api.structure.list.value.map((s) => s.label)).toEqual(['1CBS']);
    expect(v.managers.history.steps.value).toHaveLength(1);
  });

  it('records no step for what is already so, or refused', async () => {
    const v = await create();
    await v.api.load.files({ files: [crambin()] });
    const same = await v.api.structure.model({ index: 0 });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    const refused = await v.api.structure.model({ index: 9 });
    expect(refused.ok).toBe(false);
    expect(v.managers.history.steps.value).toHaveLength(1);
  });

  it("joins an agent's calls into its turn's entry", async () => {
    const v = await create();
    const agent = v.api.as({ kind: 'agent', turn: 't1' });
    const end = v.managers.history.beginGroup('Load two', { kind: 'agent', turn: 't1' });
    await agent.load.files({ files: [crambin()] });
    await agent.load.files({ files: [cbs()] });
    end();
    expect(v.managers.history.steps.value.map((s) => s.label)).toEqual(['Load two']);
    const records: string[] = [];
    v.api.events.finished.subscribe((r) =>
      records.push(`${r.origin.kind}:${(r.origin as { turn?: string }).turn ?? ''}`),
    );
    await agent.structure.get();
    expect(records).toEqual(['agent:t1']);
  });

  it('makes a group all or nothing, with its calls and what they changed', async () => {
    const v = await create();
    const ok = await v.api.group('Two', async (api) => {
      await api.load.files({ files: [crambin()] });
      await api.load.files({ files: [cbs()] });
    });
    expect(ok.ok).toBe(true);
    expect(ok.runs).toHaveLength(2);
    expect(ok.ok && ok.changed.nodes).toHaveLength(2);

    // A call inside that fails, ignored by the work, still fails the group.
    const bad = await v.api.group('Ignores a failure', async (api) => {
      await api.session.clear();
      await api.structure.assembly({ id: '1', structure: '1crn' });
    });
    // The first failure inside is the group's: 1CRN was cleared before it was asked for.
    expect(bad).toMatchObject({ ok: false, code: 'not-found' });
    expect(v.api.structure.list.value).toHaveLength(2);

    // A call refused before it ran (invalid values: no history step of its own) fails the group too, and what came
    // before it in the group is undone.
    const refused = await v.api.group('Refused', async (api) => {
      await api.session.clear();
      await api.structure.preset({ preset: 'no-such-preset' as never });
    });
    expect(refused.ok).toBe(false);
    expect(v.api.structure.list.value).toHaveLength(2);
    // So does a call to an operation that does not exist.
    const unknown = await v.api.group('Unknown', async (api) => {
      await api.session.clear();
      await api.run('no.such-operation', {});
    });
    expect(unknown).toMatchObject({ ok: false, code: 'not-found' });
    expect(v.api.structure.list.value).toHaveLength(2);
    // With onFailure 'continue', a refused call is passed over.
    const goesOn = await v.api.group(
      'Goes on',
      async (api) => {
        await api.run('no.such-operation', {});
      },
      { onFailure: 'continue' },
    );
    expect(goesOn.ok).toBe(true);
  });

  it('reads Object members on the API as what they are, not as missing areas', async () => {
    const v = await create();
    expect(() => String(v.api)).not.toThrow();
    expect((v.api as unknown as { hasOwnProperty: unknown }).hasOwnProperty).toBeTypeOf('function');
  });

  it('runs viewer.run as one opaque step, its API calls inside it', async () => {
    const v = await create();
    await v.run('Load and preset', async (_plugin, api) => {
      await api.load.files({ files: [crambin()] });
      await api.structure.preset({ preset: 'atomic-detail' });
    });
    const steps = v.managers.history.steps.value;
    expect(steps.map((s) => s.label)).toEqual(['Load and preset']);
    expect(steps[0].runs.length).toBeGreaterThan(1);
  });

  it('stops a run through its signal and leaves nothing behind', async () => {
    const v = await create();
    gates.current = gate();
    const abort = new AbortController();
    const running = demoApi(v).slow({}, { signal: abort.signal });
    await new Promise((r) => setTimeout(r, 20));
    abort.abort();
    gates.current.open();
    // The work finished, but the call was stopped: its result says so and nothing is kept.
    const r = await running;
    expect(r.ok).toBe(true);
    const stopped = await demoApi(v).slow({}, { signal: abort.signal });
    expect(stopped).toMatchObject({ ok: false, code: 'aborted' });
  });
});

describe('areas', () => {
  it('describes operations, and watches their definitions', async () => {
    const v = await create();
    expect(v.api.describe('demo.slow')?.available).toBe(true);
    expect(v.api.describe('camera.reset')?.available).toMatch(/canvas/);
    expect(v.api.definitions.value.has('demo.slow')).toBe(true);
  });

  it('refuses an extension id taken by a core area, and an operation outside its area', async () => {
    await expect(
      createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), extensions: [{ id: 'canvas', title: 'X' }] } }),
    ).rejects.toThrow(/core API area/);
    class Wrong extends ApiArea {
      readonly id = 'wrong';
      readonly op = this.operation({ id: 'other.op', title: 'X', params: {}, run: async () => ApiResult.ok('x') });
    }
    await expect(
      createPrismViewerHeadless({
        spec: {
          ...DefaultPrismViewerSpec(),
          extensions: [{ id: 'wrong', title: 'W', api: { wrong: (c) => new Wrong(c) } }],
        },
      }),
    ).rejects.toThrow(/must start with its area's id/);
  });

  it("lets the viewer's features add operations to core areas, and no one else", async () => {
    const extra = {
      id: 'load.extra',
      title: 'Extra',
      params: { n: PD.Numeric(1) },
      run: async () => ApiResult.ok('Extra.', { nodes: [] }),
    };
    const feature: ViewerExtension = {
      id: 'loaders',
      title: 'Loaders',
      builtIn: true,
      contributes: () => [{ area: 'load', name: 'extra', operation: extra }],
    };
    const v = await create([feature]);
    expect((await v.api.run('load.extra', { n: 2 })).ok).toBe(true);
    expect(typeof (v.api.load as unknown as Record<string, unknown>).extra).toBe('function');
    await v.dispose();
    viewer = undefined;
    await expect(create([{ ...feature, builtIn: undefined }])).rejects.toThrow(/only the viewer's features/);
  });
});

describe('pointer events', () => {
  it('say what is under the pointer as a target', async () => {
    const v = await create();
    await v.api.load.files({ files: [cbs()] });
    const r = await (await import('../../sources/targets.js')).createTargetResolver(v.plugin).resolve('A:30');
    if (!r.ok) throw new Error(r.message);
    const events: ViewerPointerEvent[] = [];
    v.events.hover.subscribe((e) => events.push(e));
    v.plugin.behaviors.interaction.hover.next({
      current: { loci: r.parts[0].loci },
      buttons: 0,
      button: 0,
      modifiers: { shift: false, control: false, alt: false, meta: false },
    } as never);
    expect(events.at(-1)?.target).toBe('A:30');
    expect(events.at(-1)?.structure).toBe('1CBS');
  });
});

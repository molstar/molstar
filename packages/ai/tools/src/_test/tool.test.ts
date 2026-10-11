/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineTool, directRun, render } from '../index.js';
import { call, toy, tool } from './helpers.js';

const rename = tool({
  id: 'rename',
  description: 'Rename a scene.',
  input: {
    kind: 'object',
    fields: {
      scene: { kind: 'string', description: 'Scene id.' },
      title: { kind: 'string', description: 'Title.' },
      repr: { kind: 'string', enum: ['cartoon', 'ball-and-stick'], optional: true, description: 'Representation.' },
      size: { kind: 'number', optional: true },
      label: { kind: 'boolean', optional: true },
      frame: { kind: 'array', item: { kind: 'string' }, optional: true },
    },
  },
  modules: ['state', 'log'],
  run: ({ args, api }) => {
    // The arguments are typed as the shape says.
    const repr: 'cartoon' | 'ball-and-stick' | undefined = args.repr;
    const size: number | undefined = args.size;
    const label: boolean | undefined = args.label;
    const frame: string[] | undefined = args.frame;
    // @ts-expect-error the shape has no such field
    void args.colour;
    const scene: string = args.scene;
    // @ts-expect-error an optional field may be missing
    const sure: number = args.size;
    void [scene, sure];
    api.log.lines.push(JSON.stringify({ scene: args.scene, title: args.title, repr, size, label, frame }));
    if (args.title === 'same') return { ok: true, did: `Renamed ${args.scene}: already so.` };
    api.state.bump();
    return { ok: true, did: `Renamed ${args.scene}.`, changed: [`§4 ${args.scene}`] };
  },
});

describe('tools over field shapes (tools-and-session.md §1)', () => {
  it('normalizes what small models send before checking it', async () => {
    const api = toy();
    const r = await call(rename, api, {
      scene: ' pocket ',
      title: 'T',
      repr: 'Cartoon',
      size: '2',
      label: 'true',
      frame: 'A:315, site.pocket',
      note: '',
    });
    expect(r.status).toBe('succeeded');
    expect(JSON.parse(api.log.lines[0] ?? '')).toEqual({
      scene: 'pocket',
      title: 'T',
      repr: 'cartoon',
      size: 2,
      label: true,
      frame: ['A:315', 'site.pocket'],
    });
  });

  it('defines the JSON Schema the model sees, and refuses what it cannot use, saying why', async () => {
    expect(rename.definition).toMatchObject({
      tool: { id: 'rename' },
      name: 'rename',
      description: 'Rename a scene.',
      inputSchema: {
        type: 'object',
        properties: { scene: { type: 'string', description: 'Scene id.' } },
        required: ['scene', 'title'],
      },
    });
    const api = toy();
    expect(await call(rename, api, { scene: 'pocket' })).toMatchObject({
      status: 'failed',
      text: 'title: required (Title.).',
    });
    expect((await call(rename, api, { scene: 'pocket', title: 'x', repr: 'spheres' })).text).toBe(
      'repr: expected one of cartoon, ball-and-stick (got "spheres").',
    );
    expect((await call(rename, api, { scene: 'pocket', title: 'x', colour: 'red' })).text).toBe(
      'Renamed pocket.\nChanged: §4 pocket.\nWarnings (fix if it affects the request, else tell the user): Ignored colour: rename takes scene, title, repr, size, label, frame.',
    );
    expect(await call(rename, api, { scene: 'pocket', title: 'The pocket' })).toMatchObject({
      status: 'succeeded',
      text: 'Renamed pocket.\nChanged: §4 pocket.',
    });
    expect((await call(rename, api, { scene: 'pocket', title: 'same' })).text).toBe('Renamed pocket: already so.');
  });

  it('a prototype name is not an argument', async () => {
    const r = await call(rename, toy(), { scene: 'a', title: 'b', constructor: 'x', toString: 'y' });
    expect(r.text).toContain('Ignored constructor, toString');
  });

  it('gives a tool only the modules it declared, frozen', async () => {
    let seen: Record<string, unknown> = {};
    const t = tool({
      id: 'peek',
      description: 'Looks.',
      input: { kind: 'object', fields: {} },
      modules: ['state'],
      read: true,
      run: ({ api }) => {
        seen = api;
        return { ok: true, did: 'Looked.' };
      },
    });
    await call(t, toy(), {});
    expect(Object.keys(seen)).toEqual(['state']);
    expect(Object.isFrozen(seen)).toBe(true);
  });

  it('renders every result in one order', () => {
    expect(
      render({
        ok: true,
        did: 'Loaded 1IEP.',
        changed: ['struct.1iep', 'chain A'],
        view: 'chain A in focus',
        warnings: ['Two chains.', 'Waters left out.'],
        body: 'Chains: A, B.',
        next: 'Draw it with set_look.',
      }),
    ).toEqual({
      ok: true,
      text: 'Loaded 1IEP.\nChanged: struct.1iep, chain A.\nView: chain A in focus.\nWarnings (fix if it affects the request, else tell the user): Two chains. Waters left out.\nChains: A, B.\nNext: Draw it with set_look.',
    });
    expect(render({ ok: true, did: 'Selected A:315.', view: 'selected 1 residue.' })).toMatchObject({
      text: 'Selected A:315.\nView: selected 1 residue.',
    });
    // Nothing about a view when there is none.
    expect(render({ ok: true, did: 'Done.', view: '' })).toMatchObject({ text: 'Done.' });
  });

  it('a failure is the message and its hint, and says whether it may be retried', async () => {
    const t = tool({
      id: 'no',
      description: 'Fails.',
      input: { kind: 'object', fields: {} },
      modules: [],
      run: () => ({ ok: false, code: 'no.way', message: 'No.', hint: 'Try yes.', transient: true }),
    });
    expect(await call(t, toy(), {})).toMatchObject({
      status: 'failed',
      text: 'No. Try yes.',
      result: { transient: true },
    });
  });
});

describe('a long call', () => {
  it('goes out as progress frames before the result, in order', async () => {
    const slow = tool({
      id: 'slow',
      description: 'Takes a while.',
      input: { kind: 'object', fields: {} },
      modules: [],
      read: true,
      run: async ({ progress }) => {
        progress('submitted');
        await new Promise((r) => setTimeout(r, 5));
        progress('running, 1 s');
        await new Promise((r) => setTimeout(r, 5));
        return { ok: true, did: 'Done.' };
      },
    });
    const r = await call(slow, toy(), {});
    expect(r.frames.map((f) => (f.kind === 'progress' ? f.value : f.kind === 'result' ? f.status : f.kind))).toEqual([
      'submitted',
      'running, 1 s',
      'succeeded',
    ]);
  });

  it('a call that throws still ends with its failure, after what it said', async () => {
    const broken = tool({
      id: 'broken',
      description: 'Fails.',
      input: { kind: 'object', fields: {} },
      modules: [],
      read: true,
      run: async ({ progress }) => {
        progress('trying');
        throw new Error('It broke.');
      },
    });
    const r = await call(broken, toy(), {});
    expect(r.progress).toEqual(['trying']);
    expect(r).toMatchObject({ status: 'failed', text: 'It broke.' });
    expect(r.frames.at(-1)?.kind).toBe('result');
  });

  it('a throw that says it is transient is a transient failure', async () => {
    const flaky = tool({
      id: 'flaky',
      description: 'Fails.',
      input: { kind: 'object', fields: {} },
      modules: [],
      read: true,
      run: () => {
        throw Object.assign(new Error('Try later.'), { transient: true });
      },
    });
    expect((await call(flaky, toy(), {})).result).toMatchObject({ status: 'failed', transient: true });
  });
});

describe('a tool called from code (direct run)', () => {
  it('checks the arguments as a call is checked and returns the result as data', async () => {
    const api = toy();
    const run = directRun(rename);
    expect(run).toBeDefined();
    expect(await run?.({ scene: ' a ', title: 'b', size: '3' }, api, new AbortController().signal)).toEqual({
      ok: true,
      did: 'Renamed a.',
      changed: ['§4 a'],
    });
    expect(api.state.value).toBe(1);
    expect(await run?.({ scene: 'a' }, api, new AbortController().signal)).toMatchObject({
      ok: false,
      code: 'tool.arguments',
      message: 'title: required (Title.).',
    });
    // Arguments the tool does not take are refused here, not dropped.
    expect(await run?.({ scene: 'a', title: 'b', colour: 'red' }, api, new AbortController().signal)).toMatchObject({
      ok: false,
      message: expect.stringContaining('rename takes scene, title, repr, size, label, frame, not colour.'),
    });
  });

  it('is only for tools made by defineTool', () => {
    expect(directRun({ ...rename })).toBeUndefined();
  });

  it('defineTool without a host type works for tools that use no module', async () => {
    const t = defineTool({
      id: 'plain',
      description: 'Plain.',
      input: { kind: 'object', fields: { n: { kind: 'integer' } } },
      modules: [],
      run: ({ args }) => ({ ok: true, did: `n is ${args.n}.` }),
    });
    expect((await call(t, {}, { n: '4' })).text).toBe('n is 4.');
  });
});

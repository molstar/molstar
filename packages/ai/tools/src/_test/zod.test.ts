/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { z } from 'zod';
import { shapeInput } from '../index.js';
import { verbatim, zodInput } from '../zod.js';
import { call, toy, tool } from './helpers.js';

const input = z.object({
  scene: z.string().describe('Scene id.'),
  repr: z.enum(['cartoon', 'ball-and-stick']).optional().describe('Representation.'),
  size: z.number().optional().describe('Size factor.'),
  count: z.int().min(1).optional(),
  label: z.boolean().optional().describe('A label.'),
  frame: z.array(z.string()).optional().describe('Targets.'),
  residues: z.string().optional(),
});

describe('tools declared with zod (tools-and-session.md §1)', () => {
  const t = tool({
    id: 'zodded',
    description: 'Zod.',
    input: zodInput(input),
    modules: ['log'],
    run: ({ args, api }) => {
      // z.output
      const size: number | undefined = args.size;
      api.log.lines.push(JSON.stringify({ ...args, size }));
      return { ok: true, did: 'Done.' };
    },
  });

  it('gives the JSON Schema zod makes, input side, draft 2020-12', () => {
    expect(t.definition.inputSchema).toEqual({
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      type: 'object',
      properties: {
        scene: { type: 'string', description: 'Scene id.' },
        repr: { type: 'string', enum: ['cartoon', 'ball-and-stick'], description: 'Representation.' },
        size: { type: 'number', description: 'Size factor.' },
        count: { type: 'integer', minimum: 1 },
        label: { type: 'boolean', description: 'A label.' },
        frame: { type: 'array', items: { type: 'string' }, description: 'Targets.' },
        residues: { type: 'string' },
      },
      required: ['scene'],
    });
  });

  it('normalizes what small models send, from the shape read off the zod object', async () => {
    const api = toy();
    const r = await call(t, api, {
      scene: ' pocket ',
      repr: 'Cartoon',
      size: '2',
      count: '3',
      label: 'true',
      frame: 'A:315, site.pocket',
      residues: ['A:315', 'A:318'],
      note: '',
    });
    expect(r.status).toBe('succeeded');
    expect(JSON.parse(api.log.lines[0] ?? '')).toEqual({
      scene: 'pocket',
      repr: 'cartoon',
      size: 2,
      count: 3,
      label: true,
      frame: ['A:315', 'site.pocket'],
      residues: 'A:315, A:318',
    });
  });

  it('refuses with one message: the path, what it expects, what came, the description', async () => {
    expect((await call(t, toy(), {})).text).toBe('scene: required (Scene id.).');
    expect((await call(t, toy(), { scene: 's', repr: 'spheres' })).text).toBe(
      'repr: expected one of cartoon, ball-and-stick (got "spheres").',
    );
    expect((await call(t, toy(), { scene: 's', size: 'big' })).text).toBe(
      'size: expected number (got "big") (Size factor.).',
    );
  });

  it('verbatim marks a field taken as sent, wrapped or not', async () => {
    const code = z.object({
      find: verbatim(z.string()),
      replace: verbatim(z.string()),
      note: verbatim(z.string()).optional(),
      keep: z.string().optional(),
    });
    const shape = zodInput(code).shape;
    expect(shape.fields.find).toEqual({ kind: 'string', verbatim: true });
    const c = tool({
      id: 'code',
      description: 'Code.',
      input: zodInput(code),
      modules: ['log'],
      run: ({ args, api }) => {
        api.log.lines.push(JSON.stringify(args));
        return { ok: true, did: 'Done.' };
      },
    });
    const api = toy();
    await call(c, api, { find: '  x();\n', replace: '', keep: ' y ' });
    expect(JSON.parse(api.log.lines[0] ?? '')).toEqual({ find: '  x();\n', replace: '', keep: 'y' });
  });

  it('is read as the same field shape a host would declare', () => {
    const shape = zodInput(
      z.object({
        a: z.string().describe('A.'),
        b: z.number().min(1).max(5).optional(),
        c: z.int().default(2),
        d: z.array(z.object({ x: z.boolean() })),
        e: z.discriminatedUnion('type', [
          z.object({ type: z.literal('box'), size: z.number() }),
          z.object({ type: z.literal('tag'), text: z.string() }),
        ]),
        f: z.literal('on'),
        g: z.any(),
      }),
    ).shape;
    expect(shape).toEqual({
      kind: 'object',
      fields: {
        a: { kind: 'string', description: 'A.' },
        b: { kind: 'number', min: 1, max: 5, optional: true },
        c: { kind: 'integer', default: 2 },
        d: { kind: 'array', item: { kind: 'object', fields: { x: { kind: 'boolean' } } } },
        e: {
          kind: 'variant',
          tag: 'type',
          variants: {
            box: { kind: 'object', fields: { size: { kind: 'number' } } },
            tag: { kind: 'object', fields: { text: { kind: 'string' } } },
          },
        },
        f: { kind: 'string', enum: ['on'] },
        g: { kind: 'json' },
      },
    });
  });

  it('a variant is normalized by the variant its tag names', async () => {
    const v = tool({
      id: 'variant',
      description: 'Variant.',
      input: zodInput(
        z.object({
          mark: z.discriminatedUnion('type', [
            z.object({ type: z.literal('box'), size: z.number() }),
            z.object({ type: z.literal('tag'), text: z.string() }),
          ]),
        }),
      ),
      modules: ['log'],
      run: ({ args, api }) => {
        api.log.lines.push(JSON.stringify(args.mark));
        return { ok: true, did: 'Done.' };
      },
    });
    const api = toy();
    expect((await call(v, api, { mark: { type: 'Box', size: '2' } })).status).toBe('succeeded');
    expect(api.log.lines).toEqual(['{"type":"box","size":2}']);
    expect((await call(v, api, { mark: { type: 'cone' } })).text).toBe(
      'mark.type: expected one of box, tag (got "cone").',
    );
  });

  it('keeps what a field shape cannot say: zod checks the original', async () => {
    const e = tool({
      id: 'email',
      description: 'Email.',
      input: zodInput(z.object({ to: z.string().regex(/^\S+@\S+$/, 'not an address') })),
      modules: [],
      run: () => ({ ok: true, did: 'Sent.' }),
    });
    expect((await call(e, toy(), { to: 'nobody' })).text).toBe('to: not an address.');
    expect((await call(e, toy(), { to: 'a@b' })).status).toBe('succeeded');
  });

  it('a field shape and the zod object it equals give the same JSON Schema', () => {
    const fromShape = shapeInput({
      kind: 'object',
      fields: {
        scene: { kind: 'string', description: 'Scene id.' },
        repr: { kind: 'string', enum: ['cartoon', 'ball-and-stick'], optional: true, description: 'Representation.' },
        count: { kind: 'integer', min: 1, optional: true },
      },
    });
    const fromZod = zodInput(
      z.object({
        scene: z.string().describe('Scene id.'),
        repr: z.enum(['cartoon', 'ball-and-stick']).optional().describe('Representation.'),
        count: z.int().min(1).optional(),
      }),
    );
    expect(fromShape.jsonSchema).toEqual(fromZod.jsonSchema);
  });
});

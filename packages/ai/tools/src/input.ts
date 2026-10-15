/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject } from '@molstar/ai-agent';
import { z } from 'zod';
import type { ArgsOf, FieldShape, ObjectShape } from './shape.js';

/**
 * A tool's input as `ai-tools` holds it (tools-and-session.md §1): the shape its arguments are normalized against, the
 * JSON Schema the model sees, and the check. Made once, from a field shape (`shapeInput`) or a zod object
 * (`zodInput`, in `@molstar/ai-tools/zod`); zod is this package's own business, so no type here names it.
 */
export interface ToolInput<A = unknown> {
  readonly shape: ObjectShape;
  /** JSON Schema draft 2020-12, the input side. */
  readonly jsonSchema: JsonObject;
  /** The arguments (normalized first: `normalize`) checked: the typed value, or one message the model can act on. */
  parse(args: Readonly<Record<string, unknown>>): { ok: true; value: A } | { ok: false; message: string };
}

/** What a tool's `input` may be: a field shape, or an input made elsewhere (zod). */
export type InputSpec = ObjectShape | ToolInput;

/** The arguments an input checks to. */
export type ArgsFrom<I extends InputSpec> = I extends ObjectShape
  ? ArgsOf<I>
  : I extends ToolInput<infer A>
    ? A
    : never;

/** The input of a spec: a shape is turned into its input; an input is itself. */
export const inputOf = (spec: InputSpec): ToolInput => ('fields' in spec ? shapeInput(spec) : spec);

// ------------------------------------------------------------------ field shape -> zod

const clip = (s: string) => (s.length > 100 ? `${s.slice(0, 99)}…` : s);

function fieldSchema(f: FieldShape): z.ZodType {
  let s: z.ZodType;
  switch (f.kind) {
    case 'string':
      s = f.enum?.length ? z.enum(f.enum as [string, ...string[]]) : z.string();
      break;
    case 'number':
    case 'integer': {
      let n = f.kind === 'integer' ? z.number().int() : z.number();
      if (f.min !== undefined) n = n.min(f.min);
      if (f.max !== undefined) n = n.max(f.max);
      s = n;
      break;
    }
    case 'boolean':
      s = z.boolean();
      break;
    case 'array': {
      const a = z.array(fieldSchema(f.item));
      s = f.max !== undefined ? a.max(f.max) : a;
      break;
    }
    case 'object':
      s = objectSchema(f);
      break;
    case 'variant': {
      const options = Object.entries(f.variants).map(([key, v]) =>
        z.object({ [f.tag]: z.literal(key), ...objectSchema(v).shape }),
      );
      if (!options.length) throw new Error(`The variant on "${f.tag}" has no variants.`);
      s = z.discriminatedUnion(f.tag, options as [(typeof options)[0], ...typeof options]);
      break;
    }
    case 'json':
      s = z.any();
      break;
  }
  if (f.default !== undefined) s = s.default(f.default as never);
  else if (f.optional) s = s.optional();
  const note = [f.description, f.kind === 'string' && f.format ? `Format: ${f.format}.` : undefined]
    .filter(Boolean)
    .join(' ');
  return note ? s.describe(note) : s;
}

function objectSchema(shape: ObjectShape): z.ZodObject {
  return z.object(Object.fromEntries(Object.entries(shape.fields).map(([k, f]) => [k, fieldSchema(f)])));
}

// ------------------------------------------------------------------ the check

/** The value a failing path points at in what was checked. */
function valueAt(data: unknown, path: readonly PropertyKey[]): unknown {
  let cur = data;
  for (const p of path) {
    if (typeof cur !== 'object' || cur === null) return undefined;
    cur = (cur as Record<PropertyKey, unknown>)[p];
  }
  return cur;
}

/** A failed parse as one message the model can act on: the field, what it expects, what came, its description. */
function argumentsMessage(shape: ObjectShape, error: z.ZodError, data: unknown): string {
  const issue = error.issues[0];
  const path = issue?.path ?? [];
  // The whole path, as the model wrote it: claims[2].page.
  const field = path.map((p, i) => (typeof p === 'number' ? `[${p}]` : i ? `.${String(p)}` : String(p))).join('');
  const describe = path.length === 1 ? shape.fields[String(path[0])]?.description : undefined;
  const about = describe ? ` (${describe})` : '';
  const given = valueAt(data, path);
  const got = given === undefined ? '' : ` (got ${clip(JSON.stringify(given as Json))})`;
  const choices =
    issue?.code === 'invalid_value'
      ? issue.values
      : issue?.code === 'invalid_union' && 'options' in issue
        ? (issue.options as unknown[])
        : undefined;
  if (choices && given === undefined) return `${field}: required, one of ${choices.join(', ')}${about}.`;
  if (choices) return `${field}: expected one of ${choices.join(', ')}${got}.`;
  if (given === undefined && (issue?.code === 'invalid_type' || issue?.code === 'custom'))
    return `${field}: required${about}.`;
  if (issue?.code === 'invalid_type') return `${field}: expected ${issue.expected}${got}${about}.`;
  return `${field ? `${field}: ` : ''}${issue?.message ?? 'invalid arguments'}${about}.`;
}

const JSON_SCHEMA = {
  target: 'draft-2020-12',
  io: 'input',
  unrepresentable: 'throw',
  cycles: 'throw',
  reused: 'inline',
  // An integer carries the safe-integer range as bounds: noise for the model, so only the bounds that were asked for.
  override: ({ jsonSchema }: { jsonSchema: Record<string, unknown> }) => {
    if (jsonSchema.type !== 'integer') return;
    if (jsonSchema.maximum === Number.MAX_SAFE_INTEGER) delete jsonSchema.maximum;
    if (jsonSchema.minimum === Number.MIN_SAFE_INTEGER) delete jsonSchema.minimum;
  },
} as const;

/**
 * An input over a zod object and the field shape that reads it: the JSON Schema the model sees (draft 2020-12, the
 * input side, nothing a schema cannot say faithfully), and the check.
 */
export function zodBacked<A>(shape: ObjectShape, schema: z.ZodObject): ToolInput<A> {
  // The dialect is the providers' own: naming it costs every tool tokens on every request.
  const { $schema: _dialect, ...jsonSchema } = z.toJSONSchema(schema, JSON_SCHEMA) as JsonObject;
  return {
    shape,
    jsonSchema,
    parse(args) {
      const parsed = schema.safeParse(args);
      return parsed.success
        ? { ok: true, value: parsed.data as A }
        : { ok: false, message: argumentsMessage(shape, parsed.error, args) };
    },
  };
}

/** An input from a field shape (tools-and-session.md §2), validated and described through zod. */
export function shapeInput<const S extends ObjectShape>(shape: S): ToolInput<ArgsOf<S>> {
  return zodBacked<ArgsOf<S>>(shape, objectSchema(shape));
}

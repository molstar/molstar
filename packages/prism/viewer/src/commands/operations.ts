/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandDefinition, CommandValues, IconName, KeyChord } from '@molstar/prism-ui/core';
import type { OperationFn, ViewerApi } from '../api/area.js';
import { filePickers, filesOf } from '../api/files.js';
import type { ViewerCommandContext } from '../api/managers.js';

// Commands over operations (operations.md §3): one command per operation by default, with the operation's id, title
// and params (resolved as the API describes them, so registry-driven params are current); its form is the params'
// form, with file pickers for file values. It runs the operation through the API bound to the person, says why it
// cannot run when the operation says so, and its result line is what the status column says.

export interface OperationCommandOptions {
  readonly category: string;
  readonly icon?: IconName;
  readonly alias?: string | readonly string[];
  readonly shortcut?: KeyChord;
  /** The apply button's label ("Set", "Open"). */
  readonly action?: string;
  readonly group?: string;
  /** Another id than the operation's (a generated command with a value fixed). */
  readonly id?: string;
  readonly title?: string;
  /** Values fixed for this command; the form shows the rest. */
  readonly values?: Readonly<Record<string, unknown>>;
}

/** Runs an operation by id as the person; throws the failure's message and hint (the command menu shows it). */
export async function runAsCommand(
  api: ViewerApi,
  id: string,
  values: Readonly<Record<string, unknown>>,
): Promise<{ message: string; ops: string }> {
  const r = await api.run(id, values);
  if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
  return { message: r.did, ops: `${id}(${formatValues(values)})` };
}

export function operationCommand(
  api: ViewerApi,
  op: OperationFn<any, any>,
  o: OperationCommandOptions,
): CommandDefinition<ViewerCommandContext> {
  const def = op.definition;
  const fixed = o.values ?? {};
  const form = () => {
    const params = api.describe(def.id)?.params ?? {};
    const open = Object.keys(params).filter((k) => !(k in fixed));
    return open.length ? formParams(Object.fromEntries(open.map((k) => [k, params[k]]))) : undefined;
  };
  // Target fields read and complete as targets (`sources/target-format.ts`).
  const fields = Object.fromEntries(
    Object.entries(api.describe(def.id)?.params ?? {})
      .filter(([k, p]) => p.type === 'value' && typeof p.defaultValue === 'string' && !(k in fixed))
      .map(([k]) => [k, { format: 'target' }]),
  );
  const files = new Set(
    Object.entries(api.describe(def.id)?.params ?? {})
      .filter(([, p]) => p.type === 'value' && (Array.isArray(p.defaultValue) || p.defaultValue === null))
      .map(([k]) => k),
  );
  return {
    id: o.id ?? def.id,
    title: o.title ?? def.title,
    category: o.category,
    ...(o.icon ? { icon: o.icon } : {}),
    ...(o.alias ? { alias: o.alias } : {}),
    ...(o.shortcut ? { shortcut: o.shortcut } : {}),
    ...(o.action ? { action: o.action } : {}),
    ...(o.group ? { group: o.group } : {}),
    ...(form() ? { params: () => form() ?? {}, defaults: () => formDefaults(form()) } : {}),
    ...(Object.keys(fields).length ? { fields } : {}),
    when: () => api.describe(def.id)?.available ?? `There is no operation '${def.id}'.`,
    async run(args) {
      // File pickers give Mol*'s file assets; the operation takes DOM files.
      const values: Record<string, unknown> = { ...args, ...fixed };
      for (const k of files) {
        if (!(k in values)) continue;
        const list = filesOf(values[k]);
        values[k] = Array.isArray(api.describe(def.id)?.params[k]?.defaultValue) ? list : (list[0] ?? null);
      }
      const r = await runAsCommand(api, def.id, values);
      return { ...r, undoable: !def.read && def.history !== false && def.undoable !== false };
    },
  };
}

/**
 * A form for an operation's params: file pickers for files, text fields for targets (asked for). A text field whose
 * default is empty means "the default" (the only structure, …) when left empty, so it is optional.
 */
function formParams(params: PD.Params): PD.Params {
  const out = filePickers(params);
  for (const [k, p] of Object.entries(out)) {
    if (p.type === 'value' && typeof p.defaultValue === 'string') {
      out[k] = PD.Text(p.defaultValue, { label: p.label, description: p.description });
    } else if (p.type === 'text' && p.defaultValue === '' && !p.isOptional) {
      out[k] = { ...p, isOptional: true };
    }
  }
  return out;
}

/** The values a form starts with: its choices' defaults (selects, switches, numbers); texts and files start empty. */
function formDefaults(params: PD.Params | undefined): CommandValues {
  const out: Record<string, string | number | boolean> = {};
  for (const [k, p] of Object.entries(params ?? {})) {
    if (p.type === 'select' || p.type === 'boolean' || p.type === 'number') out[k] = p.defaultValue;
  }
  return out;
}

function formatValues(values: Readonly<Record<string, unknown>>): string {
  const parts: string[] = [];
  for (const [k, v] of Object.entries(values)) {
    if (v === undefined || typeof v === 'function' || (typeof File !== 'undefined' && v instanceof File)) continue;
    if (Array.isArray(v) && v.some((x) => typeof File !== 'undefined' && x instanceof File)) continue;
    parts.push(`${k}: ${JSON.stringify(v)}`);
  }
  return parts.length ? `{ ${parts.join(', ')} }` : '';
}

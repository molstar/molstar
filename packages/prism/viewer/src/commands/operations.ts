/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, IconName, KeyChord } from '@molstar/prism-ui';
import type { OperationFn, ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';

// Commands over operations (operations.md §3): one command per operation by default, with the operation's id, title
// and params; its form is the params' form. It runs the operation through the API bound to the person, and its result
// line is what the status column says. A failure is the command's failure, with the hint.

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
): Promise<{ message: string; undoable?: boolean; ops: string }> {
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
  const params = typeof def.params === 'function' ? undefined : def.params;
  const open = params ? Object.keys(params).filter((k) => !(k in fixed)) : [];
  return {
    id: o.id ?? def.id,
    title: o.title ?? def.title,
    category: o.category,
    ...(o.icon ? { icon: o.icon } : {}),
    ...(o.alias ? { alias: o.alias } : {}),
    ...(o.shortcut ? { shortcut: o.shortcut } : {}),
    ...(o.action ? { action: o.action } : {}),
    ...(o.group ? { group: o.group } : {}),
    ...(params && open.length ? { params: () => Object.fromEntries(open.map((k) => [k, params[k]])) } : {}),
    async run(args) {
      const r = await runAsCommand(api, def.id, { ...args, ...fixed });
      return { ...r, undoable: !def.read && def.history !== false };
    },
  };
}

function formatValues(values: Readonly<Record<string, unknown>>): string {
  const parts: string[] = [];
  for (const [k, v] of Object.entries(values)) {
    if (
      v === undefined ||
      typeof v === 'function' ||
      (typeof v === 'object' && v !== null && !Array.isArray(v) && 'name' in v && 'size' in v)
    )
      continue;
    parts.push(`${k}: ${JSON.stringify(v)}`);
  }
  return parts.length ? `{ ${parts.join(', ')} }` : '';
}

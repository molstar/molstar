/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { readSettings, SettingsMapping, type SettingsGroupId, type SettingsGroupState } from '../sources/settings.js';
import { type ApiContext, ApiArea, ApiResult } from './area.js';

// `settings` (operations.md §2, parity.md §5, plan 3.3): the Settings tab's groups over the plugin's real params
// (Canvas3D, its context, the interactivity manager, the behaviors) and the viewer's own (Undo). `settings.get` reads
// them; `settings.set` changes one group's values (partial values are merged into the group's current ones) as a
// history step: its checkpoint holds the canvas props and the behaviors, so Undo restores them. The canvas context's
// props (Performance) and the picking level are not in the checkpoint (sequence.md §7.2): the result says so.

/** What `settings.get` gives back for one group. */
export interface SettingsInfo {
  readonly id: SettingsGroupId;
  readonly title: string;
  readonly description: string;
  readonly values: Record<string, unknown>;
}

const isPlain = (x: unknown): x is Record<string, any> =>
  !!x && typeof x === 'object' && !Array.isArray(x) && Object.getPrototypeOf(x) === Object.prototype;

/** A patch checked against the params it is for (deep: groups and mapped params), by the path of the value. */
function check(params: PD.Params, patch: Record<string, unknown>, path: string): string | undefined {
  for (const [key, value] of Object.entries(patch)) {
    const at = path ? `${path}.${key}` : key;
    const p = params[key];
    if (!p) return `Unknown setting '${at}'. It takes: ${Object.keys(params).join(', ')}.`;
    const problem = checkParam(p, value, at);
    if (problem) return problem;
  }
  return undefined;
}

function checkParam(p: PD.Any, value: unknown, at: string): string | undefined {
  switch (p.type) {
    case 'number':
      if (typeof value !== 'number' || Number.isNaN(value)) return `'${at}' takes a number.`;
      if ((p.min !== undefined && value < p.min) || (p.max !== undefined && value > p.max)) {
        return `'${at}' takes a number from ${p.min ?? '…'} to ${p.max ?? '…'}.`;
      }
      return undefined;
    case 'boolean':
      return typeof value === 'boolean' ? undefined : `'${at}' takes true or false.`;
    case 'text':
      return typeof value === 'string' ? undefined : `'${at}' takes text.`;
    case 'color':
      return typeof value === 'number' ? undefined : `'${at}' takes a color as a number (0xrrggbb).`;
    case 'select': {
      const options = p.options.map((o) => o[0]);
      return options.includes(value as never) ? undefined : `'${at}' takes: ${options.join(', ')}.`;
    }
    case 'group':
      return isPlain(value) ? check(p.params, value, at) : `'${at}' takes an object.`;
    case 'mapped': {
      if (!isPlain(value)) return `'${at}' takes { name, params }.`;
      const options = p.select.options.map((o) => o[0]);
      if (value.name !== undefined && !options.includes(value.name))
        return `'${at}.name' takes: ${options.join(', ')}.`;
      if (value.params === undefined) return undefined;
      const inner = p.map(value.name);
      return inner.type === 'group' && isPlain(value.params)
        ? check(inner.params, value.params, `${at}.params`)
        : undefined;
    }
    default:
      return undefined;
  }
}

/** A mapped value patched: another name takes that option's defaults, the same one keeps its params. */
function mergeMapped(p: PD.Mapped<any>, current: any, patch: Record<string, any>) {
  const name = patch.name ?? current?.name;
  const same = current?.name === name;
  const params = same ? current.params : structuredClone(p.map(name).defaultValue);
  if (patch.params === undefined) return { name, params };
  return { name, params: isPlain(params) && isPlain(patch.params) ? { ...params, ...patch.params } : patch.params };
}

/** A group's values with a (checked) patch merged in: groups one by one, anything else replaced. */
function mergeValues(params: PD.Params, current: Record<string, any>, patch: Record<string, any>): Record<string, any> {
  const out: Record<string, any> = { ...current };
  for (const [key, value] of Object.entries(patch)) {
    const p = params[key];
    if (p.type === 'group' && isPlain(value)) out[key] = mergeValues(p.params, current[key] ?? {}, value);
    else if (p.type === 'mapped' && isPlain(value)) out[key] = mergeMapped(p, current[key], value);
    else out[key] = value;
  }
  return out;
}

/** The top-level params whose value differs. */
function changedKeys(group: SettingsGroupState, next: Record<string, unknown>): string[] {
  return Object.keys(group.params).filter((k) => !PD.isParamEqual(group.params[k], group.values[k], next[k]));
}

const subset = (next: Record<string, unknown>, keys: readonly string[], among: readonly string[]) =>
  Object.fromEntries(keys.filter((k) => among.includes(k)).map((k) => [k, next[k]]));

export class SettingsApi extends ApiArea {
  readonly id = 'settings';

  readonly get = this.operation({
    id: 'settings.get',
    title: 'Read the settings',
    params: { group: PD.Text('', { description: 'A group: appearance, camera, …; empty: every group' }) },
    read: true,
    async run({ group }, ctx) {
      const groups = readSettings(ctx.plugin, ctx.managers.history.checkpointing.value);
      const wanted = group
        ? groups.filter((g) => g.id === group || g.title.toLowerCase() === group.toLowerCase())
        : groups;
      if (group && !wanted.length) {
        return ApiResult.fail(
          groups.some((g) => g.id === group) ? 'failed' : 'not-found',
          `There is no settings group '${group}' here.`,
          `Groups: ${groups.map((g) => g.id).join(', ')}.`,
        );
      }
      const value: SettingsInfo[] = wanted.map(({ id, title, description, values }) => ({
        id,
        title,
        description,
        values,
      }));
      return ApiResult.ok(
        `Settings: ${value.map((g) => `${g.title} (${g.id}): ${Object.keys(g.values).join(', ')}`).join('; ')}.`,
        {},
        { value },
      );
    },
  });

  readonly set = this.operation({
    id: 'settings.set',
    title: 'Change settings',
    params: {
      group: PD.Select<SettingsGroupId>(
        'appearance',
        SettingsMapping.map((g) => [g.id, g.title] as [SettingsGroupId, string]),
        { description: 'The group of settings' },
      ),
      values: PD.Value<Record<string, unknown>>(
        {},
        { description: "The group's values to change, as `settings.get` gives them (some of them: the rest stay)" },
      ),
    },
    // `modules` is per operation, and the Behaviors and Viewer groups work without a canvas: the canvas groups check.
    async run({ group, values }, ctx) {
      const groups = readSettings(ctx.plugin, ctx.managers.history.checkpointing.value);
      const current = groups.find((g) => g.id === group);
      const title = SettingsMapping.find((g) => g.id === group)?.title ?? group;
      if (!current) {
        const canvas = group !== 'behaviors' && group !== 'viewer';
        return canvas
          ? ApiResult.fail(
              'unavailable',
              `${title} settings need a canvas.`,
              ctx.plugin.canvas3d ? 'This host has no canvas context.' : 'This host has no canvas.',
            )
          : ApiResult.fail('not-found', `No behavior has settings to change.`);
      }
      const problem = check(current.params, values, '');
      if (problem) return ApiResult.fail('invalid', problem, `Read the values with settings.get.`);
      const next = mergeValues(current.params, current.values, values);
      const changed = changedKeys(current, next);
      if (!changed.length) return ApiResult.already(`The ${title} settings are already so.`);
      return apply(ctx, current, next, changed, title);
    },
  });
}

async function apply(
  ctx: ApiContext,
  group: SettingsGroupState,
  next: Record<string, unknown>,
  changed: string[],
  title: string,
): Promise<ApiResult> {
  const { plugin } = ctx;
  const warnings: string[] = [];
  const nodes: string[] = [];
  let view = false;
  const canvas = subset(next, changed, group.canvas);
  if (Object.keys(canvas).length) {
    plugin.setCanvas3DProps(canvas as never);
    view = true;
  }
  const context = subset(next, changed, group.context);
  if (Object.keys(context).length) {
    plugin.setCanvas3DContextProps(context as never);
    view = true;
    warnings.push('The canvas context settings (resolution, transparency) are not restored by Undo.');
  }
  const interactivity = subset(next, changed, group.interactivity);
  if (Object.keys(interactivity).length) {
    plugin.managers.interactivity.setProps(interactivity as never);
    warnings.push('The picking level is not restored by Undo.');
  }
  if (group.behaviors) {
    const state = plugin.state.behaviors;
    const tree = state.build();
    const refs: string[] = [];
    for (const key of changed) {
      const ref = group.behaviors[key];
      const cell = ref ? state.cells.get(ref) : undefined;
      if (!ref || !cell) continue;
      tree.to(ref).update({ ...(cell.transform.params as object), ...(next[key] as object) } as never);
      refs.push(ref);
    }
    await plugin.runTask(state.updateTree(tree));
    const failed = refs.map((r) => state.cells.get(r)).find((c) => c?.status === 'error');
    if (failed)
      return ApiResult.fail('failed', `A behavior could not take the change: ${failed.errorText ?? 'error'}.`);
    nodes.push(...refs);
  }
  if (group.viewer) ctx.managers.history.setCheckpointing(next.checkpointing as boolean);
  return ApiResult.ok(
    `Changed the ${title} settings: ${changed.join(', ')}.`,
    { ...(view ? { view } : {}), ...(nodes.length ? { nodes } : {}) },
    warnings.length ? { warnings } : {},
  );
}

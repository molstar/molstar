/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import {
  animateKinemage,
  applyViewSnapshot,
  getLoadedKinemages,
  isKinemageVisible,
  type KinemageVisibilityTarget,
  type LoadedKinemage,
  setKinemageVisibility,
} from '@molstar/kinemage-extension/model';
import type { Kinemage } from '@molstar/kinemage-extension/reader/schema';
import { type ApiAreaContext, ApiArea, ApiResult, type ApiFailure } from '@molstar/prism-viewer/api/area';

// The `kinemage` area (parity.md §6, README decision 27): what the v5 right-panel controls did, as operations. Which
// groups, subgroups and masters are shown, the views and the animate frames are the kinemage's own state in the plugin
// (its visibility controller), so showing and hiding is a step the history restores; a view only moves the camera.

export interface KinemageGroupInfo {
  readonly name: string;
  readonly visible: boolean;
  /** Marked `nobutton`: the v5 panel does not list it; it can still be shown and hidden. */
  readonly nobutton: boolean;
  /** In the animate (1) or 2animate (2) frames. */
  readonly animate?: 1 | 2;
  readonly subgroups: readonly { readonly name: string; readonly visible: boolean; readonly nobutton: boolean }[];
}

export interface KinemageInfo {
  /** 1-based. */
  readonly index: number;
  readonly title: string;
  readonly ref: string;
  readonly views: readonly { readonly key: string; readonly name: string }[];
  readonly groups: readonly KinemageGroupInfo[];
  readonly masters: readonly { readonly name: string; readonly visible: boolean }[];
  readonly animate: readonly string[];
  readonly animate2: readonly string[];
}

const KinemageParam = () =>
  PD.Text('', { description: 'A kinemage by title, number (from 1) or state ref; empty: the only one' });

const Kinds: [KinemageVisibilityTarget['type'], string][] = [
  ['group', 'Group'],
  ['subgroup', 'Subgroup'],
  ['master', 'Master'],
];

const titleOf = (k: Kinemage) =>
  k.pdbfile || (k.caption ? k.caption.trim().split('\n')[0].slice(0, 60) : '') || 'Kinemage';

const ifLoaded = (ctx: ApiAreaContext): true | string =>
  getLoadedKinemages(ctx.plugin).length ? true : 'No kinemage is loaded.';

/** The kinemage a value names: by state ref, number, then title (any case); empty is the only one. */
function resolve(ctx: ApiAreaContext, ref: string): (LoadedKinemage & { index: number }) | ApiFailure {
  const all = getLoadedKinemages(ctx.plugin).map((k, i) => ({ ...k, index: i + 1 }));
  if (!all.length) return ApiResult.fail('not-found', 'No kinemage is loaded.', 'Open a .kin file (load.files).');
  const text = ref.trim();
  if (!text) {
    if (all.length === 1) return all[0];
    return ApiResult.fail(
      'ambiguous',
      `${all.length} kinemages are loaded.`,
      `Say which: ${all.map((k) => `${k.index}. ${titleOf(k.kinData)}`).join('; ')}.`,
    );
  }
  const found =
    all.find((k) => k.ref === text) ??
    all.find((k) => String(k.index) === text) ??
    all.filter((k) => titleOf(k.kinData).toLowerCase() === text.toLowerCase());
  if (Array.isArray(found)) {
    if (found.length === 1) return found[0];
    if (found.length > 1) {
      return ApiResult.fail('ambiguous', `${found.length} kinemages are titled '${text}'.`, 'Give its number.');
    }
    return ApiResult.fail('not-found', `There is no kinemage '${text}'.`, 'List them with kinemage.get.');
  }
  return found;
}

function describe(ctx: ApiAreaContext, k: LoadedKinemage, index: number): KinemageInfo {
  const d = k.kinData;
  const shown = (type: KinemageVisibilityTarget['type'], key: string) =>
    isKinemageVisible(ctx.plugin, k.ref, { type, key });
  const subgroups = Object.entries(d.subgroupDict ?? {});
  return {
    index,
    title: titleOf(d),
    ref: k.ref,
    views: Object.entries(d.viewDict ?? {}).map(([key, v]) => ({ key, name: v.name || `View ${key}` })),
    groups: Object.entries(d.groupDict ?? {}).map(([name, info]) => ({
      name,
      visible: shown('group', name),
      nobutton: !!info.nobutton,
      ...(d.groupsAnimate?.includes(name)
        ? { animate: 1 as const }
        : d.groupsAnimate2?.includes(name)
          ? { animate: 2 as const }
          : {}),
      subgroups: subgroups
        .filter(([key]) => key.startsWith(`${name}:`))
        .map(([key, s]) => ({
          name: key.slice(name.length + 1),
          visible: shown('subgroup', key),
          nobutton: !!s.nobutton,
        })),
    })),
    masters: Object.keys(d.masterDict ?? {}).map((name) => ({ name, visible: shown('master', name) })),
    animate: d.groupsAnimate ?? [],
    animate2: d.groupsAnimate2 ?? [],
  };
}

/** An item by exact key, then (subgroups) by the name after its group, then any case; ambiguity is reported. */
function findKey(keys: string[], name: string, kind: KinemageVisibilityTarget['type']): string | ApiFailure {
  const text = name.trim();
  if (!text) return ApiResult.fail('invalid', `Which ${kind}?`, 'Give its name (kinemage.get lists them).');
  if (keys.includes(text)) return text;
  const lower = text.toLowerCase();
  const tail = (k: string) => (kind === 'subgroup' && k.includes(':') ? k.slice(k.indexOf(':') + 1) : k);
  const matches = keys.filter((k) => k.toLowerCase() === lower || tail(k).toLowerCase() === lower);
  if (matches.length === 1) return matches[0];
  if (matches.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `${matches.length} ${kind}s are named '${text}'.`,
      `Use: ${matches.join('; ')}.`,
    );
  }
  return ApiResult.fail('not-found', `There is no ${kind} '${text}'.`, 'List them with kinemage.get.');
}

export class KinemageApi extends ApiArea {
  readonly id = 'kinemage';

  readonly get = this.operation({
    id: 'kinemage.get',
    title: 'List the kinemages',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = getLoadedKinemages(ctx.plugin).map((k, i) => describe(ctx, k, i + 1));
      if (!list.length) return ApiResult.ok('No kinemage is loaded.', {}, { value: { kinemages: list } });
      const lines = list.flatMap((k) => [
        `${k.index}. ${k.title}`,
        ...(k.views.length ? [`   views: ${k.views.map((v) => v.name).join(', ')}`] : []),
        ...k.groups.flatMap((g) => [
          `   group ${g.name}: ${g.visible ? 'shown' : 'hidden'}${g.animate ? ` (animate${g.animate === 2 ? '2' : ''})` : ''}`,
          ...g.subgroups.map((s) => `     subgroup ${s.name}: ${s.visible ? 'shown' : 'hidden'}`),
        ]),
        ...k.masters.map((m) => `   master ${m.name}: ${m.visible ? 'shown' : 'hidden'}`),
      ]);
      return ApiResult.ok(lines.join('\n'), {}, { value: { kinemages: list } });
    },
  });

  readonly visibility = this.operation({
    id: 'kinemage.visibility',
    title: 'Show or hide a kinemage group, subgroup or master',
    params: {
      kind: PD.Select<KinemageVisibilityTarget['type']>('group', Kinds, { description: 'What to show or hide' }),
      name: PD.Text('', { description: 'Its name (kinemage.get lists them); a subgroup as it or Group:Subgroup' }),
      visible: PD.Boolean(true, { description: 'Show it (false: hide it)' }),
      kinemage: KinemageParam(),
    },
    when: ifLoaded,
    async run({ kind, name, visible, kinemage }, ctx) {
      const k = resolve(ctx, kinemage);
      if (ApiResult.isFailure(k)) return k;
      const d = k.kinData;
      const dict = kind === 'group' ? d.groupDict : kind === 'subgroup' ? d.subgroupDict : d.masterDict;
      const key = findKey(Object.keys(dict ?? {}), name, kind);
      if (typeof key !== 'string') return key;
      const target = { type: kind, key };
      const label = `${kind} '${key}'`;
      if (isKinemageVisible(ctx.plugin, k.ref, target) === visible) {
        return ApiResult.already(`The ${label} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      await setKinemageVisibility(ctx.plugin, k.ref, d, target, visible);
      return ApiResult.ok(
        `${visible ? 'Showed' : 'Hid'} the ${label}.`,
        { nodes: [] },
        { value: { kind, name: key, visible } },
      );
    },
  });

  readonly view = this.operation({
    id: 'kinemage.view',
    title: 'Apply a kinemage view',
    params: {
      view: PD.Text('', { description: 'A view by name or number; empty: the first' }),
      kinemage: KinemageParam(),
    },
    // A step: a view sets the background and the scene's radius too, which Undo should bring back (G7).
    modules: ['canvas'],
    when: ifLoaded,
    async run({ view, kinemage }, ctx) {
      const k = resolve(ctx, kinemage);
      if (ApiResult.isFailure(k)) return k;
      const views = Object.entries(k.kinData.viewDict ?? {}).map(([key, v]) => ({
        key,
        name: v.name || `View ${key}`,
      }));
      if (!views.length) return ApiResult.fail('not-found', `'${titleOf(k.kinData)}' has no views.`);
      const text = view.trim().toLowerCase();
      const found = text
        ? (views.find((v) => v.key === text) ?? views.find((v) => v.name.toLowerCase() === text))
        : views[0];
      if (!found) {
        return ApiResult.fail(
          'not-found',
          `There is no view '${view}'.`,
          `Views: ${views.map((v) => v.name).join(', ')}.`,
        );
      }
      const snap = (k.kinData.viewSnapshots as Record<string, Partial<Camera.Snapshot>> | undefined)?.[found.key];
      if (!snap) return ApiResult.fail('failed', `The view '${found.name}' has no camera.`);
      await applyViewSnapshot(ctx.plugin, snap);
      return ApiResult.ok(`Applied the view '${found.name}'.`, { view: true }, { value: found });
    },
  });

  readonly animate = this.operation({
    id: 'kinemage.animate',
    title: 'Step a kinemage animation',
    params: {
      sequence: PD.Select<'animate' | '2animate'>(
        'animate',
        [
          ['animate', 'Animate'],
          ['2animate', 'Animate 2'],
        ],
        { description: "Which of the kinemage's frame sequences" },
      ),
      kinemage: KinemageParam(),
    },
    when: ifLoaded,
    async run({ sequence, kinemage }, ctx) {
      const k = resolve(ctx, kinemage);
      if (ApiResult.isFailure(k)) return k;
      const frames = sequence === 'animate' ? k.kinData.groupsAnimate : k.kinData.groupsAnimate2;
      if (!frames?.length) {
        return ApiResult.fail(
          'invalid',
          `'${titleOf(k.kinData)}' has no ${sequence} frames.`,
          'kinemage.get lists them.',
        );
      }
      const shown = await animateKinemage(ctx.plugin, k.ref, k.kinData, sequence);
      return ApiResult.ok(
        `Showing '${shown}' of ${frames.length} ${sequence} frames.`,
        { nodes: [] },
        { value: { shown } },
      );
    },
  });
}

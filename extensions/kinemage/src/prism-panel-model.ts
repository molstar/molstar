/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { KinemageInfo } from './operations.js';

// What the Kinemage panel (prism-panel.tsx) draws from `kinemage.list()`: its layout (the names, which change only when a
// kinemage is loaded or removed) apart from what is shown (which every visibility step changes), so a step rewrites the
// boxes it flips and rebuilds nothing.

export type KinemageItemKind = 'group' | 'subgroup' | 'master';

/** A row of the panel (its eye): what `kinemage.visibility` gets for it. */
export interface KinemageItem {
  readonly kind: KinemageItemKind;
  /** As the operation finds it: a subgroup as `Group:Subgroup`, never ambiguous. */
  readonly key: string;
  readonly label: string;
}

export interface KinemagePanelGroup extends KinemageItem {
  readonly animate?: 1 | 2;
  readonly subgroups: readonly KinemageItem[];
}

export interface KinemagePanelLayout {
  readonly ref: string;
  readonly index: number;
  readonly title: string;
  readonly views: readonly { readonly key: string; readonly name: string }[];
  readonly masters: readonly KinemageItem[];
  /** Without the `nobutton` ones: the v5 panel (and KiNG) leave them out; the operation still reaches them. */
  readonly groups: readonly KinemagePanelGroup[];
  readonly animate: number;
  readonly animate2: number;
}

/** The panel's layout of each kinemage: names only, compared by value (`sameLayout`). */
export function kinemagePanelLayout(list: readonly KinemageInfo[]): readonly KinemagePanelLayout[] {
  return list.map((k) => ({
    ref: k.ref,
    index: k.index,
    title: k.title,
    views: k.views.map((v) => ({ key: v.key, name: v.name })),
    masters: k.masters.map((m) => ({ kind: 'master' as const, key: m.name, label: m.name })),
    groups: k.groups
      .filter((g) => !g.nobutton)
      .map((g) => ({
        kind: 'group' as const,
        key: g.name,
        label: g.name,
        ...(g.animate ? { animate: g.animate } : {}),
        subgroups: g.subgroups
          .filter((s) => !s.nobutton)
          .map((s) => ({ kind: 'subgroup' as const, key: `${g.name}:${s.name}`, label: s.name })),
      })),
    animate: k.animate.length,
    animate2: k.animate2.length,
  }));
}

export function sameLayout(a: readonly KinemagePanelLayout[], b: readonly KinemagePanelLayout[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (x.ref !== y.ref || x.index !== y.index || x.title !== y.title) return false;
    if (x.animate !== y.animate || x.animate2 !== y.animate2) return false;
    if (!sameItems(x.views, y.views, (v, w) => v.key === w.key && v.name === w.name)) return false;
    if (!sameItems(x.masters, y.masters, sameItem)) return false;
    const same = (g: KinemagePanelGroup, h: KinemagePanelGroup) =>
      sameItem(g, h) && g.animate === h.animate && sameItems(g.subgroups, h.subgroups, sameItem);
    if (!sameItems(x.groups, y.groups, same)) return false;
  }
  return true;
}

const sameItem = (a: KinemageItem, b: KinemageItem) => a.kind === b.kind && a.key === b.key && a.label === b.label;

function sameItems<T>(a: readonly T[], b: readonly T[], same: (x: T, y: T) => boolean) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (!same(a[i], b[i])) return false;
  return true;
}

/** A row's key in `kinemageShown`: its kinemage, kind and key. */
export const shownKey = (ref: string, item: Pick<KinemageItem, 'kind' | 'key'>) => `${ref}\n${item.kind}\n${item.key}`;

/** Whether each group, subgroup and master is shown, by `shownKey`: read once a change, then looked up per box. */
export function kinemageShown(list: readonly KinemageInfo[]): ReadonlyMap<string, boolean> {
  const shown = new Map<string, boolean>();
  for (const k of list) {
    for (const m of k.masters) shown.set(shownKey(k.ref, { kind: 'master', key: m.name }), m.visible);
    for (const g of k.groups) {
      shown.set(shownKey(k.ref, { kind: 'group', key: g.name }), g.visible);
      for (const s of g.subgroups) {
        shown.set(shownKey(k.ref, { kind: 'subgroup', key: `${g.name}:${s.name}` }), s.visible);
      }
    }
  }
  return shown;
}

export function sameShown(a: ReadonlyMap<string, boolean>, b: ReadonlyMap<string, boolean>): boolean {
  if (a.size !== b.size) return false;
  for (const [k, v] of a) if (b.get(k) !== v) return false;
  return true;
}

const count = (n: number, one: string, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

/** A kinemage's row: what it holds ("3 groups · 5 masters · 4 views"). */
export function kinemageDescription(k: KinemagePanelLayout): string {
  const parts = [count(k.groups.length, 'group')];
  if (k.masters.length) parts.push(count(k.masters.length, 'master'));
  if (k.views.length) parts.push(count(k.views.length, 'view'));
  return parts.join(' · ');
}

/** The Masters row: how many are on ("4 of 5 on"). */
export function mastersDescription(k: KinemageInfo | undefined): string {
  if (!k) return '';
  return `${k.masters.filter((m) => m.visible).length} of ${k.masters.length} on`;
}

/** A group's row: its subgroups ("3 subgroups"). */
export function groupDescription(group: KinemagePanelGroup): string {
  return group.subgroups.length ? count(group.subgroups.length, 'subgroup') : '';
}

/** An animate sequence's row: the frame shown ("CA trace · 2 of 2"), the first of its groups that is; or none. */
export function framesDescription(k: KinemageInfo | undefined, sequence: 'animate' | '2animate'): string {
  if (!k) return '';
  const frames = sequence === 'animate' ? k.animate : k.animate2;
  const shown = new Set(k.groups.filter((g) => g.visible).map((g) => g.name));
  const i = frames.findIndex((f) => shown.has(f));
  return i < 0 ? `none of ${frames.length} shown` : `${frames[i]} · ${i + 1} of ${frames.length}`;
}

/** What Next frame steps: each kinemage's animate and 2animate sequences (one: the button steps it). */
export function animateMenu(list: readonly KinemagePanelLayout[]) {
  const many = list.length > 1;
  const items: { id: string; label: string; sequence: 'animate' | '2animate'; kinemage: string }[] = [];
  for (const k of list) {
    for (const [sequence, n, name] of [
      ['animate', k.animate, 'Frames'],
      ['2animate', k.animate2, 'Frames (2animate)'],
    ] as const) {
      if (n)
        items.push({
          id: `${k.ref}:${sequence}`,
          label: many ? `${k.title} › ${name}` : name,
          sequence,
          kinemage: k.ref,
        });
    }
  }
  return items;
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewRef } from './layout.js';
import type { RailItem, ShellViewDef } from './shell-model.js';

/**
 * The rail's steps (shell §3.2): its buttons step down as its height shrinks. The sizes are CSS (`@container
 * (max-height: …)` in `shell.module.css`, kept equal to this table by `_test/rail.test.ts`); what fits needs the same
 * numbers. `min` is the least height a step is used at.
 */
export interface RailStep {
  readonly name: 'normal' | 'compact' | 'small';
  readonly min: number;
  /** A button's height (its width is the rail's). */
  readonly button: number;
  /** A button with a caption under its icon. */
  readonly captioned: number;
}

export const RailSteps: readonly RailStep[] = [
  { name: 'normal', min: 620, button: 32, captioned: 44 },
  { name: 'compact', min: 520, button: 28, captioned: 40 },
  { name: 'small', min: 0, button: 24, captioned: 36 },
];

/** The rail's padding at its top and bottom, and between sections (a hairline and its room). */
export const RailPad = 4;
export const RailSectionGap = 8;

export function railStep(height: number): RailStep {
  return RailSteps.find((s) => height >= s.min) ?? RailSteps[RailSteps.length - 1];
}

/** One button of the rail. */
export interface RailEntry {
  readonly id: string;
  /**
   * The section it sits in: `head`, `strip`, a rail item's group, `foot`, `tools` (items without a group), `corner`
   * (items marked `corner`: after the rail's own Hide button, the shell's bottom-right corner).
   */
  readonly section: string;
  /** Higher goes under More sooner. */
  readonly priority: number;
  readonly captioned: boolean;
  readonly view?: { readonly ref: ViewRef; readonly def: ShellViewDef };
  readonly item?: RailItem;
}

/**
 * The rail's buttons in order: the head views, the strip's, the items by group, the foot's views, the other items, the
 * corner's items. A corner item without a priority stays when the rail runs short (it is the corner's one button).
 */
export function railEntries(
  views: readonly { readonly ref: ViewRef; readonly def: ShellViewDef }[],
  items: readonly RailItem[],
): RailEntry[] {
  const all: RailEntry[] = [];
  const view = (v: { ref: ViewRef; def: ShellViewDef }, section: string): RailEntry => ({
    id: `view:${v.ref}`,
    section,
    priority: v.def.priority ?? Number.NaN,
    captioned: false,
    view: v,
  });
  const item = (i: RailItem, section: string): RailEntry => ({
    id: `item:${i.id}`,
    section,
    priority: i.priority ?? Number.NaN,
    captioned: !!i.caption,
    item: i,
  });
  for (const v of views) if (v.def.place === 'head') all.push(view(v, 'head'));
  for (const v of views) if (!v.def.place || v.def.place === 'strip') all.push(view(v, 'strip'));
  const groups: string[] = [];
  for (const i of items) if (i.group !== undefined && !i.corner && !groups.includes(i.group)) groups.push(i.group);
  const byOrder = (a: RailItem, b: RailItem) => (a.order ?? 0) - (b.order ?? 0);
  for (const g of groups)
    for (const i of items.filter((x) => x.group === g && !x.corner).sort(byOrder)) all.push(item(i, `group:${g}`));
  for (const v of views) if (v.def.place === 'foot') all.push(view(v, 'foot'));
  for (const i of items.filter((x) => x.group === undefined && !x.corner).sort(byOrder)) all.push(item(i, 'tools'));
  for (const i of items.filter((x) => x.corner).sort(byOrder))
    all.push({ ...item(i, 'corner'), priority: i.priority ?? -1 });
  // Without a priority: the position in the rail (the foot goes first).
  return all.map((e, i) => (Number.isNaN(e.priority) ? { ...e, priority: i } : e));
}

export interface RailFit {
  readonly step: RailStep;
  /** The entries left out (under More). */
  readonly hidden: ReadonlySet<string>;
}

const NoneHidden: ReadonlySet<string> = new Set();

/**
 * What fits in a rail `height` px tall at its step: everything, or, when even the small step is too tall, the
 * entries of the highest priority numbers (the later ones when equal) leave until the rest and More fit. A height of 0
 * (not measured yet) leaves nothing out. The rail's own button (More) counts as one.
 */
export function fitRail(entries: readonly RailEntry[], height: number, extra = 1): RailFit {
  const step = railStep(height);
  if (height <= 0) return { step, hidden: NoneHidden };
  const total = (kept: readonly RailEntry[], more: boolean) => {
    const sections = new Set(kept.map((e) => e.section));
    let sum = 2 * RailPad + extra * step.button;
    for (const e of kept) sum += e.captioned ? step.captioned : step.button;
    if (more) sum += step.button;
    return sum + Math.max(0, sections.size - 1) * RailSectionGap;
  };
  if (total(entries, false) <= height) return { step, hidden: NoneHidden };
  const order = entries
    .map((e, i) => ({ e, i }))
    .sort((a, b) => b.e.priority - a.e.priority || b.i - a.i)
    .map((x) => x.e);
  const hidden = new Set<string>();
  let kept = entries;
  for (const drop of order) {
    hidden.add(drop.id);
    kept = kept.filter((e) => e !== drop);
    if (total(kept, true) <= height) break;
  }
  return { step, hidden };
}

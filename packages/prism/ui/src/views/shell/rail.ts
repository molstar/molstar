/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable, of, switchMap } from 'rxjs';
import { observableState, type ObservableState } from '../../reactive/observable-state.js';
import type { CanvasToolbarItem, DockModel, ShellModel, ShellTab } from '../../shell/shell.js';

/**
 * The rail's sizing steps (shell.md, Rail): the buttons step down as the rail's height shrinks. The sizes themselves are
 * CSS (`@container (max-height: …)` in `dock.module.css`, which `_test/rail.test.ts` keeps equal to this table); the
 * overflow below needs the same numbers to know what fits, so they are here too. `min` is the smallest height the step
 * is used at (the step before it starts at `min - 1`).
 */
export interface RailStep {
  name: 'normal' | 'compact' | 'small';
  min: number;
  /** A button's height (its width is the rail's). */
  button: number;
  /** A button with a caption under its icon (the command menu). */
  captioned: number;
}

export const RailSteps: readonly RailStep[] = [
  { name: 'normal', min: 620, button: 40, captioned: 46 },
  { name: 'compact', min: 530, button: 34, captioned: 40 },
  { name: 'small', min: 0, button: 30, captioned: 36 },
];

/** The step used at a rail height. */
export function railStep(height: number): RailStep {
  return RailSteps.find((s) => height >= s.min) ?? RailSteps[RailSteps.length - 1];
}

/** One button of the rail as the overflow sees it. */
export interface RailEntry {
  id: string;
  /** The section it sits in (a hairline between sections). */
  section: string;
  /** Lower stays longer. */
  priority: number;
  captioned: boolean;
  /** What a click does. */
  kind: 'tab' | 'tool';
  tab?: ShellTab;
  item?: CanvasToolbarItem;
}

export interface RailFit {
  step: RailStep;
  /** Ids of the entries that do not fit and go under the More button. */
  hidden: ReadonlySet<string>;
}

const NoneHidden: ReadonlySet<string> = new Set();

/**
 * What fits in a rail `height` px tall, at the step that height uses: all of it, or, when even the smallest step is too
 * tall, the entries of the highest `priority` numbers (the later in the rail when they are equal) are left out until the
 * rest and the More button fit. A height of 0 (not measured yet) hides nothing.
 */
export function fitRail(entries: readonly RailEntry[], height: number): RailFit {
  const step = railStep(height);
  if (height <= 0) return { step, hidden: NoneHidden };
  const total = (kept: readonly RailEntry[], more: boolean) => {
    const sections = new Set(kept.map((e) => e.section));
    let sum = 0;
    for (const e of kept) sum += e.captioned ? step.captioned : step.button;
    if (more) {
      sum += step.button;
      sections.add('more');
    }
    // A hairline between sections.
    return sum + Math.max(0, sections.size - 1);
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

/** A toolbar's items grouped in sections, the ungrouped ones (the command menu) first. */
export interface ToolSection {
  group: string | undefined;
  items: readonly CanvasToolbarItem[];
}

export function sectionsOf(items: readonly CanvasToolbarItem[]): ToolSection[] {
  const sections: ToolSection[] = [];
  for (const item of items) {
    const group = item.group;
    const section = sections.find((s) => s.group === group);
    if (section) (section.items as CanvasToolbarItem[]).push(item);
    else sections.push({ group, items: [item] });
  }
  // Items without a section (the command menu) lead.
  return sections.sort((a, b) => (a.group === undefined ? -1 : b.group === undefined ? 1 : 0));
}

/** The rail's buttons in order: head, tabs, the toolbar's sections, the foot's tabs, the ungrouped tools. */
export function railEntries(tabs: readonly ShellTab[], tools: readonly CanvasToolbarItem[]): RailEntry[] {
  const tabEntry = (t: ShellTab, section: string) =>
    ({ id: `tab:${t.id}`, section, priority: t.priority, captioned: false, kind: 'tab', tab: t }) as RailEntry;
  const toolEntry = (i: CanvasToolbarItem, section: string) =>
    ({ id: `tool:${i.id}`, section, priority: i.priority, captioned: !!i.caption, kind: 'tool', item: i }) as RailEntry;
  const all: RailEntry[] = [];
  for (const t of tabs.filter((t) => t.place === 'head')) all.push(tabEntry(t, 'head'));
  for (const t of tabs.filter((t) => !t.place || t.place === 'strip')) all.push(tabEntry(t, 'strip'));
  const sections = sectionsOf(tools);
  for (const s of sections)
    if (s.group !== undefined) for (const i of s.items) all.push(toolEntry(i, `group:${s.group}`));
  for (const t of tabs.filter((t) => t.place === 'foot')) all.push(tabEntry(t, 'foot'));
  for (const s of sections) if (s.group === undefined) for (const i of s.items) all.push(toolEntry(i, 'lead'));
  // Without a priority: the position in the rail.
  return all.map((e, i) => (e.priority === undefined ? { ...e, priority: i } : e));
}

/**
 * The rail's layout: which buttons are left out for the rail's height. Fed the height by the rail's ResizeObserver
 * (`height`), it is recomputed when the height or the set of buttons changes, nothing else; the CSS sizes the steps.
 */
export class RailLayout {
  /** The rail's height in px (0 until measured). */
  readonly height: ObservableState<number> = observableState(0);
  readonly entries: Observable<readonly RailEntry[]>;
  readonly hidden: Observable<ReadonlySet<string>>;
  /**
   * The open menu beside the rail: its bottom edge's distance from the dock's bottom, and the button it belongs to; More's
   * (the left-out buttons) or a rail button's own (`item`: the theme's choices).
   */
  readonly menu: ObservableState<{ bottom: number; anchor: Element; item?: CanvasToolbarItem } | undefined> =
    observableState<{ bottom: number; anchor: Element; item?: CanvasToolbarItem } | undefined>(undefined);
  /** The left-out buttons, in rail order: the More menu's rows. */
  readonly overflow: Observable<readonly RailEntry[]>;

  constructor(shell: ShellModel, dock: DockModel) {
    const tools = shell.toolbar.pipe(
      map((items) => items.filter((i) => i.place !== 'canvas')),
      switchMap((items) =>
        items.length === 0
          ? of([] as CanvasToolbarItem[])
          : combineLatest(items.map((i) => i.visible ?? of(true))).pipe(
              map((flags) => items.filter((_, k) => flags[k])),
            ),
      ),
    );
    this.entries = combineLatest([dock.tabs, tools]).pipe(map(([tabs, items]) => railEntries(tabs, items)));
    const fit = combineLatest([this.entries, this.height.asReadonly]).pipe(
      map(([entries, h]) => ({ entries, fit: fitRail(entries, h) })),
    );
    this.hidden = fit.pipe(
      map((f) => f.fit.hidden),
      distinctUntilChanged((a, b) => a.size === b.size && [...a].every((id) => b.has(id))),
    );
    this.overflow = fit.pipe(
      map((f) => f.entries.filter((e) => f.fit.hidden.has(e.id))),
      distinctUntilChanged((a, b) => a.length === b.length && a.every((e, i) => e === b[i])),
    );
  }

  /** Observes the rail element: its height is the only thing measured, and only when it changes. */
  /**
   * Opens a menu beside the rail at a button, its bottom edge level with the button's (one read, on the click), or
   * closes it when that button's is open: More's (no `item`) or the button's own.
   */
  toggleMenu(anchor: Element, item?: CanvasToolbarItem): void {
    const open = this.menu.value;
    if (open && open.anchor === anchor) return this.menu.set(undefined);
    const root = this.#rail?.parentElement;
    const bottom = root ? root.getBoundingClientRect().bottom - anchor.getBoundingClientRect().bottom : 0;
    this.menu.set({ bottom, anchor, ...(item ? { item } : {}) });
  }

  #rail: Element | undefined;

  observe(rail: Element): () => void {
    this.#rail = rail;
    const ro = new ResizeObserver((entries) => {
      const h = Math.round(entries[entries.length - 1].contentRect.height);
      if (h !== this.height.value) this.height.set(h);
    });
    ro.observe(rail);
    return () => ro.disconnect();
  }
}

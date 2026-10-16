/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

// The UI's layout as state (spec README §3.3, library.md §1): the docks' tabs, sizes and float, the shell's mode and its
// rail, which toolbar buttons show. DOM-free and serializable: what a layout preference keeps and restores. What is
// open is `navigation`'s; `ShellModel` registers the parts and reads both.

export type DockId = 'left' | 'right' | 'top' | 'bottom';
export const DockIds: readonly DockId[] = ['left', 'right', 'top', 'bottom'];

export type ShellMode = 'desktop' | 'compact' | 'embedded';
const ShellModes: readonly ShellMode[] = ['desktop', 'compact', 'embedded'];

/** One dock's layout: what it holds and how it sits. */
export interface DockLayout {
  /** Tab ids in this dock, in order. */
  readonly tabs: readonly string[];
  /** Width (left/right) or height (top/bottom) in pixels. */
  readonly size: number;
  /**
   * Floats over the canvas instead of taking space beside it (desktop mode only; compact mode keeps the bottom sheet).
   * The right and top docks float by default: the canvas runs under them and they read as part of the visualization.
   */
  readonly floating: boolean;
  /**
   * A floating side dock pinned open: the canvas gives it room instead of running under it (some workflows keep a
   * panel open the whole time and want the whole molecule beside it). Only while floating; off by default.
   */
  readonly pinned: boolean;
}

export interface ViewLayout {
  readonly docks: Readonly<Record<DockId, DockLayout>>;
  readonly mode: 'auto' | ShellMode;
  /** The rail column hidden behind a button in the canvas's corner (`ShellOptions.compactRail`). */
  readonly compactRail: boolean;
  /** Which toolbar buttons show, by item id; ids not listed show. */
  readonly buttons: Readonly<Record<string, boolean>>;
}

export const DefaultDockSizes: Readonly<Record<DockId, number>> = { left: 300, right: 340, top: 88, bottom: 220 };
// Side docks are resized by their splitter; a top or bottom dock is sized by its view (the sequence follows its content).
export const MinDockSizes: Readonly<Record<DockId, number>> = { left: 180, right: 180, top: 40, bottom: 40 };

export function defaultViewLayout(): ViewLayout {
  const docks = {} as Record<DockId, DockLayout>;
  for (const id of DockIds) {
    docks[id] = { tabs: [], size: DefaultDockSizes[id], floating: id === 'right' || id === 'top', pinned: false };
  }
  return { docks, mode: 'auto', compactRail: false, buttons: {} };
}

/** The serialized form, versioned: a stored layout of another version is dropped (journal.md §7). */
export interface StoredViewLayout {
  readonly version: 1;
  readonly layout: ViewLayout;
}

export class ViewStateManager {
  readonly #layout = observableState<ViewLayout>(defaultViewLayout());

  constructor(initial?: Partial<Omit<ViewLayout, 'docks'>> & { docks?: Partial<Record<DockId, Partial<DockLayout>>> }) {
    if (initial) this.#layout.set(mergeLayout(this.#layout.value, initial));
  }

  /** The docks, the mode, the rail and the buttons. */
  get layout(): ReadonlyState<ViewLayout> {
    return this.#layout.asReadonly;
  }

  setDock(id: DockId, partial: Partial<DockLayout>): void {
    const current = this.#layout.value;
    const dock = current.docks[id];
    if (!hasChange(dock, partial)) return;
    const next: DockLayout = { ...dock, ...partial };
    if (next.size !== dock.size) (next as { size: number }).size = Math.max(MinDockSizes[id], Math.round(next.size));
    this.#layout.set({ ...current, docks: { ...current.docks, [id]: next } });
  }

  setMode(mode: ViewLayout['mode']): void {
    if (this.#layout.value.mode !== mode) this.#layout.set({ ...this.#layout.value, mode });
  }

  setCompactRail(compactRail: boolean): void {
    if (this.#layout.value.compactRail !== compactRail) this.#layout.set({ ...this.#layout.value, compactRail });
  }

  setButtons(buttons: Readonly<Record<string, boolean>>): void {
    const current = this.#layout.value;
    if (!hasChange(current.buttons, buttons)) return;
    this.#layout.set({ ...current, buttons: { ...current.buttons, ...buttons } });
  }

  /** The layout to keep (a preference), with its version. */
  save(): StoredViewLayout {
    return { version: 1, layout: this.#layout.value };
  }

  /**
   * Restores a kept layout; a value of another version or shape is dropped (false). Tab ids that are not registered
   * stay in the layout, so a tab that registers later lands where it was left.
   */
  restore(stored: unknown): boolean {
    if (!isStoredLayout(stored)) return false;
    const docks: Partial<Record<DockId, Partial<DockLayout>>> = {};
    for (const id of DockIds) {
      const d = stored.layout.docks?.[id];
      // Sizes keep setDock's floor, as a stored value did not come through it.
      if (d)
        docks[id] = { ...d, ...(d.size === undefined ? {} : { size: Math.max(MinDockSizes[id], Math.round(d.size)) }) };
    }
    this.#layout.set(mergeLayout(this.#layout.value, { ...stored.layout, docks }));
    return true;
  }
}

function hasChange<T extends object>(current: T, partial: Partial<T>): boolean {
  for (const key of Object.keys(partial) as (keyof T)[]) {
    if (!Object.is(current[key], partial[key])) return true;
  }
  return false;
}

function mergeLayout(
  base: ViewLayout,
  over: Partial<Omit<ViewLayout, 'docks'>> & { docks?: Partial<Record<DockId, Partial<DockLayout>>> },
): ViewLayout {
  const docks = { ...base.docks };
  for (const id of DockIds) {
    const d = over.docks?.[id];
    if (d) docks[id] = { ...docks[id], ...d };
  }
  return {
    docks,
    mode: over.mode ?? base.mode,
    compactRail: over.compactRail ?? base.compactRail,
    buttons: { ...base.buttons, ...over.buttons },
  };
}

function isStoredLayout(value: unknown): value is StoredViewLayout {
  if (!value || typeof value !== 'object') return false;
  const v = value as Partial<StoredViewLayout>;
  if (v.version !== 1 || !v.layout || typeof v.layout !== 'object') return false;
  const docks = v.layout.docks as Record<string, unknown> | undefined;
  if (!docks || typeof docks !== 'object') return false;
  const { mode, compactRail, buttons } = v.layout;
  if (mode !== undefined && !(mode === 'auto' || (ShellModes as readonly unknown[]).includes(mode))) return false;
  if (compactRail !== undefined && typeof compactRail !== 'boolean') return false;
  if (buttons !== undefined && (!buttons || typeof buttons !== 'object')) return false;
  return DockIds.every((id) => {
    const d = docks[id] as Partial<DockLayout> | undefined;
    if (!d) return true;
    const flag = (f: unknown) => f === undefined || typeof f === 'boolean';
    return (
      Array.isArray(d.tabs) &&
      (d.size === undefined || (typeof d.size === 'number' && Number.isFinite(d.size))) &&
      flag(d.floating) &&
      flag(d.pinned)
    );
  });
}

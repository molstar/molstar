/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import { type DockId, DockIds } from './view-state.js';

// What is open in the UI (spec README §3.3, §3.6): each dock's active tab and whether it is collapsed, the card over
// the canvas and its section, the page over the app, and a tab shown at a section. DOM-free and serializable; the view
// operations (`view.open`) and a light restore write it, `ShellModel` and the views read it.

export interface DockNavigation {
  readonly active: string | undefined;
  /** Collapsed docks show only their tab strip; a collapsed floating dock shows nothing but its rail. */
  readonly collapsed: boolean;
}

export interface NavigationState {
  readonly docks: Readonly<Record<DockId, DockNavigation>>;
  /** The card over the canvas (an overlay), and the section to bring into view. */
  readonly card: { readonly id: string; readonly section?: string } | undefined;
  /** The page over the app, and what to show in it. */
  readonly page: { readonly id: string; readonly focus?: string } | undefined;
  /** A tab shown at a section (the Selection tab at Measure); each request a new object, so asking again scrolls. */
  readonly focus: { readonly tab: string; readonly section: string } | undefined;
  /** In a compact rail: whether its column is shown. */
  readonly railOpen: boolean;
}

export function defaultNavigation(): NavigationState {
  const docks = {} as Record<DockId, DockNavigation>;
  for (const id of DockIds) docks[id] = { active: undefined, collapsed: false };
  return { docks, card: undefined, page: undefined, focus: undefined, railOpen: false };
}

export interface StoredNavigation {
  readonly version: 1;
  readonly navigation: Pick<NavigationState, 'docks' | 'card' | 'page'>;
}

/** The main dock: the panel the rail opens (`tab` is its open tab). */
const MainDock: DockId = 'right';

export class NavigationManager {
  readonly #state = observableState<NavigationState>(defaultNavigation());
  readonly #tab: ReadonlyState<string | undefined>;
  readonly #card: ReadonlyState<string | undefined>;
  readonly #section: ReadonlyState<string | undefined>;
  readonly #page: ReadonlyState<string | undefined>;

  constructor(initial?: { docks?: Partial<Record<DockId, Partial<DockNavigation>>>; railOpen?: boolean }) {
    if (initial) {
      const current = this.#state.value;
      const docks = { ...current.docks };
      for (const id of DockIds) if (initial.docks?.[id]) docks[id] = { ...docks[id], ...initial.docks[id] };
      this.#state.set({ ...current, docks, railOpen: initial.railOpen ?? current.railOpen });
    }
    this.#tab = derived(this.#state.asReadonly, (s) => openTab(s, MainDock));
    this.#card = derived(this.#state.asReadonly, (s) => s.card?.id);
    this.#section = derived(this.#state.asReadonly, (s) => s.card?.section ?? s.focus?.section);
    this.#page = derived(this.#state.asReadonly, (s) => s.page?.id);
  }

  get state(): ReadonlyState<NavigationState> {
    return this.#state.asReadonly;
  }

  /** The main panel's open tab; undefined while it is collapsed. */
  get tab(): ReadonlyState<string | undefined> {
    return this.#tab;
  }

  /** The open card, by overlay id. */
  get card(): ReadonlyState<string | undefined> {
    return this.#card;
  }

  /** The section asked for: the open card's, or the focused tab's. */
  get section(): ReadonlyState<string | undefined> {
    return this.#section;
  }

  /** The open page, by id. */
  get page(): ReadonlyState<string | undefined> {
    return this.#page;
  }

  /** The open tab of a dock; undefined while it is collapsed. */
  openTab(dock: DockId): string | undefined {
    return openTab(this.#state.value, dock);
  }

  setDock(id: DockId, partial: Partial<DockNavigation>): void {
    const current = this.#state.value;
    const dock = current.docks[id];
    const same =
      (!('active' in partial) || partial.active === dock.active) &&
      (partial.collapsed === undefined || partial.collapsed === dock.collapsed);
    if (same) return;
    this.#state.set({ ...current, docks: { ...current.docks, [id]: { ...dock, ...partial } } });
  }

  setCard(card: NavigationState['card']): void {
    const current = this.#state.value;
    if (current.card?.id === card?.id && current.card?.section === card?.section && !card?.section) return;
    this.#state.set({ ...current, card });
  }

  setPage(page: NavigationState['page']): void {
    const current = this.#state.value;
    if (current.page?.id === page?.id && current.page?.focus === page?.focus) return;
    this.#state.set({ ...current, page });
  }

  setFocus(focus: NavigationState['focus']): void {
    const current = this.#state.value;
    if (current.focus === focus) return;
    this.#state.set({ ...current, focus });
  }

  setRailOpen(railOpen: boolean): void {
    if (this.#state.value.railOpen !== railOpen) this.#state.set({ ...this.#state.value, railOpen });
  }

  save(): StoredNavigation {
    const { docks, card, page } = this.#state.value;
    return { version: 1, navigation: { docks, card, page } };
  }

  /** Restores what was open; a value of another version or shape is dropped (false). */
  restore(stored: unknown): boolean {
    if (!isStoredNavigation(stored)) return false;
    const current = this.#state.value;
    const docks = { ...current.docks };
    for (const id of DockIds) {
      const d = stored.navigation.docks[id];
      if (d) docks[id] = { active: d.active, collapsed: !!d.collapsed };
    }
    this.#state.set({ ...current, docks, card: stored.navigation.card, page: stored.navigation.page });
    return true;
  }
}

function openTab(s: NavigationState, dock: DockId): string | undefined {
  const d = s.docks[dock];
  return d.collapsed ? undefined : d.active;
}

function derived<T>(source: ReadonlyState<NavigationState>, pick: (s: NavigationState) => T): ReadonlyState<T> {
  const out = source.pipe(map(pick), distinctUntilChanged()) as ReadonlyState<T>;
  Object.defineProperty(out, 'value', { get: () => pick(source.value) });
  return out;
}

function isStoredNavigation(value: unknown): value is StoredNavigation {
  if (!value || typeof value !== 'object') return false;
  const v = value as Partial<StoredNavigation>;
  return v.version === 1 && !!v.navigation && typeof v.navigation.docks === 'object' && !!v.navigation.docks;
}

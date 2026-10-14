/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React UI's context and hooks (seqspec:interaction#the-default-ui): parts read the smallest state they render
// through `useSyncExternalStore`, so a change re-renders only the parts that read it; pointer-rate values are read in
// subscriptions that write the DOM, never in render.

import {
  type ComponentType,
  createContext,
  type ReactElement,
  type ReactNode,
  type Ref,
  type RefObject,
  useCallback,
  useContext,
  useLayoutEffect,
  useRef,
  useSyncExternalStore,
} from 'react';
import { type Observable, skip } from 'rxjs';
import type { MenuItem } from '../chrome.js';
import { landFocus, type SurfaceFactory } from '../dom/core.js';
import type { LegendMode, Live, SearchFn, SeqState, SeqViewerContext } from '../index.js';
import { sameShallow } from '../same.js';

// The surface contract is the DOM UI's: one declaration for both.
export type { MountedSurface, SurfaceFactory } from '../dom/core.js';

export interface PartProps {
  readonly context: SeqViewerContext;
}
export interface SectionProps extends PartProps {
  readonly section: string;
}
export interface TrackProps extends SectionProps {
  readonly track: string;
  /** The header line's layout row, read when the menu is built (a split row's index; −1: the track's own line). */
  readonly row?: () => number;
}
export interface RowProps extends SectionProps {
  readonly displayRow: number;
  readonly dataRow: number;
}
export interface HintProps {
  readonly label: string;
  readonly children: ReactElement;
}
export interface ScrollAreaProps {
  readonly children: ReactNode;
  readonly viewportRef?: Ref<HTMLElement>;
  readonly className?: string;
}
/**
 * Controlled, as Base UI's: the part that owns the menu holds `open`; seqspec:interaction#slots-and-chrome has its
 * contract.
 */
export interface MenuProps {
  /** The button that opens it; the menu adds its own handlers and ARIA state. */
  readonly trigger: ReactElement<{ onClick?: () => void; 'aria-expanded'?: boolean }>;
  readonly items: readonly MenuItem[];
  readonly open: boolean;
  /** Asked for on the trigger, Escape, Tab, an item chosen and a press outside; the owner answers with `open`. */
  onOpenChange(open: boolean): void;
}

/**
 * Every part, and the primitives they are drawn with (a host maps those to its own). Each part is also a component a
 * host places itself, anywhere under the provider (seqspec:interaction#slots-and-chrome).
 */
export interface SlotComponents {
  ViewTabs: ComponentType<PartProps>;
  /** A section's label, its filter chip, its `HiddenChips` and its `SectionSwitches`. */
  SectionHeader: ComponentType<SectionProps>;
  /** A section's switches (numbering, copies, ligand form, rows, style), each an action. */
  SectionSwitches: ComponentType<SectionProps>;
  TrackActions: ComponentType<TrackProps>;
  HiddenChips: ComponentType<SectionProps>;
  RowControl: ComponentType<RowProps>;
  Tooltip: ComponentType<PartProps>;
  /** Item actions on hover (seqspec:interaction#the-hover-card); inert without `UI.ItemActions`. */
  HoverCard: ComponentType<PartProps>;
  Legend: ComponentType<PartProps>;
  SnapshotPlayer: ComponentType<PartProps>;
  Description: ComponentType<PartProps>;
  /** Only when the UI spec gives `search`. */
  Search: ComponentType<PartProps>;
  Export: ComponentType<PartProps>;
  Issues: ComponentType<PartProps>;
  Busy: ComponentType<PartProps>;
  EmptyState: ComponentType<PartProps>;
  ErrorState: ComponentType<PartProps>;
  /** Hints, never a `title` attribute. */
  Hint: ComponentType<HintProps>;
  /** DOM scroll regions. */
  ScrollArea: ComponentType<ScrollAreaProps>;
  Menu: ComponentType<MenuProps>;
}
export type SlotName = keyof SlotComponents;

export interface SeqViewerUISpec {
  /**
   * A component, "none" (the part is left out wherever it is placed), or "default" (its default, also placed under
   * `clean` chrome).
   */
  readonly components?: { readonly [K in SlotName]?: SlotComponents[K] | 'none' | 'default' };
  /** On the root, where a host scopes its token mapping. */
  readonly className?: string;
  /** The parts `SeqViewerLayout` places (`ChromeCompositions`); unset: the context's `UI.Chrome`. Never written to it. */
  readonly chrome?: 'full' | 'clean';
  /**
   * This UI's legend (seqspec:interaction#the-legend): collapsed, expanded or hidden; unset: the context's `UI.Legend`.
   */
  readonly legend?: LegendMode;
  /** The canvas engine's `mountSection` by default. A stable function: a new one remounts every surface. */
  readonly surface?: SurfaceFactory;
  /** An application's search; absent: no `Search` part. A stable function, as `surface`. */
  readonly search?: SearchFn;
}
export interface ResolvedUISpec {
  /** Each slot's part; null: "none" (or `RowControl` without one). */
  readonly components: { readonly [K in SlotName]: SlotComponents[K] | null };
  /** The parts the default layout places (`placedParts`): the chrome's composition and the host's slots. */
  readonly placed: ReadonlySet<string>;
  readonly className: string;
  /** The spec's `legend`, if given. */
  readonly legend?: LegendMode;
  readonly surface: SurfaceFactory;
  readonly search?: SearchFn;
}

export interface ViewerScope {
  readonly context: SeqViewerContext;
  readonly ui: ResolvedUISpec;
  /** Mounted surface elements by section ref (the tooltip places itself from them); one map per provider instance. */
  readonly surfaces: Map<string, HTMLElement>;
}
export const ScopeContext = createContext<ViewerScope | null>(null);

export function useScope(): ViewerScope {
  const s = useContext(ScopeContext);
  if (!s) throw new Error('Seq* parts render under SeqViewerProvider');
  return s;
}
/** The provider's scope, or null: a part placed without one draws with the default primitives. */
export const useScopeMaybe = (): ViewerScope | null => useContext(ScopeContext);
export function useSeqViewer(): { context: SeqViewerContext; ui: ResolvedUISpec } {
  return useScope();
}

type Subscribable = Observable<unknown>;
/**
 * `useSyncExternalStore`'s subscribe over states: a state replays its value on subscribe (as a BehaviorSubject), which
 * is not a change, so the first emission of each is skipped; calling back on it made React check the store once more.
 */
export const subscriber = (states: readonly Subscribable[]) => (cb: () => void) => {
  const subs = states.map((s) => s.pipe(skip(1)).subscribe(() => cb()));
  return () => {
    for (const s of subs) s.unsubscribe();
  };
};

/** The state's value; re-renders when it changes (a state's replay on subscribe is not a change). */
export function useSeqState<T>(s: SeqState<T>): T {
  const subscribe = useCallback((cb: () => void) => subscriber([s])(cb), [s]);
  return useSyncExternalStore(subscribe, () => s.value);
}

/**
 * Re-renders on the live record's version; read `l.current` in render. A live record is one object changed in place
 * (the hover's, the viewport's), so what a render reads is what it is *now*: read the hover's loci
 * (`ctx.hover.toLoci()`), the way `<seq-hover>` carries them, in the subscription or handler that saw the change, never
 * in a later render, where the hover has moved on.
 */
export function useLiveVersion(l: Live<unknown>): number {
  return useSeqState(l.version);
}

export interface DerivedOptions<T> {
  /** When the new value is the same as the last (shallow by default; `sameData` for a list of records, or your own). */
  readonly compare?: (a: T, b: T) => boolean;
  /** What `compute` reads besides the states (props, a context): it is computed again when one changes. */
  readonly deps?: readonly unknown[];
}
const NO_DEPS: readonly unknown[] = [];
const sameList = (a: readonly unknown[], b: readonly unknown[]) =>
  a.length === b.length && a.every((x, i) => x === b[i]);
/**
 * A value derived from states, computed once and again only when one of the states changed or one of `deps` did (not on
 * each render, which React's snapshot reads would make of it); the previous object is kept while the new one is the same
 * by `compare` (shallow by default), so a part re-renders only when what it shows changed. `compute` reads the states
 * and `deps` and nothing else. Edge-rate values only.
 */
export function useDerived<T>(
  states: readonly SeqState<unknown>[],
  compute: () => T,
  { compare = sameShallow, deps = NO_DEPS }: DerivedOptions<T> = {},
): T {
  const fn = useRef(compute);
  fn.current = compute;
  const memo = useRef<{
    states: readonly SeqState<unknown>[];
    deps: readonly unknown[];
    seen: unknown[];
    value: { v: T } | null;
    dirty: boolean;
  }>({ states, deps, seen: [], value: null, dirty: true });
  const m = memo.current;
  if (!sameList(m.states, states) || !sameList(m.deps, deps)) {
    [m.states, m.deps, m.dirty] = [states, deps, true];
  }
  const watched = (cb: () => void) =>
    subscriber(states)(() => {
      memo.current.dirty = true;
      cb();
    });
  // biome-ignore lint/correctness/useExhaustiveDependencies: the states' identities are the subscription.
  const subscribe = useCallback(watched, states);
  const get = () => {
    // A state may have changed between a render and the subscription (which skips the replay): its value shows it.
    if (!m.dirty) for (let i = 0; i < states.length; i++) if (states[i].value !== m.seen[i]) m.dirty = true;
    if (m.dirty || !m.value) {
      m.dirty = false;
      m.seen = states.map((s) => s.value);
      const v = fn.current();
      if (!m.value || !compare(m.value.v, v)) m.value = { v };
    }
    return m.value.v;
  };
  return useSyncExternalStore(subscribe, get);
}

/**
 * Focus inside the element when React removes it lands where seqspec:interaction#accessibility says (`landFocus`): its
 * group's next or previous node, the section's surface, the viewer's root. The cleanup runs before the node leaves the
 * document.
 * `near` (a stable function) is where the landing is looked for, when it is not the element itself.
 */
export function useLanding(ref: RefObject<HTMLElement | null>, near?: (el: HTMLElement) => Node) {
  useLayoutEffect(() => {
    const el = ref.current;
    return () => {
      if (el) landFocus(el, near?.(el) ?? el);
    };
  }, [ref, near]);
}

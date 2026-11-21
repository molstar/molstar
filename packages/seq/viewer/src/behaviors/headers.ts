/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { HeaderPartsOf, Hit as RenderHit } from '@molstar/seqx-core/render';
import type { CommandDefinition } from '@molstar/ui/commands/types';
import { describerFor } from '@molstar/seqx-viewer/describe';
import {
  type HeaderActionEntry,
  type HeaderEntry,
  type HeaderFacts,
  type HeaderGeometry,
  type HeaderMenu,
  type HeaderMenuItem,
  type HeaderMenuPlace,
  type HeaderMenuPort,
  type HeaderPartsPort,
  type HeaderPartsState,
  type HeaderTarget,
  type HeadersParams,
} from '@molstar/seqx-viewer/header-entries';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { HeaderExtras } from './header-extras.js';
import type { ThemingBehavior } from './theming.js';

// The `headers` behavior (_spec/README.md "Behaviors", Headers): the header entries (the registry's `header-action`
// entries and a host's at run time, with their undo) made into the layout's parts once a layout (`partsOf`, the frame's
// `headers`), a click on a header's part, the secondary click under `contextClick` (Shift+secondary click always the
// browser's) and the long press, the menu key once per ask (EC-156), the geometry (`viewer.headers`) and a host's DOM
// parts (`headerParts`), each published on change only. A hover redraws the actions with no layout and nothing
// published. What the hover path does not need (the menu's data, a run, the geometry, the words, the cursor on
// headers) is the late parts' (`header-extras.ts`, a lazy chunk).

/** The menu shown: what, where, and whether the keyboard asked (its focus goes back to the canvas). */
export interface ShownHeaderMenu {
  readonly menu: HeaderMenu;
  readonly place: HeaderMenuPlace;
}

export interface HeadersBehavior {
  readonly commands: readonly CommandDefinition<unknown>[];
  /**
   * @internal The layout's parts (the frame's `headers`): a new function when the entries change (a new layout);
   * undefined until the late parts land (core's built-in states drawn meanwhile).
   */
  readonly partsOf: HeaderPartsOf | undefined;
  /** The entries in order: the registry's, then a host's registered at run time. */
  readonly entries: readonly HeaderEntry[];
  /** Every header's geometry (`viewer.headers`), published when a layout changes it. */
  readonly geometry: ReadonlyState<readonly HeaderGeometry[]>;
  /** The menu shown, undefined when none. */
  readonly menu: ReadonlyState<ShownHeaderMenu | undefined>;
  /** A host's DOM parts (`headerParts`), published on change while a parts port is given; undefined otherwise. */
  readonly parts: ReadonlyState<HeaderPartsState | undefined>;
  /** Registers entries at run time, atomically (a conflict throws before any is added); returns their undo. */
  register(...entries: HeaderEntry[]): () => void;
  /** A header's facts now (from the resolved view and the last layout); undefined when it is not shown. */
  facts(target: HeaderTarget): HeaderFacts | undefined;
  /** The menu's data for a target (DOM-free; once the late parts landed). */
  menuOf(target: HeaderTarget): HeaderMenu | undefined;
  /** The parts' items of every shown header (DOM-free; what the port is given). */
  partsState(): HeaderPartsState;
  /** Runs an entry for a target (a slot's press, a DOM part's): its command's message, undefined when none ran. */
  run(name: string, target: HeaderTarget): string | Promise<string> | undefined;
  /** Runs a menu item's command (a choice). */
  runItem(item: HeaderMenuItem): string | Promise<string> | undefined;
  /** @internal Previews a menu item while it is hovered (a scheme), when its command has a preview; false otherwise. */
  preview(item: HeaderMenuItem): boolean;
  /** Opens the menu of a target at a place (`only`: an entry's own items, a chip's list or a slot's submenu). */
  openMenu(target: HeaderTarget, place: HeaderMenuPlace, only?: string): boolean;
  /** Closes the menu; a preview it made ends in nothing. */
  closeMenu(): void;
  /** @internal A click on a header's part at a canvas point (the input's, taken). */
  click(x: number, y: number, touch: boolean): void;
  /** @internal A secondary click at a canvas point: true when the menu took it (the browser's menu then prevented). */
  context(x: number, y: number, e?: { readonly shiftKey?: boolean }): boolean;
  composed(viewer: SeqViewer): void;
  /** @internal Resolves once its lazy chunks (the late parts, the icons, the gestures, the ports') landed. */
  loaded(): Promise<void>;
  /** @internal A frame was drawn (each frame): the geometry and the parts follow it; a hover makes nothing. */
  framed(): void;
  attach(container: HTMLElement, viewer: SeqViewer): () => void;
  dispose(): void;
}

/** The command that opens the menu at the keyboard cursor (the menu key, Shift+F10; Enter on a header). */
export const HeadersMenuCommand = 'headers.menu';

/** @internal What the late parts read and set of the behavior (`header-extras.ts`). */
export interface HeaderShared {
  readonly behavior: HeadersBehavior;
  readonly context: BehaviorContext;
  readonly params: HeadersParams;
  viewer: SeqViewer;
  disposed: boolean;
  /** A menu item's preview is showing (a choice or a close ends it). */
  previewed: boolean;
  /** A host's DOM parts draw the parts (their port's chunk landed). */
  hosted(): boolean;
  entries(): readonly HeaderEntry[];
  /** A lazy chunk asked for (`loaded` waits for it; ignored once retired). */
  lazy<M>(load: Promise<M>, landed: (m: M) => void): void;
  /** The headers laid out again (the icons landed). */
  relayout(): void;
  geometry(): readonly HeaderGeometry[];
  /** Shows a menu (its state and its port). */
  show(menu: HeaderMenu, place: HeaderMenuPlace): void;
}

function headersBehavior(context: BehaviorContext, params: HeadersParams = {}): HeadersBehavior {
  const dom = typeof params.parts === 'function' ? params.parts : undefined;
  /**
   * The lazy chunks asked for (`loaded`, what `settle` waits for). A chunk landing after the viewer is retired is
   * ignored, its failure too (a torn-down test environment); before, a failure is the host's to see.
   */
  const loading: Promise<void>[] = [];
  const lazy = <M>(load: Promise<M>, landed: (m: M) => void): void => {
    loading.push(
      load.then(
        (m) => void (shared.disposed || landed(m)),
        (e) => {
          if (!shared.disposed) throw e;
        },
      ),
    );
  };
  /** The registry's entries and a host's, in order; a new `partsOf` each change. */
  let entries: readonly HeaderEntry[] = [];
  /** The layout's parts: the late parts' (until they land, core's built-in states). */
  let partsOf: HeaderPartsOf | undefined;
  const geometry = new State<readonly HeaderGeometry[]>([]);
  const shownMenu = new State<ShownHeaderMenu | undefined>(undefined);
  const parts = new State<HeaderPartsState | undefined>(undefined);
  let menuPort: HeaderMenuPort | undefined;
  let menuPortMounted = false;
  let partsPort: HeaderPartsPort | undefined;
  /** A host's DOM parts draw the slots, the "⋯" and the chips: once their port's chunk landed (until then the canvas). */
  let hosted = false;
  /** A click before the late parts landed: run once they land (one at most, the last). */
  let queued: [RenderHit, number, number, boolean] | undefined;
  const stops: (() => void)[] = [];
  /** The late parts (`header-extras.ts`, a lazy chunk asked for once composed). */
  let extras: HeaderExtras | undefined;

  const commandOf = (id: string | undefined) => (id ? shared.viewer?.commands.find((c) => c.id === id) : undefined);

  /** Checks entries against those kept and the viewer's commands: every conflict in one error. */
  const check = (adding: readonly HeaderEntry[], kept: readonly HeaderEntry[]) => {
    const problems: string[] = [];
    const names = new Set(kept.map((e) => e.name));
    for (const e of adding) {
      if (!/^[\w-]+:[\w.-]+$/.test(e.name)) problems.push(`${e.name}: not namespaced (<id>:<name>)`);
      if (names.has(e.name)) problems.push(`${e.name}: registered twice`);
      names.add(e.name);
      if (e.command !== undefined && shared.viewer && !commandOf(e.command))
        problems.push(`${e.name}: no command ${e.command} (seq.viewer.header.command)`);
    }
    if (problems.length) throw new Error(`The header entries are refused: ${problems.join('; ')}`);
  };
  /** The headers laid out again: a new function (part of the layout's key), the late parts' once they landed. */
  const relayout = () => {
    const x = extras;
    partsOf = x ? (f) => x.layoutParts(f) : undefined;
    shared.viewer?.schedule();
  };
  const setEntries = (next: readonly HeaderEntry[]) => {
    entries = next;
    relayout();
  };

  // ---- the menu ----

  const port = (): HeaderMenuPort | undefined => {
    if (!menuPort && params.menu && !menuPortMounted) menuPort = params.menu(undefined);
    return menuPort;
  };
  const closeMenu = () => {
    if (shared.previewed) {
      shared.previewed = false;
      (shared.viewer?.behaviors.get('legend') as ThemingBehavior | undefined)?.cancelPreview();
    }
    if (!shownMenu.value) return;
    shownMenu.set(undefined);
    menuPort?.hide();
  };

  // ---- geometry ----

  // the last layout published and the headers whose actions showed (the parts: a hover moves only these)
  let publishedFor: unknown;
  let showingKey = '';
  const published = () => {
    const viewer = shared.viewer;
    const layout = viewer?.frame?.layout;
    if (!layout) return;
    const relaid = layout !== publishedFor;
    if (relaid && extras) {
      publishedFor = layout;
      geometry.set(extras.published(layout));
    }
    if (!partsPort || !extras) return;
    // the parts follow the hover only when the header whose actions show changed (a key of its address)
    const c = viewer.selection.cursor.value;
    const m = shownMenu.value?.menu.target;
    const key = `${viewer.hover.current.layer ?? ''}|${c?.track ?? ''}|${m ? (m.track ?? m.section) : ''}`;
    if (!relaid && key === showingKey) return focusPending();
    showingKey = key;
    const next = extras.partsState();
    parts.set(next);
    partsPort.update(next);
    focusPending();
  };
  /** The header the keyboard cursor came onto: its parts take the focus once they are drawn for it. */
  let pendingFocus: HeaderTarget | undefined;
  const focusPending = () => {
    const t = pendingFocus;
    pendingFocus = undefined;
    if (t) partsPort?.focus?.(t);
  };

  // ---- a click ----

  const click = (x: number, y: number, touch: boolean): void => {
    const hit = shared.viewer?.hit(x, y, touch);
    const t = hit?.target;
    if (!hit || (t !== 'header-slot' && t !== 'header-state' && t !== 'header-label' && t !== 'section-header')) return;
    // a header's parts are never a selection's: run by the late parts; before they land, once they land
    if (extras) extras.click(hit, x, y, touch);
    else queued = [hit, x, y, touch];
  };

  const commands: CommandDefinition<unknown>[] = [
    {
      id: HeadersMenuCommand,
      title: 'Open the actions of the track at the cursor',
      category: 'Headers',
      scope: 'canvas',
      when: () => (shared.disposed ? 'The viewer is retired' : extras?.cursorTarget() ? true : 'No cursor'),
      run: () => extras?.menuAtCursor() ?? { message: '' },
    },
  ];

  const behavior: HeadersBehavior = {
    commands,
    get partsOf() {
      return partsOf;
    },
    get entries() {
      return entries;
    },
    geometry: geometry.asReadonly,
    menu: shownMenu.asReadonly,
    parts: parts.asReadonly,
    register(...add) {
      check(add, entries);
      setEntries([...entries, ...add]);
      let undone = false;
      return () => {
        if (undone) return;
        undone = true;
        setEntries(entries.filter((e) => !add.includes(e)));
      };
    },
    facts: (target) => extras?.factsOf(target),
    menuOf: (target) => extras?.menuOf(target),
    partsState: () => extras?.partsState() ?? { headers: [] },
    run: (name, target) => extras?.run(name, target),
    runItem: (item) => extras?.runItem(item),
    preview: (item) => !!extras?.preview(item),
    openMenu: (target, place, only) => !!extras?.openMenu(target, place, only),
    closeMenu,
    framed() {
      extras?.keepCursor();
      published();
    },
    click,
    context: (x, y, e) => !!extras?.context(x, y, e),
    async loaded() {
      // a port's chunk asked for meanwhile (the menu's at its first open): waited for too
      for (let n = -1; n !== loading.length; ) {
        n = loading.length;
        await Promise.all([...loading, menuPort?.loaded?.(), partsPort?.loaded?.()]);
      }
    },
    composed(v) {
      shared.viewer = v;
      const listed = (v.registry ?? []).flatMap((e) =>
        e.kind === 'header-action' ? [(e as HeaderActionEntry).entry] : [],
      );
      check(listed, []);
      setEntries(listed);
      lazy(import('./header-extras.js'), (m) => {
        const x = m.headerExtras(shared);
        extras = x;
        relayout();
        describerFor(context).headers = x.words;
        publishedFor = undefined;
        published();
        const q = queued;
        queued = undefined;
        if (q) x.click(...q);
      });
      // a host's DOM parts take the focus when the cursor comes onto a header
      // (redrawn for that header first: a redraw replaces the buttons, and with them the focus)
      stops.push(
        context.selection.cursor.listen((c) => {
          if (c?.keyboard && c.track !== undefined && partsPort?.focus) {
            // the next frame draws the header's actions (the emphasis follows the cursor there), `published` then
            // redraws the parts and the focus goes in
            pendingFocus = c.track === c.section ? { section: c.section } : { section: c.section, track: c.track };
          }
        }),
      );
      // a new layout, an undo or a load: the menu was of the old picture; a preview's own (a scheme hovered in the
      // menu) is the menu's, and keeps it
      let preview = context.spec.preview.value;
      stops.push(
        context.view.resolved.listen(() => {
          const p = context.spec.preview.value;
          if (p !== preview) preview = p;
          else closeMenu();
        }),
      );
    },
    attach(container) {
      // the menu's port over the canvas (a headless one given up first)
      if (params.menu) {
        if (menuPort && menuPort.disposeWithViewer !== false) menuPort.dispose?.();
        menuPort = params.menu(container);
        menuPortMounted = true;
      }
      if (dom) {
        const p = dom(container, behavior);
        partsPort = p;
        showingKey = '\u0000';
        publishedFor = undefined;
        published();
        // the canvas draws the parts until the port's DOM is in place, then nothing in their boxes
        lazy(p.loaded?.() ?? Promise.resolve(), () => {
          if (partsPort !== p) return;
          hosted = true;
          relayout();
        });
      }
      // the secondary click (Shift: the browser's) and the long press (touch's context menu, once per ask)
      // (a lazy chunk, loaded at the mount)
      let stopGestures: (() => void) | undefined;
      let detached = false;
      lazy(import('@molstar/seqx-viewer/header-gestures'), (m) => {
        if (detached) return;
        stopGestures = m.headerGestures(
          container,
          params.longPressMs ?? 500,
          (x, y, e) => !!extras?.context(x, y, e),
          (x, y) => extras?.pressed(x, y),
        );
      });
      return () => {
        detached = true;
        stopGestures?.();
        closeMenu();
        if (menuPort && menuPort.disposeWithViewer !== false) menuPort.dispose?.();
        menuPort = undefined;
        menuPortMounted = false;
        if (partsPort && partsPort.disposeWithViewer !== false) partsPort.dispose?.();
        partsPort = undefined;
        parts.set(undefined);
        if (hosted) {
          hosted = false;
          relayout();
        }
      };
    },
    dispose() {
      shared.disposed = true;
      closeMenu();
      for (const s of stops.splice(0)) s();
      if (menuPort && menuPort.disposeWithViewer !== false) menuPort.dispose?.();
      menuPort = undefined;
      const d = describerFor(context);
      if (extras && d.headers === extras.words) d.headers = undefined;
    },
  };

  const shared: HeaderShared = {
    behavior,
    context,
    params,
    viewer: undefined as unknown as SeqViewer,
    disposed: false,
    previewed: false,
    hosted: () => hosted,
    entries: () => entries,
    lazy,
    relayout,
    geometry: () => geometry.value,
    show(menu, place) {
      shownMenu.set({ menu, place });
      port()?.show(menu, place, behavior);
    },
  };
  return behavior;
}

/** The `headers` behavior (_spec/README.md "Behaviors", Headers): slots, chips, the menu and the parts of the headers. */
export const HeadersProvider: BehaviorProvider<HeadersParams | undefined, HeadersBehavior> = {
  name: 'headers',
  create: (context, params) => headersBehavior(context, params ?? {}),
};

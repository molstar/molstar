/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { resolveAddress } from '@molstar/core/repr/tree';
import {
  actionColumns,
  type HeaderInputFacts,
  type HeaderParts,
  type HeaderPartsOf,
  type HeaderSlotPart,
  type HeaderStatePart,
} from '@molstar/seqx-core/render';
import { hiddenTracksOf } from '@molstar/seqx-core/render/headers';
import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import { formatChord, isMac } from '@molstar/ui/commands/keys';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom/host';
import { parseChord } from '@molstar/seqx-viewer/bindings';
import { describerFor } from '@molstar/seqx-viewer/describe';
import {
  type HeaderActionEntry,
  type HeaderEntry,
  type HeaderFacts,
  type HeaderGeometry,
  type HeaderLayerFacts,
  type HeaderMenu,
  type HeaderMenuItem,
  type HeaderMenuPlace,
  type HeaderMenuPort,
  type HeaderPartsPort,
  type HeaderPartsState,
  type HeaderTarget,
  type HeadersParams,
  sameGeometry,
} from '@molstar/seqx-viewer/header-entries';
import type { HeaderIcon } from '@molstar/seqx-core/render';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { HeaderExtras } from './header-extras.js';
import type { HoverCardBehavior } from './hover-card.js';
import type { ThemingBehavior } from './theming.js';

// The `headers` behavior (_spec/README.md "Behaviors", Headers): the header entries (the registry's `header-action`
// entries and a host's at run time, with their undo) made into the layout's parts once a layout (`partsOf`, the frame's
// `headers`), a slot's or a chip's click run, the menu (`headerMenu`, shown by its port), the secondary click under
// `contextClick` (Shift+secondary click always the browser's) and the long press, the menu key once per ask (EC-156),
// the geometry (`viewer.headers`) and a host's DOM parts (`headerParts`), each published on change only. A hover redraws
// the actions with no layout and nothing published (the geometry is the layout's, compared by identity first).

/** The menu shown: what, where, and whether the keyboard asked (its focus goes back to the canvas). */
export interface ShownHeaderMenu {
  readonly menu: HeaderMenu;
  readonly place: HeaderMenuPlace;
}

export interface HeadersBehavior {
  readonly commands: readonly CommandDefinition<unknown>[];
  /** The layout's parts (the frame's `headers`): a new function when the entries change (a new layout). */
  readonly partsOf: HeaderPartsOf;
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
  /** The menu's data for a target (DOM-free). */
  menuOf(target: HeaderTarget): HeaderMenu | undefined;
  /** The parts' items of every shown header (DOM-free; what the port is given). */
  partsState(): HeaderPartsState;
  /** Runs an entry for a target (a slot's press, a DOM part's): its command's message, undefined when none ran. */
  run(name: string, target: HeaderTarget): string | Promise<string> | undefined;
  /** Runs a menu item's command (a choice). */
  runItem(item: HeaderMenuItem): string | Promise<string> | undefined;
  /** Previews a menu item while it is hovered (a scheme), when its command has a preview; false otherwise. */
  preview(item: HeaderMenuItem): boolean;
  /** Opens the menu of a target at a place (`only`: an entry's own items, a chip's list or a slot's submenu). */
  openMenu(target: HeaderTarget, place: HeaderMenuPlace, only?: string): boolean;
  /** Closes the menu; a preview it made ends in nothing. */
  closeMenu(): void;
  /** A click on a header's part at a canvas point (the input's): true when taken. */
  click(x: number, y: number, touch: boolean): boolean;
  /** A secondary click at a canvas point: true when the menu took it (the browser's menu then prevented). */
  context(x: number, y: number, e?: { readonly shiftKey?: boolean }): boolean;
  composed(viewer: SeqViewer): void;
  /** Resolves once its lazy chunks (the late parts, the icons, the gestures) landed: what `settle` waits for. */
  loaded(): Promise<void>;
  /** A frame was drawn (the viewer's frame, each frame): the geometry and the parts follow it; a hover makes nothing. */
  framed(): void;
  attach(container: HTMLElement, viewer: SeqViewer): () => void;
  dispose(): void;
}

/** The command that opens the menu at the keyboard cursor (the menu key, Shift+F10; Enter on a header). */
export const HeadersMenuCommand = 'headers.menu';

/** How long a secondary click following the menu key is the same ask (EC-156), ms. */
const AskMs = 1000;
/** The layer types a header's "Color ▸" reads schemes for (their color channel's subjects). */
const Colored = new Set(['letters', 'ranges', 'values', 'logo', 'variants']);

export function headersBehavior(context: BehaviorContext, params: HeadersParams = {}): HeadersBehavior {
  const contextClick = params.contextClick ?? 'everywhere';
  const dom = typeof params.parts === 'function' ? params.parts : undefined;
  const mac = isMac();
  let viewer: SeqViewer | undefined;
  let disposed = false;
  /**
   * The lazy chunks asked for (`loaded`, what `settle` waits for; three at most). A chunk landing after the viewer is
   * retired is ignored, its failure too (a torn-down test environment); before, a failure is the host's to see.
   */
  const loading: Promise<void>[] = [];
  const lazy = <M>(load: Promise<M>, landed: (m: M) => void): void => {
    loading.push(
      load.then(
        (m) => void (disposed || landed(m)),
        (e) => {
          if (!disposed) throw e;
        },
      ),
    );
  };
  /** The registry's entries and a host's, in order; a new `partsOf` each change. */
  let entries: readonly HeaderEntry[] = [];
  let partsOf: HeaderPartsOf = () => undefined;
  /** The layout's facts of each header by its address (a track's or a section's), from the last layout. */
  const lastInput = new Map<string, HeaderInputFacts>();
  const geometry = new State<readonly HeaderGeometry[]>([], sameGeometry);
  const shownMenu = new State<ShownHeaderMenu | undefined>(undefined);
  const parts = new State<HeaderPartsState | undefined>(undefined);
  let menuPort: HeaderMenuPort | undefined;
  let menuPortMounted = false;
  let previewed = false;
  /** When the keyboard (or a long press) last opened the menu: the secondary click after it is the same ask. */
  let askedAt = Number.NEGATIVE_INFINITY;
  const stops: (() => void)[] = [];

  const commandOf = (id: string | undefined) => (id ? viewer?.commands.find((c) => c.id === id) : undefined);

  /** Checks entries against those kept and the viewer's commands: every conflict in one error. */
  const check = (adding: readonly HeaderEntry[], kept: readonly HeaderEntry[]) => {
    const problems: string[] = [];
    const names = new Set(kept.map((e) => e.name));
    for (const e of adding) {
      if (!/^[\w-]+:[\w.-]+$/.test(e.name)) problems.push(`${e.name}: not namespaced (<id>:<name>)`);
      if (names.has(e.name)) problems.push(`${e.name}: registered twice`);
      names.add(e.name);
      if (e.command !== undefined && viewer && !commandOf(e.command))
        problems.push(`${e.name}: no command ${e.command} (seq.viewer.header.command)`);
    }
    if (problems.length) throw new Error(`The header entries are refused: ${problems.join('; ')}`);
  };
  /**
   * The entries' icons as paths: `@molstar/ui`'s Lucide markup and its reader are a lazy chunk (`header-icons.ts`),
   * loaded once an entry has an icon; until then (and for an icon no path reads, `seq.viewer.header.icon`) a slot draws
   * its label's first letter. Loaded, the headers are laid out again.
   */
  let iconOf: ((name: string | undefined, markup?: string) => HeaderIcon | undefined) | undefined;
  let iconsAsked = false;
  const setEntries = (next: readonly HeaderEntry[]) => {
    entries = next;
    if (!iconOf && !iconsAsked && next.some((e) => e.icon !== undefined || e.iconMarkup !== undefined)) {
      iconsAsked = true;
      lazy(import('@molstar/seqx-viewer/header-icons'), (m) => {
        iconOf = m.headerIcon;
        setEntries(entries);
      });
    }
    // a new function: part of the layout's key, so the next frame lays the headers out again
    partsOf = (f) => layoutParts(f);
    viewer?.schedule();
  };

  // ---- facts ----

  const sectionOfAddress = (address: string): ResolvedSection | undefined =>
    context.view.resolved.value?.sections.find((s) => s.address === address);
  const trackOf = (section: ResolvedSection, address: string | undefined): ResolvedTrack | undefined =>
    address === undefined ? undefined : section.tracks.find((t) => t.address === address);

  /** A header's facts: from the resolved view, the spec and the layout's own facts of it (filtered rows). */
  const factsOf = (target: HeaderTarget, input?: HeaderInputFacts): HeaderFacts | undefined => {
    const section = sectionOfAddress(target.section);
    if (!section) return undefined;
    const track = trackOf(section, target.track);
    if (target.track !== undefined && !track) return undefined;
    const model = context.spec.model.value;
    const at = input ?? lastInput.get(target.track ?? target.section);
    const legend = viewer?.behaviors.get('legend') as ThemingBehavior | undefined;
    const specNode = resolveAddress(context.spec.root.value, target.section.split('/').filter(Boolean));
    const layers: HeaderLayerFacts[] = [];
    for (const l of track?.layers ?? []) {
      const address = l.address ?? `${track!.address}:${l.index}`;
      const filtered = at?.filtered.find((x) => x.layer === address)?.count ?? 0;
      // a layer has no label of its own: its id, else its type
      const label = l.address?.split('/').pop() ?? l.type;
      layers.push({
        address,
        type: l.type,
        label,
        hidden: l.params.hidden === true,
        filtered,
        ...(l.params.color !== undefined ? { color: l.params.color } : {}),
        schemes: l.address && legend && Colored.has(l.type) ? legend.fittingSchemes(l.address) : [],
      });
    }
    const label = track
      ? track.label
      : typeof section.params.label === 'string'
        ? section.params.label
        : section.address.split('/').pop()!;
    return {
      on: target.element ? 'element' : track ? 'track' : 'section',
      address: track?.address ?? section.address,
      label,
      section: section.address,
      system: section.axis.system,
      linear: !section.axis.cyclic,
      tracks: (specNode?.children ?? []).flatMap((c) => (c.kind === 'track' && c.id !== undefined ? [c.id] : [])),
      layers,
      hiddenTracks: at?.hiddenTracks ?? hiddenTracksOf(model, section),
      ...(target.element ? { element: target.element } : {}),
    };
  };

  /** What applies to a header: the entries on its kind whose `when` holds. */
  const applying = (f: HeaderFacts, on: HeaderEntry['on']) =>
    entries.filter((e) => e.on === on && (e.when === undefined || e.when(f)));
  const bySlot = (a: HeaderEntry, b: HeaderEntry) =>
    (a.order ?? Number.MAX_SAFE_INTEGER) - (b.order ?? Number.MAX_SAFE_INTEGER) || (a.name < b.name ? -1 : 1);
  const slotted = (e: HeaderEntry) =>
    (e.kind ?? 'action') === 'action' && (e.slot ?? (e.icon ?? commandOf(e.command)?.icon) !== undefined);
  const argsOf = (e: HeaderEntry, f: HeaderFacts) => e.args?.(f) ?? { address: f.address };
  const whyNot = (command: string | undefined, args: Readonly<Record<string, unknown>>): string | undefined => {
    if (command === undefined) return undefined;
    const c = commandOf(command);
    if (!c) return `No command ${command}`;
    const can = c.when?.(args);
    return can === undefined || can === true ? undefined : typeof can === 'string' ? can : 'Not available';
  };
  const keyOf = (command: string | undefined): string | undefined => {
    if (command === undefined) return undefined;
    const listed = viewer?.config.view.keys[command];
    const chord = listed?.length ? parseChord(listed[0]) : commandOf(command)?.shortcut;
    return chord ? formatChord(chord, mac) : undefined;
  };

  /** An entry's menu item for a header's facts. */
  const itemOf = (e: HeaderEntry, f: HeaderFacts, slot = false): HeaderMenuItem => {
    const command = commandOf(e.command);
    const args = argsOf(e, f);
    const subs = e.items?.(f);
    // an action with a submenu opens it: enabled when any of its items runs
    const why = subs?.length && (e.kind ?? 'action') === 'action' ? undefined : whyNot(e.command, args);
    const label = typeof e.label === 'function' ? e.label(f) : (e.label ?? command?.title ?? e.name);
    const text = e.text?.(f);
    const icon = e.icon ?? command?.icon;
    return {
      name: e.name,
      kind: e.kind ?? 'action',
      label,
      ...(icon ? { icon } : {}),
      ...(text !== undefined ? { text } : {}),
      ...(keyOf(e.command) ? { key: keyOf(e.command)! } : {}),
      enabled: why === undefined,
      ...(why !== undefined ? { why } : {}),
      ...(slot ? { slot } : {}),
      ...(e.command !== undefined ? { command: e.command, args } : {}),
      ...(subs?.length
        ? {
            items: subs.map((s): HeaderMenuItem => {
              const w = whyNot(s.command, s.args);
              return {
                name: e.name,
                kind: 'action',
                label: s.label,
                enabled: w === undefined,
                ...(w !== undefined ? { why: w } : {}),
                ...(s.checked !== undefined ? { checked: s.checked } : {}),
                command: s.command,
                args: s.args,
              };
            }),
          }
        : {}),
    };
  };

  /** A header's items: its slotted actions first (in slot order), then its other actions, then its states. */
  const itemsOf = (f: HeaderFacts, on: HeaderEntry['on']): HeaderMenuItem[] => {
    const all = applying(f, on).sort(bySlot);
    const n = viewer ? Math.max(0, actionColumns(viewer.context.renderConfig).n) : 0;
    const slots = all.filter(slotted).slice(0, n);
    return [
      ...slots.map((e) => itemOf(e, f, true)),
      ...all.filter((e) => !slots.includes(e) && (e.kind ?? 'action') === 'action').map((e) => itemOf(e, f)),
      ...all.filter((e) => e.kind === 'state').map((e) => itemOf(e, f)),
    ];
  };

  /** The layout's parts of a header (`HeaderPartsOf`): its slots (canvas only), its states, whether it has a menu. */
  const layoutParts = (input: HeaderInputFacts): HeaderParts | undefined => {
    const target: HeaderTarget = {
      section: input.section.address,
      ...(input.track ? { track: input.track.address } : {}),
    };
    lastInput.set(input.track?.address ?? input.section.address, input);
    const f = factsOf(target, input);
    if (!f) return undefined;
    const items = itemsOf(f, input.track ? 'track' : 'section');
    let slots: HeaderSlotPart[] | undefined;
    let states: HeaderStatePart[] | undefined;
    for (const it of items) {
      if (it.kind === 'state') {
        if (it.text === undefined) continue;
        const n = Number.parseInt(it.text, 10);
        const e = entries.find((x) => x.name === it.name);
        (states ??= []).push({
          name: it.name,
          text: it.text,
          ...(Number.isFinite(n) ? { count: n } : {}),
          // what the state is about when not the header's own target (the filtered layer)
          ...(e?.command === 'layer.filter' && typeof it.args?.address === 'string'
            ? { address: it.args.address }
            : {}),
        });
      } else if (it.slot && !dom) {
        const icon = iconOf?.(it.icon, entries.find((x) => x.name === it.name)?.iconMarkup);
        (slots ??= []).push({ name: it.name, label: it.label, ...(icon ? { icon } : {}), enabled: it.enabled });
      }
    }
    const menu = !dom && items.some((i) => i.kind === 'action' || i.command !== undefined);
    if (!slots && !states && !menu) return undefined;
    return { ...(slots ? { slots } : {}), ...(states ? { states } : {}), ...(menu ? { menu } : {}) };
  };

  // ---- the menu ----

  const menuOf = (target: HeaderTarget, only?: string) => extras?.menuOf(target, only);

  const port = (): HeaderMenuPort | undefined => {
    if (!menuPort && params.menu && !menuPortMounted) menuPort = params.menu(undefined);
    return menuPort;
  };
  const closeMenu = () => {
    if (previewed) {
      previewed = false;
      (viewer?.behaviors.get('legend') as ThemingBehavior | undefined)?.cancelPreview();
    }
    if (!shownMenu.value) return;
    shownMenu.set(undefined);
    menuPort?.hide();
  };
  const openMenu = (target: HeaderTarget, place: HeaderMenuPlace, only?: string): boolean => {
    if (disposed || !viewer) return false;
    const menu = menuOf(target, only);
    if (!menu) return false;
    // the card goes: the menu is what the pointer reads now
    (viewer.behaviors.get('hoverCard') as HoverCardBehavior | undefined)?.interrupt();
    shownMenu.set({ menu, place });
    port()?.show(menu, place, behavior);
    return true;
  };

  const runItem = (item: HeaderMenuItem): string | Promise<string> | undefined => {
    if (!viewer || item.command === undefined) return undefined;
    // a preview ends in the choice's patch (the command cancels it first), not in nothing
    previewed = false;
    return viewer.execute(item.command, item.args ?? {});
  };

  // ---- geometry ----

  const partsState = (): HeaderPartsState => extras?.partsState() ?? { headers: [] };

  // the last layout published and the headers whose actions showed (the parts: a hover moves only these)
  let publishedFor: unknown;
  let showingKey = '';
  const published = () => {
    if (!viewer) return;
    const layout = viewer.frame?.layout;
    if (!layout) return;
    const relaid = layout !== publishedFor;
    if (relaid && extras) {
      publishedFor = layout;
      geometry.set(extras.geometryOf(layout));
    }
    if (!partsPort) return;
    // the parts follow the hover only when the header whose actions show changed (a key of its address)
    const c = viewer.selection.cursor.value;
    const key = `${viewer.hover.current.layer ?? ''}|${c?.track ?? ''}`;
    if (!relaid && key === showingKey) return;
    showingKey = key;
    const next = partsState();
    parts.set(next);
    partsPort.update(next);
  };
  let partsPort: HeaderPartsPort | undefined;

  // ---- a click ----

  const click = (x: number, y: number, touch: boolean): boolean => {
    if (!viewer) return false;
    const hit = viewer.hit(x, y, touch);
    const t = hit?.target;
    if (!hit || (t !== 'header-slot' && t !== 'header-state' && t !== 'header-label' && t !== 'section-header'))
      return false;
    // a header's parts are never a selection's: taken, run by the late parts; before they land a chip's click is the
    // input's (the "n filtered" chip clears its filter in one step), the rest taken and doing nothing
    if (!extras) return t !== 'header-state';
    extras.click(hit, x, y, touch);
    return true;
  };

  const run = (name: string, target: HeaderTarget): string | Promise<string> | undefined => {
    const e = entries.find((x) => x.name === name);
    const f = factsOf(target);
    if (!e || !f || !viewer || e.command === undefined) return undefined;
    return viewer.execute(e.command, argsOf(e, f));
  };

  const contextAt = (x: number, y: number, e: { readonly shiftKey?: boolean } = {}): boolean => {
    // Shift+secondary click is always the browser's; the menu key's own secondary click is the same ask (EC-156)
    if (e.shiftKey || contextClick === 'none' || !viewer) return false;
    if (DomHost.now() - askedAt < AskMs) {
      askedAt = Number.NEGATIVE_INFINITY;
      return true;
    }
    return !!extras?.context(x, y, contextClick === 'header');
  };

  // ---- the menu key ----

  /** The menu at the keyboard cursor (the late parts'): the browser's secondary click after it is this ask (EC-156). */
  const menuAtCursor = (): CommandResult => {
    const ok = !!extras?.menuAtCursor();
    if (ok) askedAt = DomHost.now();
    return { message: ok ? 'Actions' : '' };
  };
  const commands: CommandDefinition<unknown>[] = [
    {
      id: HeadersMenuCommand,
      title: 'Open the actions of the track at the cursor',
      category: 'Headers',
      scope: 'canvas',
      when: () => (disposed ? 'The viewer is retired' : extras?.cursorTarget() ? true : 'No cursor'),
      run: () => menuAtCursor(),
    },
  ];

  /**
   * The late parts (`header-extras.ts`, a lazy chunk asked for once composed): the menu's data, the geometry and the
   * parts, what a secondary click opens, the words, the cursor on headers.
   */
  let extras: HeaderExtras | undefined;

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
    facts: (target) => factsOf(target),
    menuOf: (target) => menuOf(target),
    partsState,
    run,
    runItem,
    preview(item) {
      const legend = viewer?.behaviors.get('legend') as ThemingBehavior | undefined;
      const a = item.args;
      if (!legend || item.command !== 'layer.channel' || typeof a?.address !== 'string') return false;
      previewed = legend.preview(a.address, typeof a.channel === 'string' ? a.channel : 'color', a.value as never);
      return previewed;
    },
    openMenu,
    closeMenu,
    framed() {
      extras?.keepCursor();
      published();
    },

    click,
    context: contextAt,
    async loaded() {
      await Promise.all(loading);
    },
    composed(v) {
      viewer = v;
      const listed = (v.registry ?? []).flatMap((e) =>
        e.kind === 'header-action' ? [(e as HeaderActionEntry).entry] : [],
      );
      check(listed, []);
      setEntries(listed);
      lazy(import('./header-extras.js'), (m) => {
        extras = m.headerExtras({
          viewer: v,
          openMenu,
          run,
          geometry: () => geometry.value,
          hiddenChipMenu: params.hiddenChip === 'menu',
          dom: !!dom,
          entries: () => entries,
          factsOf: (t) => factsOf(t),
          itemOf,
          itemsOf,
          menuAtCursor,
        });
        describerFor(context).headers = extras.words;
        publishedFor = undefined;
        published();
      });
      // a host's DOM parts take the focus when the cursor comes onto a header
      stops.push(
        context.selection.cursor.listen((c) => {
          if (c?.keyboard && c.track !== undefined && partsPort?.focus)
            partsPort.focus(c.track === c.section ? { section: c.section } : { section: c.section, track: c.track });
        }),
      );
      // a new layout, an undo or a load: the menu was of the old picture
      stops.push(context.view.resolved.listen(() => closeMenu()));
    },
    attach(container, v) {
      // the menu's port over the canvas (a headless one given up first)
      if (params.menu) {
        if (menuPort && menuPort.disposeWithViewer !== false) menuPort.dispose?.();
        menuPort = params.menu(container);
        menuPortMounted = true;
      }
      if (dom) {
        partsPort = dom(container, behavior);
        showingKey = '\u0000';
        publishedFor = undefined;
        published();
      }
      // the secondary click (Shift: the browser's) and the long press (touch's context menu, once per ask)
      // (a lazy chunk, loaded at the mount)
      let stopGestures: (() => void) | undefined;
      let detached = false;
      lazy(import('@molstar/seqx-viewer/header-gestures'), (m) => {
        if (detached) return;
        stopGestures = m.headerGestures(container, params.longPressMs ?? 500, contextAt, (x, y) => {
          if (!shownMenu.value && contextAt(x, y)) askedAt = DomHost.now();
        });
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
      };
    },
    dispose() {
      disposed = true;
      closeMenu();
      for (const s of stops.splice(0)) s();
      if (menuPort && menuPort.disposeWithViewer !== false) menuPort.dispose?.();
      menuPort = undefined;
      const d = describerFor(context);
      if (extras && d.headers === extras.words) d.headers = undefined;
    },
  };
  return behavior;
}

/** The `headers` behavior (_spec/README.md "Behaviors", Headers): slots, chips, the menu and the parts of the headers. */
export const HeadersProvider: BehaviorProvider<HeadersParams | undefined, HeadersBehavior> = {
  name: 'headers',
  create: (context, params) => headersBehavior(context, params ?? {}),
};

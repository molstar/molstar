/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { resolveAddress } from '@molstar/core/repr/tree';
import {
  actionColumns,
  type Box,
  type HeaderInputFacts,
  type HeaderLine,
  HeaderMenuName,
  type HeaderParts,
  type HeaderPartsOf,
  type HeaderSlotPart,
  type HeaderStatePart,
  type Hit as RenderHit,
  type SectionLayout,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import { actionsShown } from '@molstar/seqx-core/render/headers-draw';
import { hiddenTracksOf } from '@molstar/seqx-core/render/headers';
import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import { formatChord, isMac } from '@molstar/ui/commands/keys';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom';
import { parseChord } from '@molstar/seqx-viewer/bindings';
import { describerFor, type Place } from '@molstar/seqx-viewer/describe';
import type {
  HeaderActionEntry,
  HeaderEntry,
  HeaderFacts,
  HeaderGeometry,
  HeaderLayerFacts,
  HeaderMenu,
  HeaderMenuItem,
  HeaderPartsHeader,
  HeaderPartsState,
  HeaderTarget,
} from '@molstar/seqx-viewer/header-entries';
import { headerIcon } from '@molstar/seqx-viewer/header-icons';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { HoverCardBehavior } from './hover-card.js';
import type { ThemingBehavior } from './theming.js';

// The `headers` behavior (_spec/README.md "Behaviors", Headers): the header entries (the registry's `header-action`
// entries and a host's at run time, with their undo) made into the layout's parts once a layout (`partsOf`, the frame's
// `headers`), a slot's or a chip's click run, the menu (`headerMenu`, shown by its port), the secondary click under
// `contextClick` (Shift+secondary click always the browser's) and the long press, the menu key once per ask (EC-156),
// the geometry (`viewer.headers`) and a host's DOM parts (`headerParts`), each published on change only. A hover redraws
// the actions with no layout and nothing published (the geometry is the layout's, compared by identity first).

/** Where a menu is shown: a point in the canvas (CSS px from its top left). */
export interface HeaderMenuPlace {
  readonly x: number;
  readonly y: number;
}

/** The menu's port (_spec/README.md "Ports", Header menu): shows a menu at a place, one at a time. */
export interface HeaderMenuPort {
  show(menu: HeaderMenu, place: HeaderMenuPlace, behavior: HeadersBehavior): void;
  hide(): void;
  dispose?(): void;
  /** false: the host owns it and the viewer never disposes it. */
  readonly disposeWithViewer?: boolean;
}

/** Made once a viewer: with the canvas's element when mounted, else none (headless). */
export type HeaderMenuPortFactory = (container: HTMLElement | undefined) => HeaderMenuPort;

/** A host's DOM parts over the headers' boxes (_spec/README.md "Ports", Header parts). */
export interface HeaderPartsPort {
  update(parts: HeaderPartsState): void;
  /** The keyboard cursor came onto a header: its parts take the focus. */
  focus?(target: HeaderTarget): void;
  dispose?(): void;
  readonly disposeWithViewer?: boolean;
}

/** Made at the mount, with the canvas's element and the behavior its presses run through. */
export type HeaderPartsPortFactory = (container: HTMLElement, behavior: HeadersBehavior) => HeaderPartsPort;

export interface HeadersParams {
  /** The menu's adapter (the default spec gives `@molstar/ui`'s context menu; none: the menu's state moves alone). */
  readonly menu?: HeaderMenuPortFactory;
  /** `canvas` (default): slots drawn; a factory: a host's DOM parts over the boxes (the canvas draws the states). */
  readonly parts?: 'canvas' | HeaderPartsPortFactory;
  /** What a secondary click opens: `everywhere` (default), `header` (the label column), `none` (the browser's). */
  readonly contextClick?: 'everywhere' | 'header' | 'none';
}

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
  attach(container: HTMLElement, viewer: SeqViewer): () => void;
  dispose(): void;
}

/** The command that opens the menu at the keyboard cursor (the menu key, Shift+F10; Enter on a header). */
export const HeadersMenuCommand = 'headers.menu';

/** How long a secondary click following the menu key is the same ask (EC-156), ms. */
const AskMs = 1000;
/** A long press's slop, CSS px. */
const PressSlop = 8;
/** The layer types a header's "Color ▸" reads schemes for (their color channel's subjects). */
const Colored = new Set(['letters', 'ranges', 'values', 'logo', 'variants']);

const sameBox = (a: Box | undefined, b: Box | undefined) =>
  a === b || (!!a && !!b && a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height);

function sameGeometry(a: readonly HeaderGeometry[], b: readonly HeaderGeometry[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (
      x.target.section !== y.target.section ||
      x.target.track !== y.target.track ||
      x.shown !== y.shown ||
      !sameBox(x.line, y.line) ||
      !sameBox(x.label, y.label) ||
      !sameBox(x.well, y.well) ||
      x.slots.length !== y.slots.length ||
      x.states.length !== y.states.length
    )
      return false;
    for (let k = 0; k < x.slots.length; k++)
      if (x.slots[k].name !== y.slots[k].name || !sameBox(x.slots[k].box, y.slots[k].box)) return false;
    for (let k = 0; k < x.states.length; k++)
      if (
        x.states[k].name !== y.states[k].name ||
        x.states[k].text !== y.states[k].text ||
        !sameBox(x.states[k].box, y.states[k].box)
      )
        return false;
  }
  return true;
}

export function headersBehavior(context: BehaviorContext, params: HeadersParams = {}): HeadersBehavior {
  const contextClick = params.contextClick ?? 'everywhere';
  const dom = typeof params.parts === 'function' ? params.parts : undefined;
  const mac = isMac();
  let viewer: SeqViewer | undefined;
  let disposed = false;
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
      if (e.icon !== undefined && !headerIcon(e.icon))
        problems.push(`${e.name}: icon ${e.icon} not readable as paths (seq.viewer.header.icon)`);
    }
    if (problems.length) throw new Error(`The header entries are refused: ${problems.join('; ')}`);
  };
  const setEntries = (next: readonly HeaderEntry[]) => {
    entries = next;
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
      const label = typeof l.params.label === 'string' ? l.params.label : (l.address?.split('/').pop() ?? l.type);
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
        const icon = it.icon ? headerIcon(it.icon) : undefined;
        (slots ??= []).push({ name: it.name, label: it.label, ...(icon ? { icon } : {}), enabled: it.enabled });
      }
    }
    const menu = !dom && items.some((i) => i.kind === 'action' || i.command !== undefined);
    if (!slots && !states && !menu) return undefined;
    return { ...(slots ? { slots } : {}), ...(states ? { states } : {}), ...(menu ? { menu } : {}) };
  };

  // ---- the menu ----

  const menuOf = (target: HeaderTarget, only?: string): HeaderMenu | undefined => {
    const groups: HeaderMenuItem[][] = [];
    if (only !== undefined) {
      const f = factsOf(target);
      const e = entries.find((x) => x.name === only);
      if (!f || !e) return undefined;
      const it = itemOf(e, f);
      // a chip's list: each thing it hides, then the entry itself ("Show all"); a slot's submenu: its items
      const own =
        it.kind === 'state' ? [...(it.items ?? []), { ...it, items: undefined } as HeaderMenuItem] : (it.items ?? [it]);
      return { target, groups: [own.map(({ items: _, ...x }) => x)] };
    }
    if (target.element) {
      const f = factsOf(target);
      if (f) groups.push(itemsOf(f, 'element'));
    }
    if (target.track !== undefined) {
      const f = factsOf({ section: target.section, track: target.track });
      if (f) groups.push(itemsOf(f, 'track'));
    }
    const sf = factsOf({ section: target.section });
    if (sf) groups.push(itemsOf(sf, 'section'));
    const out = groups.filter((g) => g.length);
    return out.length ? { target, groups: out } : undefined;
  };

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

  /** A section's layout by address in the last frame. */
  const sectionLayout = (layout: ViewLayout | undefined, address: string): SectionLayout | undefined =>
    layout?.sections.find((s) => s.section.address === address);

  /** A header line's slot boxes: as laid out (canvas), or the action columns' a host's DOM parts fill. */
  const slotBoxes = (line: HeaderLine, items: readonly HeaderMenuItem[]) => {
    if (!dom) {
      const out = line.slots.map((s) => ({
        name: s.part.name,
        box: { x: s.x, y: s.y, width: s.width, height: s.height },
      }));
      if (line.menu)
        out.push({
          name: HeaderMenuName,
          box: { x: line.menu.x, y: line.menu.y, width: line.menu.width, height: line.menu.height },
        });
      return out;
    }
    const cols = actionColumns(viewer!.context.renderConfig);
    if (cols.n < 0 || line.short) return [];
    const px = viewer!.context.renderConfig.header.slotPx;
    const cy = line.label.y + line.label.height / 2;
    const y = cy - px / 2;
    const out = items
      .filter((i) => i.slot)
      .slice(0, cols.n)
      .map((i, k) => ({ name: i.name, box: { x: cols.x + k * px, y, width: px, height: px } }));
    if (items.some((i) => i.kind === 'action'))
      out.push({ name: HeaderMenuName, box: { x: cols.x + cols.n * px, y, width: px, height: px } });
    return out;
  };

  const geometryOf = (layout: ViewLayout): HeaderGeometry[] => {
    const out: HeaderGeometry[] = [];
    const h = viewer?.height ?? Number.POSITIVE_INFINITY;
    for (const s of layout.sections) {
      const lines: [HeaderLine, HeaderTarget, Box | undefined][] = [];
      if (s.header) lines.push([s.header, { section: s.section.address }, undefined]);
      const seen = new Set<string>();
      for (const t of s.tracks) {
        if (!t.header || seen.has(t.track.address)) continue;
        seen.add(t.track.address);
        lines.push([t.header, { section: s.section.address, track: t.track.address }, t.well]);
      }
      for (const [line, target, well] of lines) {
        const f = dom ? factsOf(target) : undefined;
        const items = f ? itemsOf(f, target.track ? 'track' : 'section') : [];
        const b = line.line;
        out.push({
          target,
          line: b,
          label: line.label,
          slots: slotBoxes(line, items),
          states: line.states.map((st) => ({ name: st.name, box: st.box, text: st.text, mark: st.mark })),
          ...(well ? { well } : {}),
          shown: b.y + b.height > 0 && b.y < h && b.y >= s.y - 1 && b.y + b.height <= s.y + s.height + 1,
        });
      }
    }
    return out;
  };

  /** The hovered or cursor's header whose actions show (`header.actions.show`), for the parts. */
  const showing = (g: HeaderGeometry): boolean => {
    if (!viewer) return false;
    const s = sectionLayout(viewer.frame?.layout, g.target.section);
    return (
      !!s &&
      actionsShown(
        g.target.track ?? g.target.section,
        s.id,
        viewer.emphasis,
        viewer.context.renderConfig,
        g.target.track !== undefined,
      )
    );
  };

  const partsState = (): HeaderPartsState => {
    const layout = viewer?.frame?.layout;
    const headers: HeaderPartsHeader[] = [];
    for (const g of layout ? geometryOf(layout) : []) {
      if (!g.shown) continue;
      const f = factsOf(g.target);
      if (!f) continue;
      headers.push({ geometry: g, actions: showing(g), items: itemsOf(f, g.target.track ? 'track' : 'section') });
    }
    return { headers };
  };

  // the last layout published and the headers whose actions showed (the parts: a hover moves only these)
  let publishedFor: unknown;
  let showingKey = '';
  const published = () => {
    if (!viewer) return;
    const layout = viewer.frame?.layout;
    if (!layout) return;
    const relaid = layout !== publishedFor;
    if (relaid) {
      publishedFor = layout;
      geometry.set(geometryOf(layout));
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

  // ---- the target under a point ----

  /** The header a hit is on (a track's line, a slot, a chip, a section's line): its target. */
  const headerTarget = (hit: RenderHit, layout: ViewLayout | undefined): HeaderTarget | undefined => {
    for (const s of layout?.sections ?? []) {
      if (s.section.address === hit.layer) return { section: s.section.address };
      for (const t of s.tracks)
        if (t.track.address === hit.layer) return { section: s.section.address, track: t.track.address };
    }
    return undefined;
  };
  /** What a secondary click at a point opens: an element's, a track's or a section's menu. */
  const targetAt = (x: number, y: number, hit: RenderHit | undefined): HeaderTarget | undefined => {
    const layout = viewer?.frame?.layout;
    if (!layout || !viewer) return undefined;
    if (hit) {
      const header = headerTarget(hit, layout);
      if (header) return header;
      for (const s of layout.sections)
        for (const t of s.tracks)
          if (t.layers.some((l) => l.key === hit.layer))
            return {
              section: s.section.address,
              track: t.track.address,
              ...(hit.element ? { element: { layer: hit.layer, key: hit.element.key } } : {}),
            };
    }
    // a track's band under the point, else its section
    const section = viewer.input.sectionAt(y);
    const s = section === undefined ? undefined : sectionLayout(layout, section);
    if (!s) return undefined;
    const t = s.tracks.find((t) => y >= t.box.y && y < t.box.y + t.box.height);
    return t ? { section: s.section.address, track: t.track.address } : { section: s.section.address };
  };

  /** Where a menu opens for a header: under its line's label start, else at the point. */
  const placeOf = (target: HeaderTarget): HeaderMenuPlace => {
    const g = geometry.value.find((x) => x.target.section === target.section && x.target.track === target.track);
    return g ? { x: g.label.x, y: g.line.y + g.line.height } : { x: 0, y: 0 };
  };

  const click = (x: number, y: number, touch: boolean): boolean => {
    if (!viewer) return false;
    const hit = viewer.hit(x, y, touch);
    const t = hit?.target;
    if (!hit || (t !== 'header-slot' && t !== 'header-state' && t !== 'header-label' && t !== 'section-header'))
      return false;
    const target = headerTarget(hit, viewer.frame?.layout);
    if (!target) return true;
    if (t === 'header-label' || t === 'section-header') {
      // no hover on touch: a tap on the label shows its slots (as the pointer resting there would), a tap on a slot
      // then runs it
      if (touch) viewer.input.move(x, y, { pointerType: 'touch' });
      return true;
    }
    const name = hit.state ?? '';
    const g = geometry.value.find((g) => g.target.track === target.track && g.target.section === target.section);
    const box = g?.slots.find((s) => s.name === name)?.box ?? g?.states.find((s) => s.name === name)?.box;
    const at = box ? { x: box.x, y: box.y + box.height } : { x, y };
    if (name === HeaderMenuName) return openMenu(target, at), true;
    const e = entries.find((e) => e.name === name);
    if (!e) return true;
    // an action with a submenu opens it at its slot; a chip lists what it hides under `header.hiddenChip.click: menu`
    const f = factsOf(target);
    const subs = f && e.items?.(f);
    if (subs?.length && (e.kind !== 'state' || viewer.config.header.hiddenChip.click === 'menu'))
      return openMenu(target, at, e.name), true;
    void run(e.name, target);
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
    const hit = viewer.hit(x, y);
    const c = viewer.context.renderConfig;
    const inColumn = x >= c.inset && x < c.inset + c.labelWidth;
    if (contextClick === 'header' && !inColumn) return false;
    const target = targetAt(x, y, hit);
    if (!target) return false;
    // a chip's own list (each thing it hides, then all)
    const only = hit?.target === 'header-state' && hit.state ? hit.state : undefined;
    const known = only !== undefined && entries.find((x) => x.name === only)?.items;
    return openMenu(target, { x, y }, known ? only : undefined);
  };

  // ---- the words ----

  /** The words of a header's part (a slot, a chip, a label): its entry's label and key, why not, what a click does. */
  const words = (place: Place): { title: string; text: string } | undefined => {
    if (!viewer || place.layer === undefined) return undefined;
    const t = place.target;
    const target = headerTarget({ layer: place.layer } as RenderHit, viewer.frame?.layout);
    if (!target) return undefined;
    const f = factsOf(target);
    if (!f) return undefined;
    const items = itemsOf(f, target.track ? 'track' : 'section');
    if (t === 'header-slot') {
      if (place.state === HeaderMenuName) return { title: 'Actions', text: `Actions of ${f.label}` };
      const it = items.find((i) => i.name === place.state);
      if (!it) return undefined;
      const title = it.key ? `${it.label} (${it.key})` : it.label;
      return { title, text: it.why ? `${title}, unavailable: ${it.why}` : title };
    }
    if (t === 'header-state') {
      const it = items.find((i) => i.name === place.state);
      if (!it?.text) return undefined;
      return { title: it.text, text: it.command ? `${it.text}: click ${it.label.toLowerCase()}` : it.text };
    }
    if (t === 'header-label' || t === 'section-header') {
      const states = items.filter((i) => i.kind === 'state' && i.text).map((i) => i.text!);
      const cursor = place.state === 'cursor';
      const layers = f.layers.filter((l) => !l.hidden).map((l) => l.label);
      const said = [f.label, ...states, ...(cursor ? ['Enter for actions'] : [])];
      return { title: f.label, text: (cursor || !layers.length ? said : [...said, ...layers]).join(', ') };
    }
    return undefined;
  };

  // ---- the menu key ----

  const cursorTarget = (): HeaderTarget | undefined => {
    const c = context.selection.cursor.value;
    if (!c) return undefined;
    if (c.track !== undefined)
      return c.track === c.section ? { section: c.section } : { section: c.section, track: c.track };
    const s = sectionOfAddress(c.section);
    const first = s?.tracks[0];
    return s ? { section: s.address, ...(first ? { track: first.address } : {}) } : undefined;
  };

  const commands: CommandDefinition<unknown>[] = [
    {
      id: HeadersMenuCommand,
      title: 'Open the actions of the track at the cursor',
      category: 'Headers',
      scope: 'canvas',
      when: () => (disposed ? 'The viewer is retired' : cursorTarget() ? true : 'No cursor'),
      run: (): CommandResult => {
        const target = cursorTarget()!;
        const ok = openMenu(target, placeOf(target));
        // the browser's own secondary click after the menu key is this ask (EC-156)
        if (ok) askedAt = DomHost.now();
        return { message: ok ? 'Actions' : '' };
      },
    },
  ];

  // ---- the cursor off what was hidden ----

  /** The headers' order the last layout drew, by section: where the cursor goes when its track is hidden. */
  let order = new Map<string, string[]>();
  const keepCursor = () => {
    const c = context.selection.cursor.value;
    const layout = viewer?.frame?.layout;
    if (c?.track !== undefined && c.track !== c.section) {
      const s = sectionOfAddress(c.section);
      if (s && !s.tracks.some((t) => t.address === c.track)) {
        const was = order.get(c.section) ?? [];
        const i = was.indexOf(c.track);
        const there = (a: string) => a === c.section || s.tracks.some((t) => t.address === a);
        const next = was.slice(i + 1).find(there) ?? was.slice(0, Math.max(0, i)).reverse().find(there) ?? c.section;
        context.selection.setCursor({ ...c, track: next });
      }
    }
    if (!layout) return;
    order = new Map();
    for (const s of layout.sections)
      order.set(s.section.address, [
        ...(s.header ? [s.section.address] : []),
        ...s.tracks.filter((t) => t.header?.hit).map((t) => t.track.address),
      ]);
  };

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
    click,
    context: contextAt,
    composed(v) {
      viewer = v;
      const listed = (v.registry ?? []).flatMap((e) =>
        e.kind === 'header-action' ? [(e as HeaderActionEntry).entry] : [],
      );
      check(listed, []);
      setEntries(listed);
      describerFor(context).headers = words;
      stops.push(v.drawn.listen(published));
      stops.push(v.drawn.listen(keepCursor));
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
      // the secondary click (Shift: the browser's) and the long press (touch's context menu)
      const local = (e: MouseEvent) => {
        const r = container.getBoundingClientRect();
        return [e.clientX - r.left, e.clientY - r.top] as const;
      };
      const onContext = (e: MouseEvent) => {
        const [x, y] = local(e);
        if (contextAt(x, y, e)) e.preventDefault();
      };
      let press: { x: number; y: number; timer: unknown } | undefined;
      const endPress = () => {
        if (press) DomHost.clearTimeout(press.timer as never);
        press = undefined;
      };
      const onDown = (e: PointerEvent) => {
        endPress();
        if (e.pointerType !== 'touch') return;
        const [x, y] = local(e);
        press = {
          x,
          y,
          timer: DomHost.setTimeout(() => {
            press = undefined;
            if (shownMenu.value) return;
            if (contextAt(x, y)) askedAt = DomHost.now();
          }, v.config.header.menu.longPressMs),
        };
      };
      const onMove = (e: PointerEvent) => {
        if (!press) return;
        const [x, y] = local(e);
        if (Math.max(Math.abs(x - press.x), Math.abs(y - press.y)) > PressSlop) endPress();
      };
      const listeners: [string, (e: any) => void][] = [
        ['contextmenu', onContext],
        ['pointerdown', onDown],
        ['pointermove', onMove],
        ['pointerup', endPress],
        ['pointercancel', endPress],
      ];
      for (const [type, fn] of listeners) container.addEventListener(type, fn);
      return () => {
        endPress();
        for (const [type, fn] of listeners) container.removeEventListener(type, fn);
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
      if (viewer) {
        const d = describerFor(context);
        if (d.headers === words) d.headers = undefined;
      }
    },
  };
  return behavior;
}

/** The `headers` behavior (_spec/README.md "Behaviors", Headers): slots, chips, the menu and the parts of the headers. */
export const HeadersProvider: BehaviorProvider<HeadersParams | undefined, HeadersBehavior> = {
  name: 'headers',
  create: (context, params) => headersBehavior(context, params ?? {}),
};

/** The menu's data for a target (DOM-free; _spec/README.md "Behaviors", Headers): undefined without the behavior. */
export function headerMenu(viewer: SeqViewer, target: HeaderTarget): HeaderMenu | undefined {
  return (viewer.behaviors.get('headers') as HeadersBehavior | undefined)?.menuOf(target);
}

/** Every shown header's parts (DOM-free): what a parts port is given. */
export function headerParts(viewer: SeqViewer): HeaderPartsState | undefined {
  return (viewer.behaviors.get('headers') as HeadersBehavior | undefined)?.partsState();
}

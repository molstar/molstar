/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node as SpecNode, resolveAddress } from '@molstar/core/repr/tree';
import {
  actionColumns,
  BlockPlacement,
  type Box,
  type HeaderIcon,
  type HeaderInputFacts,
  type HeaderLine,
  HeaderMenuName,
  type HeaderParts,
  type HeaderSlotPart,
  type HeaderStatePart,
  type Hit as RenderHit,
  type SectionLayout,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import { hiddenTracksOf } from '@molstar/seqx-core/render/headers';
import { actionsShown } from '@molstar/seqx-core/render/headers-draw';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import { formatChord, isMac } from '@molstar/ui/commands/keys';
import type { CommandResult } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom/host';
import { parseChord } from '@molstar/seqx-viewer/bindings';
import type { Place } from '@molstar/seqx-viewer/describe';
import type {
  HeaderEntry,
  HeaderFacts,
  HeaderGeometry,
  HeaderLayerFacts,
  HeaderMenu,
  HeaderMenuItem,
  HeaderMenuPlace,
  HeaderPartsHeader,
  HeaderPartsState,
  HeaderSubEntry,
  HeaderTarget,
} from '@molstar/seqx-viewer/header-entries';
import { headerCursor } from './header-cursor.js';
import type { HeaderShared } from './headers.js';
import type { HoverCardBehavior } from './hover-card.js';
import type { ThemingBehavior } from './theming.js';

// The `headers` behavior's late parts (_spec/README.md "Behaviors", Headers; "Composition", Injected, opt in: a listed
// part used late loads on first use): a lazy chunk the behavior asks for once composed. The menu's data, its opening and
// a run, the headers' geometry and a host's DOM parts' state, what a click and a secondary click do, the words of a
// header's parts and the cursor on headers. Until it lands the slots and chips draw (the layout's parts are the
// behavior's own) and a click waits for it; then the geometry is published and the rest works. Internal: no entry
// exports it.

/** @internal */
export interface HeaderExtras {
  /** A header's facts now (from the resolved view and the last layout); undefined when it is not shown. */
  factsOf(target: HeaderTarget): HeaderFacts | undefined;
  /** The layout's parts of a header (`HeaderPartsOf`). */
  layoutParts(input: HeaderInputFacts): HeaderParts | undefined;
  menuOf(target: HeaderTarget, only?: string): HeaderMenu | undefined;
  openMenu(target: HeaderTarget, place: HeaderMenuPlace, only?: string): boolean;
  run(name: string, target: HeaderTarget): string | Promise<string> | undefined;
  runItem(item: HeaderMenuItem): string | Promise<string> | undefined;
  preview(item: HeaderMenuItem): boolean;
  geometryOf(layout: ViewLayout): HeaderGeometry[];
  partsState(): HeaderPartsState;
  targetAt(x: number, y: number, hit: RenderHit | undefined): HeaderTarget | undefined;
  words(place: Place): { title: string; text: string } | undefined;
  /** A click on a header's part (a slot runs its entry, a chip its state, a label on touch shows the slots). */
  click(hit: RenderHit, x: number, y: number, touch: boolean): void;
  /** What the menu key opens: the cursor's header, else its section's first track's. */
  cursorTarget(): HeaderTarget | undefined;
  /** Opens the menu at the keyboard cursor (the menu key): the secondary click after it is this ask (EC-156). */
  menuAtCursor(): CommandResult;
  /** A secondary click at a point (`contextClick`; Shift: the browser's): true when the menu took it. */
  context(x: number, y: number, e?: { readonly shiftKey?: boolean }): boolean;
  /** A long press at a point: the menu there, once per ask. */
  pressed(x: number, y: number): void;
  /** A layout's geometry, the one published when no box moved (`viewer.headers` changes only then). */
  published(layout: ViewLayout): readonly HeaderGeometry[];
  /** The frame's: the cursor off a track a layout hid. */
  keepCursor(): void;
}

/** How long a secondary click following the menu key is the same ask (EC-156), ms. */
const AskMs = 1000;
/** The layer types a header's "Color ▸" reads schemes for (their color channel's subjects). */
const Colored = new Set(['letters', 'ranges', 'values', 'logo', 'variants']);

const sectionLayout = (layout: ViewLayout | undefined, address: string): SectionLayout | undefined =>
  layout?.sections.find((s) => s.section.address === address);

/** @internal The header a hit is on (a track's line, a slot, a chip, a section's line): its target. */
export function headerTarget(hit: RenderHit, layout: ViewLayout | undefined): HeaderTarget | undefined {
  for (const s of layout?.sections ?? []) {
    if (s.section.address === hit.layer) return { section: s.section.address };
    for (const t of s.tracks)
      if (t.track.address === hit.layer) return { section: s.section.address, track: t.track.address };
  }
  return undefined;
}

/**
 * A section's block facts for its menu (genome coordinates): the block a menu was opened on, whether the section is of
 * assemblies (its `include` can be set) and has an `include`.
 */
function blockFacts(
  model: StaticModel,
  section: ResolvedSection,
  node: SpecNode | undefined,
  index: number | undefined,
): Pick<HeaderFacts, 'block' | 'assemblies' | 'included'> {
  const of = node?.params?.of;
  const list = (Array.isArray(of) ? of : [of]).filter((x): x is string => typeof x === 'string');
  const assemblies = list.length > 0 && list.every((a) => model.assemblies.has(a));
  const b = index === undefined ? undefined : section.blocks?.[index];
  return {
    ...(b
      ? {
          block: {
            address: b.address,
            label: b.label,
            name: model.systems.get(b.axis.system)?.name ?? b.label,
            reversed: b.axis.dir !== section.axis.dir,
          },
        }
      : {}),
    ...(assemblies ? { assemblies: true } : {}),
    ...(node?.params?.include !== undefined ? { included: true } : {}),
  };
}

/** @internal */
export function headerExtras(inner: HeaderShared): HeaderExtras {
  const { viewer } = inner;
  const on = (t: HeaderTarget) => (t.track ? 'track' : 'section');
  const mac = isMac();
  const commandOf = (id: string | undefined) => (id ? viewer.commands.find((c) => c.id === id) : undefined);

  // ---- facts and the layout's parts ----

  /** The layout's facts of each header by its address (a track's or a section's), from the last layout. */
  const lastInput = new Map<string, HeaderInputFacts>();
  /**
   * The entries' icons as paths: their reader and the built-ins' markup are a lazy chunk (`header-icons.ts`), asked
   * for at the first layout with an icon; until then (and for an icon no path reads, `seq.viewer.header.icon`) a slot
   * draws its label's first letter. Landed, the headers are laid out again.
   */
  let iconOf: ((name: string | undefined, markup?: string) => HeaderIcon | undefined) | undefined;
  let iconsAsked = false;

  const sectionOfAddress = (address: string): ResolvedSection | undefined =>
    inner.context.view.resolved.value?.sections.find((s) => s.address === address);
  const trackOf = (section: ResolvedSection, address: string | undefined): ResolvedTrack | undefined =>
    address === undefined ? undefined : section.tracks.find((t) => t.address === address);

  /** A header's facts: from the resolved view, the spec and the layout's own facts of it (filtered rows). */
  const factsOf = (target: HeaderTarget, input?: HeaderInputFacts): HeaderFacts | undefined => {
    const section = sectionOfAddress(target.section);
    if (!section) return undefined;
    const track = trackOf(section, target.track);
    if (target.track !== undefined && !track) return undefined;
    const model = inner.context.spec.model.value;
    const at = input ?? lastInput.get(target.track ?? target.section);
    const legend = viewer.behaviors.get('legend') as ThemingBehavior | undefined;
    const specNode = resolveAddress(inner.context.spec.root.value, target.section.split('/').filter(Boolean));
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
        ...(l.type === 'rows' && l.params.order === 'tree' && l.params.tree !== undefined ? { tree: true } : {}),
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
      ...(track && layers.some((l) => l.tree)
        ? {
            tree: {
              open: inner.context.view.band(track.address)?.tree === true,
              lacks: viewer.frame?.layout.narrow?.lacks === true,
            },
          }
        : {}),
      ...(target.element ? { element: target.element } : {}),
      viewer,
      ...blockFacts(model, section, specNode, target.block),
    };
  };

  /** What applies to a header: the entries on its kind whose `when` holds, in slot order (`order`, then name). */
  const applying = (f: HeaderFacts, on: HeaderEntry['on']) =>
    inner
      .entries()
      .filter((e) => e.on === on && (e.when === undefined || e.when(f)))
      .sort(
        (a, b) =>
          (a.order ?? Number.MAX_SAFE_INTEGER) - (b.order ?? Number.MAX_SAFE_INTEGER) || (a.name < b.name ? -1 : 1),
      );
  const slotted = (e: HeaderEntry) =>
    (e.kind ?? 'action') === 'action' && (e.slot ?? (e.icon ?? commandOf(e.command)?.icon) !== undefined);
  const argsOf = (e: HeaderEntry, f: HeaderFacts) => e.args?.(f) ?? { address: f.address };
  /**
   * An entry's submenu: its `items`, else its command's choices when the command's one param is a select (its options,
   * each the entry's args with the choice by the param's name; the default checked).
   */
  const subsOf = (e: HeaderEntry, f: HeaderFacts): readonly HeaderSubEntry[] | undefined => {
    if (e.items) return e.items(f);
    const c = commandOf(e.command);
    if (!c?.params) return undefined;
    const args = argsOf(e, f);
    const ps = c.params(args, {});
    const names = Object.keys(ps);
    const p = names.length === 1 ? ps[names[0]] : undefined;
    if (p?.type !== 'select') return undefined;
    return p.options.map(([value, label]) => ({
      label: String(label),
      command: e.command!,
      args: { ...args, [names[0]]: value },
      checked: value === p.defaultValue,
    }));
  };
  const whyNot = (command: string | undefined, args: Readonly<Record<string, unknown>>): string | undefined => {
    if (command === undefined) return undefined;
    const c = commandOf(command);
    if (!c) return `No command ${command}`;
    const can = c.when?.(args);
    return can === undefined || can === true ? undefined : typeof can === 'string' ? can : 'Not available';
  };

  /**
   * The layout's parts of a header (`HeaderPartsOf`): its slots (the first applying slotted actions), its states,
   * whether it has a menu; `hosted` once a host's DOM parts draw them.
   */
  const layoutParts = (input: HeaderInputFacts): HeaderParts | undefined => {
    const target: HeaderTarget = {
      section: input.section.address,
      ...(input.track ? { track: input.track.address } : {}),
    };
    lastInput.set(input.track?.address ?? input.section.address, input);
    const f = factsOf(target, input);
    if (!f) return undefined;
    const all = applying(f, input.track ? 'track' : 'section');
    if (!iconsAsked && inner.entries().some((e) => e.icon !== undefined || e.iconMarkup !== undefined)) {
      iconsAsked = true;
      inner.lazy(import('@molstar/seqx-viewer/header-icons'), (m) => {
        iconOf = m.headerIcon;
        inner.relayout();
      });
    }
    let n = Math.max(0, actionColumns(viewer.context.renderConfig).n);
    let slots: HeaderSlotPart[] | undefined;
    let states: HeaderStatePart[] | undefined;
    let menu = false;
    for (const e of all) {
      if (e.kind === 'state') {
        const text = e.text?.(f);
        if (text === undefined) continue;
        const count = Number.parseInt(text, 10);
        (states ??= []).push({ name: e.name, text, ...(Number.isFinite(count) ? { count } : {}) });
      } else if (n > 0 && slotted(e)) {
        n--;
        const command = commandOf(e.command);
        const label = typeof e.label === 'function' ? e.label(f) : (e.label ?? command?.title ?? e.name);
        // an action with a submenu opens it: enabled when any of its items runs
        const enabled = !!subsOf(e, f)?.length || whyNot(e.command, argsOf(e, f)) === undefined;
        const icon = iconOf?.(e.icon ?? command?.icon, e.iconMarkup);
        (slots ??= []).push({ name: e.name, label, ...(icon ? { icon } : {}), enabled });
      }
      if ((e.kind ?? 'action') === 'action' || e.command !== undefined) menu = true;
    }
    if (!slots && !states && !menu) return undefined;
    return {
      ...(slots ? { slots } : {}),
      ...(states ? { states } : {}),
      ...(menu ? { menu } : {}),
      ...(inner.hosted() ? { hosted: true } : {}),
    };
  };

  const keyOf = (command: string | undefined): string | undefined => {
    if (command === undefined) return undefined;
    const listed = viewer.config.view.keys[command];
    const chord = listed?.length ? parseChord(listed[0]) : commandOf(command)?.shortcut;
    return chord ? formatChord(chord, mac) : undefined;
  };

  /** An entry's menu item for a header's facts. */
  const itemOf = (e: HeaderEntry, f: HeaderFacts, slot = false): HeaderMenuItem => {
    const command = commandOf(e.command);
    const args = argsOf(e, f);
    const subs = subsOf(e, f);
    // an action with a submenu opens it: enabled when any of its items runs
    const why = subs?.length && (e.kind ?? 'action') === 'action' ? undefined : whyNot(e.command, args);
    const label = typeof e.label === 'function' ? e.label(f) : (e.label ?? command?.title ?? e.name);
    const text = e.text?.(f);
    const icon = e.icon ?? command?.icon;
    const key = keyOf(e.command);
    return {
      name: e.name,
      kind: e.kind ?? 'action',
      label,
      ...(icon ? { icon } : {}),
      ...(text !== undefined ? { text } : {}),
      ...(key ? { key } : {}),
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
    const all = applying(f, on);
    const n = Math.max(0, actionColumns(viewer.context.renderConfig).n);
    const slots = all.filter(slotted).slice(0, n);
    return [
      ...slots.map((e) => itemOf(e, f, true)),
      ...all.filter((e) => !slots.includes(e) && (e.kind ?? 'action') === 'action').map((e) => itemOf(e, f)),
      ...all.filter((e) => e.kind === 'state').map((e) => itemOf(e, f)),
    ];
  };

  const menuOf = (target: HeaderTarget, only?: string): HeaderMenu | undefined => {
    if (only !== undefined) {
      const f = factsOf(target);
      const e = inner.entries().find((x) => x.name === only);
      if (!f || !e) return undefined;
      const it = itemOf(e, f);
      // a chip's list: each thing it hides, then the entry itself ("Show all"); a slot's submenu: its items
      const own =
        it.kind === 'state' ? [...(it.items ?? []), { ...it, items: undefined } as HeaderMenuItem] : (it.items ?? [it]);
      return { target, groups: [own.map(({ items: _, ...x }) => x)] };
    }
    const groups: HeaderMenuItem[][] = [];
    const f = target.element ? factsOf(target) : undefined;
    if (f) groups.push(itemsOf(f, 'element'));
    const tf = target.track !== undefined ? factsOf({ section: target.section, track: target.track }) : undefined;
    if (tf) groups.push(itemsOf(tf, 'track'));
    const sf = factsOf({ section: target.section });
    if (sf) groups.push(itemsOf(sf, 'section'));
    const out = groups.filter((g) => g.length);
    return out.length ? { target, groups: out } : undefined;
  };

  const openMenu = (target: HeaderTarget, place: HeaderMenuPlace, only?: string): boolean => {
    if (inner.disposed) return false;
    const menu = menuOf(target, only);
    if (!menu) return false;
    // the card goes: the menu is what the pointer reads now
    (viewer.behaviors.get('hoverCard') as HoverCardBehavior | undefined)?.interrupt();
    inner.show(menu, place);
    return true;
  };
  const run = (name: string, target: HeaderTarget): string | Promise<string> | undefined => {
    const e = inner.entries().find((x) => x.name === name);
    const f = factsOf(target);
    if (!e || !f || e.command === undefined) return undefined;
    return viewer.execute(e.command, argsOf(e, f));
  };
  const runItem = (item: HeaderMenuItem): string | Promise<string> | undefined => {
    if (item.command === undefined) return undefined;
    // a preview ends in the choice's patch (the command cancels it first), not in nothing
    inner.previewed = false;
    return viewer.execute(item.command, item.args ?? {});
  };
  // a choice hovered: its command's preview (`layer.channel`'s through the theming behavior), when it has one
  const preview = (item: HeaderMenuItem): boolean => {
    const c = commandOf(item.command);
    if (!c?.preview || !item.args || whyNot(item.command, item.args) !== undefined) return false;
    inner.previewed = c.preview({}, item.args) === true;
    return inner.previewed;
  };

  /** A header line's slot boxes as laid out (the canvas's, or those a host's DOM parts are placed over). */
  const slotBoxes = (line: HeaderLine) => {
    const box = (b: Box) => ({ x: b.x, y: b.y, width: b.width, height: b.height });
    const out = line.slots.map((s) => ({ name: s.part.name, box: box(s) }));
    if (line.menu) out.push({ name: HeaderMenuName, box: box(line.menu) });
    return out;
  };

  const geometryOf = (layout: ViewLayout): HeaderGeometry[] => {
    const out: HeaderGeometry[] = [];
    const h = viewer.height || Number.POSITIVE_INFINITY;
    for (const s of layout.sections) {
      const lines: [HeaderLine, HeaderTarget, Box | undefined, number | undefined][] = [];
      if (s.header) lines.push([s.header, { section: s.section.address }, undefined, undefined]);
      const seen = new Set<string>();
      for (const t of s.tracks) {
        if (!t.header || seen.has(t.track.address)) continue;
        seen.add(t.track.address);
        lines.push([
          t.header,
          { section: s.section.address, track: t.track.address },
          t.well,
          s.config.inset + t.labelWidth,
        ]);
      }
      for (const [line, target, well, labelEnd] of lines) {
        const b = line.line;
        out.push({
          target,
          line: b,
          label: line.label,
          slots: slotBoxes(line),
          states: line.states.map((st) => ({ name: st.name, box: st.box, text: st.text, mark: st.mark })),
          ...(well ? { well } : {}),
          ...(labelEnd !== undefined ? { labelEnd } : {}),
          ...(line.ellipsis !== undefined ? { ellipsis: line.ellipsis } : {}),
          shown: b.y + b.height > 0 && b.y < h && b.y >= s.y - 1 && b.y + b.height <= s.y + s.height + 1,
        });
      }
    }
    return out;
  };

  /** Whether a header's actions show now (`header.actions.show`: the hover, the cursor). */
  const showing = (g: HeaderGeometry): boolean => {
    const s = sectionLayout(viewer.frame?.layout, g.target.section);
    const address = g.target.track ?? g.target.section;
    const config = viewer.context.renderConfig;
    return !!s && actionsShown(address, s.id, viewer.emphasis, config);
  };

  const partsState = (): HeaderPartsState => {
    const layout = viewer.frame?.layout;
    const headers: HeaderPartsHeader[] = [];
    for (const g of layout ? geometryOf(layout) : []) {
      const f = g.shown ? factsOf(g.target) : undefined;
      if (f) headers.push({ geometry: g, actions: showing(g), items: itemsOf(f, on(g.target)) });
    }
    return { headers };
  };

  /** What a secondary click at a point opens: an element's, a track's or a section's menu. */
  const targetAt = (x: number, y: number, hit: RenderHit | undefined): HeaderTarget | undefined => {
    const layout = viewer.frame?.layout;
    if (!layout) return undefined;
    const header = hit && headerTarget(hit, layout);
    if (header) return header;
    if (hit)
      for (const s of layout.sections) {
        // a strip's layer, else a section's (an extension's ribbons), on the track that holds it
        const sl = s.sectionLays.find((l) => l.key === hit.layer);
        for (const t of s.tracks) {
          const l =
            sl && sl.layer.address?.startsWith(`${t.track.address}/`) ? sl : t.layers.find((x) => x.key === hit.layer);
          if (l)
            return {
              section: s.section.address,
              track: t.track.address,
              // the element with its layer's type: an extension's entries are scoped by it (its stem menu)
              ...(hit.element ? { element: { layer: hit.layer, key: hit.element.key, type: l.type.name } } : {}),
            };
        }
      }
    // a track's band under the point, else its section; on a row of blocks the block under it (its label's menu)
    const section = viewer.input.sectionAt(y);
    const s = section === undefined ? undefined : sectionLayout(layout, section);
    if (!s) return undefined;
    const p = s.placement;
    const b = hit?.block ?? (p instanceof BlockPlacement && !p.stacked && s.section.blocks ? p.row.blockAt(x) : -1);
    if (b >= 0 && s.section.blocks) return { section: s.section.address, block: b };
    const t = s.tracks.find((t) => y >= t.box.y && y < t.box.y + t.box.height);
    return t ? { section: s.section.address, track: t.track.address } : { section: s.section.address };
  };

  const words = (place: Place) => {
    if (place.layer === undefined) return undefined;
    const target = headerTarget({ layer: place.layer } as RenderHit, viewer.frame?.layout);
    const f = target && factsOf(target);
    return f ? headerWords(place, f, itemsOf(f, on(target!))) : undefined;
  };

  const click = (hit: RenderHit, x: number, y: number, touch: boolean): void => {
    const t = hit.target;
    const target = headerTarget(hit, viewer.frame?.layout);
    if (!target) return;
    if (t === 'header-label' || t === 'section-header') {
      // no hover on touch: a tap on the label shows its slots (as the pointer resting there would), a tap on a slot
      // then runs it
      if (touch) viewer.input.move(x, y, { pointerType: 'touch' });
      return;
    }
    const name = hit.state ?? '';
    const g = inner.geometry().find((g) => g.target.track === target.track && g.target.section === target.section);
    const box = g?.slots.find((s) => s.name === name)?.box ?? g?.states.find((s) => s.name === name)?.box;
    const at = box ? { x: box.x, y: box.y + box.height } : { x, y };
    if (name === HeaderMenuName) return void openMenu(target, at);
    const e = inner.entries().find((e) => e.name === name);
    if (!e) return;
    // an action with a submenu opens it at its slot; a chip lists what it hides under `hiddenChip: menu`
    const f = factsOf(target);
    const subs = f && subsOf(e, f);
    if (subs?.length && (e.kind !== 'state' || inner.params.hiddenChip === 'menu'))
      return void openMenu(target, at, e.name);
    void run(e.name, target);
  };

  const cursorTarget = (): HeaderTarget | undefined => {
    const c = viewer.selection.cursor.value;
    if (!c) return undefined;
    if (c.track !== undefined)
      return c.track === c.section ? { section: c.section } : { section: c.section, track: c.track };
    const s = viewer.view.resolved.value?.sections.find((x) => x.address === c.section);
    const first = s?.tracks[0];
    return s ? { section: s.address, ...(first ? { track: first.address } : {}) } : undefined;
  };
  /** When the keyboard (or a long press) last opened the menu: the secondary click after it is the same ask. */
  let askedAt = Number.NEGATIVE_INFINITY;
  /** Opens the menu of the cursor's header under its line's label start. */
  const menuAtCursor = (): CommandResult => {
    const target = cursorTarget();
    const g =
      target && inner.geometry().find((x) => x.target.section === target.section && x.target.track === target.track);
    const ok = !!target && openMenu(target, g ? { x: g.label.x, y: g.line.y + g.line.height } : { x: 0, y: 0 });
    if (ok) askedAt = DomHost.now();
    return { message: ok ? 'Actions' : '' };
  };

  const contextClick = inner.params.contextClick ?? 'everywhere';
  const context = (x: number, y: number, e: { readonly shiftKey?: boolean } = {}): boolean => {
    // Shift+secondary click is always the browser's; the menu key's own secondary click is the same ask (EC-156)
    if (e.shiftKey || contextClick === 'none') return false;
    if (DomHost.now() - askedAt < AskMs) {
      askedAt = Number.NEGATIVE_INFINITY;
      return true;
    }
    const hit = viewer.hit(x, y);
    const c = viewer.context.renderConfig;
    if (contextClick === 'header' && !(x >= c.inset && x < c.inset + c.labelWidth)) return false;
    const target = targetAt(x, y, hit);
    if (!target) return false;
    // a chip's own list (each thing it hides, then all)
    const only = hit?.target === 'header-state' && hit.state ? hit.state : undefined;
    const known = only !== undefined && inner.entries().find((x) => x.name === only)?.items;
    return openMenu(target, { x, y }, known ? only : undefined);
  };

  const pressed = (x: number, y: number) => {
    if (!inner.behavior.menu.value && context(x, y)) askedAt = DomHost.now();
  };
  const published = (layout: ViewLayout): readonly HeaderGeometry[] => {
    const g = geometryOf(layout);
    return sameGeometry(inner.geometry(), g) ? inner.geometry() : g;
  };

  // the cursor's moves onto and between headers (Enter on one opens its menu)
  const keepCursor = headerCursor(viewer, menuAtCursor);
  return {
    factsOf: (t) => factsOf(t),
    layoutParts,
    menuOf,
    openMenu,
    run,
    runItem,
    preview,
    geometryOf,
    partsState,
    targetAt,
    words,
    click,
    context,
    cursorTarget,
    menuAtCursor,
    pressed,
    published,
    keepCursor,
  };
}

/**
 * The words of a header's part (_spec/README.md "Behaviors", Headers: Words): a slot its entry's label and key, and why
 * not; a chip its text and what a click does; a label its track's label, its states, and its layers (the card) or the
 * key that opens its actions (the cursor, `state: 'cursor'`).
 */
function headerWords(
  place: { readonly target?: string | undefined; readonly state?: string | undefined },
  f: HeaderFacts,
  items: readonly HeaderMenuItem[],
): { title: string; text: string } | undefined {
  const t = place.target;
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
}

const sameBox = (a: Box | undefined, b: Box | undefined) =>
  a === b || (!!a && !!b && a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height);

/** Two geometries alike (the published one changes only when a box moved). */
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

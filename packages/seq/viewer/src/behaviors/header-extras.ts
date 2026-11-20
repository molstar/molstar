/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  actionColumns,
  type Box,
  type HeaderLine,
  HeaderMenuName,
  type Hit as RenderHit,
  type SectionLayout,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import { actionsShown } from '@molstar/seqx-core/render/headers-draw';
import type { CommandResult } from '@molstar/ui/commands/types';
import type { Place } from '@molstar/seqx-viewer/describe';
import {
  type HeaderEntry,
  type HeaderFacts,
  type HeaderGeometry,
  type HeaderMenu,
  type HeaderMenuItem,
  type HeaderPartsHeader,
  type HeaderPartsState,
  type HeaderTarget,
} from '@molstar/seqx-viewer/header-entries';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { headerCursor } from './header-cursor.js';

// The `headers` behavior's late parts (_spec/README.md "Behaviors", Headers; "Composition", Injected, opt in: a listed
// part used late loads on first use): a lazy chunk the behavior asks for once composed. The menu's data, the headers'
// geometry and a host's DOM parts' state, what a secondary click opens, the words of a header's parts and the cursor on
// headers. Until it lands the slots and chips still draw and run (the hover path is the behavior's own); then the
// geometry is published and the menu, the words and the cursor's header moves work.

/** What the late parts read of the behavior. */
export interface HeaderInner {
  readonly viewer: SeqViewer;
  /** A host's DOM parts are given (the canvas draws no slot). */
  readonly dom: boolean;
  entries(): readonly HeaderEntry[];
  factsOf(target: HeaderTarget): HeaderFacts | undefined;
  itemOf(e: HeaderEntry, f: HeaderFacts): HeaderMenuItem;
  itemsOf(f: HeaderFacts, on: HeaderEntry['on']): HeaderMenuItem[];
  /** Opens the menu at the keyboard cursor (Enter on a header). */
  menuAtCursor(): CommandResult;
  openMenu(target: HeaderTarget, place: { x: number; y: number }, only?: string): boolean;
  run(name: string, target: HeaderTarget): unknown;
  geometry(): readonly HeaderGeometry[];
  /** A count chip's click lists what it hides (`hiddenChip: menu`). */
  readonly hiddenChipMenu: boolean;
}

export interface HeaderExtras {
  menuOf(target: HeaderTarget, only?: string): HeaderMenu | undefined;
  geometryOf(layout: ViewLayout): HeaderGeometry[];
  partsState(): HeaderPartsState;
  targetAt(x: number, y: number, hit: RenderHit | undefined): HeaderTarget | undefined;
  words(place: Place): { title: string; text: string } | undefined;
  /** A click on a header's part (a slot runs its entry, a chip its state, a label on touch shows the slots). */
  click(hit: RenderHit, x: number, y: number, touch: boolean): void;
  /** What the menu key opens: the cursor's header, else its section's first track's. */
  cursorTarget(): HeaderTarget | undefined;
  /** Opens the menu at the keyboard cursor: whether it opened. */
  menuAtCursor(): boolean;
  /** A secondary click at a point (`header`: in the label column only): true when the menu took it. */
  context(x: number, y: number, header: boolean): boolean;
  /** The frame's: the cursor off a track a layout hid. */
  keepCursor(): void;
}

const sectionLayout = (layout: ViewLayout | undefined, address: string): SectionLayout | undefined =>
  layout?.sections.find((s) => s.section.address === address);

/** The header a hit is on (a track's line, a slot, a chip, a section's line): its target. */
export function headerTarget(hit: RenderHit, layout: ViewLayout | undefined): HeaderTarget | undefined {
  for (const s of layout?.sections ?? []) {
    if (s.section.address === hit.layer) return { section: s.section.address };
    for (const t of s.tracks)
      if (t.track.address === hit.layer) return { section: s.section.address, track: t.track.address };
  }
  return undefined;
}

export function headerExtras(inner: HeaderInner): HeaderExtras {
  const { viewer, dom } = inner;
  const on = (t: HeaderTarget) => (t.track ? 'track' : 'section');

  const menuOf = (target: HeaderTarget, only?: string): HeaderMenu | undefined => {
    if (only !== undefined) {
      const f = inner.factsOf(target);
      const e = inner.entries().find((x) => x.name === only);
      if (!f || !e) return undefined;
      const it = inner.itemOf(e, f);
      // a chip's list: each thing it hides, then the entry itself ("Show all"); a slot's submenu: its items
      const own =
        it.kind === 'state' ? [...(it.items ?? []), { ...it, items: undefined } as HeaderMenuItem] : (it.items ?? [it]);
      return { target, groups: [own.map(({ items: _, ...x }) => x)] };
    }
    const groups: HeaderMenuItem[][] = [];
    const f = target.element ? inner.factsOf(target) : undefined;
    if (f) groups.push(inner.itemsOf(f, 'element'));
    const tf = target.track !== undefined ? inner.factsOf({ section: target.section, track: target.track }) : undefined;
    if (tf) groups.push(inner.itemsOf(tf, 'track'));
    const sf = inner.factsOf({ section: target.section });
    if (sf) groups.push(inner.itemsOf(sf, 'section'));
    const out = groups.filter((g) => g.length);
    return out.length ? { target, groups: out } : undefined;
  };

  /** A header line's slot boxes: as laid out (canvas), or the action columns' a host's DOM parts fill. */
  const slotBoxes = (line: HeaderLine, items: readonly HeaderMenuItem[]) => {
    const box = (b: Box) => ({ x: b.x, y: b.y, width: b.width, height: b.height });
    if (!dom) {
      const out = line.slots.map((s) => ({ name: s.part.name, box: box(s) }));
      if (line.menu) out.push({ name: HeaderMenuName, box: box(line.menu) });
      return out;
    }
    const cols = actionColumns(viewer.context.renderConfig);
    if (cols.n < 0 || line.short) return [];
    const px = viewer.context.renderConfig.header.slotPx;
    const y = line.label.y + line.label.height / 2 - px / 2;
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
    const h = viewer.height || Number.POSITIVE_INFINITY;
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
        const f = dom ? inner.factsOf(target) : undefined;
        const b = line.line;
        out.push({
          target,
          line: b,
          label: line.label,
          slots: slotBoxes(line, f ? inner.itemsOf(f, on(target)) : []),
          states: line.states.map((st) => ({ name: st.name, box: st.box, text: st.text, mark: st.mark })),
          ...(well ? { well } : {}),
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
    return !!s && actionsShown(address, s.id, viewer.emphasis, config, g.target.track !== undefined);
  };

  const partsState = (): HeaderPartsState => {
    const layout = viewer.frame?.layout;
    const headers: HeaderPartsHeader[] = [];
    for (const g of layout ? geometryOf(layout) : []) {
      const f = g.shown ? inner.factsOf(g.target) : undefined;
      if (f) headers.push({ geometry: g, actions: showing(g), items: inner.itemsOf(f, on(g.target)) });
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
      for (const s of layout.sections)
        for (const t of s.tracks)
          if (t.layers.some((l) => l.key === hit.layer))
            return {
              section: s.section.address,
              track: t.track.address,
              ...(hit.element ? { element: { layer: hit.layer, key: hit.element.key } } : {}),
            };
    // a track's band under the point, else its section
    const section = viewer.input.sectionAt(y);
    const s = section === undefined ? undefined : sectionLayout(layout, section);
    if (!s) return undefined;
    const t = s.tracks.find((t) => y >= t.box.y && y < t.box.y + t.box.height);
    return t ? { section: s.section.address, track: t.track.address } : { section: s.section.address };
  };

  const words = (place: Place) => {
    if (place.layer === undefined) return undefined;
    const target = headerTarget({ layer: place.layer } as RenderHit, viewer.frame?.layout);
    const f = target && inner.factsOf(target);
    return f ? headerWords(place, f, inner.itemsOf(f, on(target!))) : undefined;
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
    if (name === HeaderMenuName) return void inner.openMenu(target, at);
    const e = inner.entries().find((e) => e.name === name);
    if (!e) return;
    // an action with a submenu opens it at its slot; a chip lists what it hides under `hiddenChip: menu`
    const f = inner.factsOf(target);
    const subs = f && e.items?.(f);
    if (subs?.length && (e.kind !== 'state' || inner.hiddenChipMenu)) return void inner.openMenu(target, at, e.name);
    void inner.run(e.name, target);
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
  /** Opens the menu of the cursor's header under its line's label start. */
  const menuAtCursor = (): boolean => {
    const target = cursorTarget();
    if (!target) return false;
    const g = inner.geometry().find((x) => x.target.section === target.section && x.target.track === target.track);
    return inner.openMenu(target, g ? { x: g.label.x, y: g.line.y + g.line.height } : { x: 0, y: 0 });
  };

  const context = (x: number, y: number, header: boolean): boolean => {
    const hit = viewer.hit(x, y);
    const c = viewer.context.renderConfig;
    if (header && !(x >= c.inset && x < c.inset + c.labelWidth)) return false;
    const target = targetAt(x, y, hit);
    if (!target) return false;
    // a chip's own list (each thing it hides, then all)
    const only = hit?.target === 'header-state' && hit.state ? hit.state : undefined;
    const known = only !== undefined && inner.entries().find((x) => x.name === only)?.items;
    return inner.openMenu(target, { x, y }, known ? only : undefined);
  };

  // the cursor's moves onto and between headers (Enter on one opens its menu)
  const keepCursor = headerCursor(viewer, inner.menuAtCursor);
  return { menuOf, geometryOf, partsState, targetAt, words, click, context, cursorTarget, menuAtCursor, keepCursor };
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

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ElementKey } from '@molstar/seqx-core/query';
import { DefaultNumbering, findNumbering, numberingIds } from '@molstar/seqx-core/query';
import {
  itemText,
  type ItemDescription,
  type LayerLayout,
  LettersType,
  type RowsLay,
  type SectionLayout,
  shownLetters,
  spanText,
  type ValuesLay,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { AddressModel, SpecModel, StoreModel } from '@molstar/seqx-viewer/models';
import type { HoverCard, HoverCardRow } from './hover-card.js';

// What a place is, in words (_spec/README.md "Input", Announcements; "Hover card"): one description that the live
// region says and the hover card shows. A place (the hover's record, the keyboard cursor's cell) goes through the
// sections in order; each adds its part: a row of the card (a label and a value) and its words in the sentence. The
// card's `text` is the sentence the announcer says, so a reader of the canvas and a reader of the live region are told
// the same things. The sections are internal (extensibility.md: chrome parts behind internal interfaces until an outside
// user needs one): not exported.

/** What the words read: the models a behavior is given, and the last frame's layout. */
export interface DescribeSource {
  readonly spec: SpecModel;
  readonly store: StoreModel;
  readonly addresses: AddressModel;
  layout(): ViewLayout | undefined;
}

/** A place to describe: the hover's record (`HoverRecord`) or a cell (the keyboard cursor's). */
export interface Place {
  readonly section: string | undefined;
  readonly layer?: string | undefined;
  readonly system: string | undefined;
  /** NaN when none. */
  readonly position: number;
  readonly element?: ElementKey | undefined;
  readonly ySystem?: string | undefined;
  /** NaN when none. */
  readonly yPosition?: number;
  readonly target?: string | undefined;
  /** The block (its index in the section) and a grid cell's y block, when the section has blocks. */
  readonly block?: number | undefined;
  readonly yBlock?: number | undefined;
  /** What a layer named within the position (a rows column cell: `title: value`). */
  readonly state?: string | undefined;
}

/** The parts a description collects. */
class Parts {
  title = '';
  readonly rows: HoverCardRow[] = [];
  readonly words: string[] = [];
  /** Set by a section that says all there is (an item: its own words): the place's sections then say nothing. */
  whole = false;

  row(label: string, value: string): void {
    this.rows.push({ label, value });
  }
}

/** What a section reads. */
interface PlaceContext {
  readonly words: Describer;
  readonly model: StaticModel;
  readonly section: SectionLayout;
  readonly layer: LayerLayout | undefined;
  readonly place: Place;
  /** The position, floored; NaN when none. */
  readonly position: number;
  readonly yPosition: number;
}

/** A part of a description (internal: extensibility.md, chrome parts behind internal interfaces). */
interface CardSection {
  readonly name: string;
  describe(at: PlaceContext, out: Parts): void;
}

/** A value as a card says it: up to three significant decimals, no trailing zeros. */
export function valueText(v: number): string {
  if (!Number.isFinite(v)) return String(v);
  if (Number.isInteger(v)) return v.toLocaleString('en-US');
  const a = Math.abs(v);
  const digits = a >= 100 ? 1 : a >= 1 ? 2 : 3;
  return Number(v.toFixed(digits)).toString();
}

// ---- the sections, in order ----

/** A block of a section of blocks (a chain of a complex): its label, and which copy when its system has several. */
const BlockSection: CardSection = {
  name: 'block',
  describe(at, out) {
    const blocks = at.section.section.blocks;
    const k = at.place.block;
    if (!blocks || k === undefined || !blocks[k]) return;
    const say = (b: number, which: string) => {
      const block = blocks[b];
      if (!block) return;
      out.row(which, block.label);
      let copy = '';
      if (block.members.length > 1) copy = `${block.members.length} copies merged`;
      else {
        const twins = blocks.filter((o) => o.axis.system === block.axis.system);
        if (twins.length > 1) copy = `copy ${twins.indexOf(block) + 1} of ${twins.length}`;
      }
      if (copy) out.row('Copy', copy);
      out.words.push(`${which.toLowerCase()} ${block.label}${copy ? `, ${copy}` : ''}`);
    };
    const y = at.place.yBlock;
    if (y !== undefined && y !== k) {
      say(k, 'Chain x');
      say(y, 'Chain y');
    } else say(k, 'Chain');
  },
};

/** A pairs layer's pair: its label, its two ends, its kind. */
const PairSection: CardSection = {
  name: 'pair',
  describe(at, out) {
    const el = at.place.element;
    const l = at.layer;
    if (!el || !l || l.type.name !== 'pairs') return;
    const owner = l.layer.data?.address;
    const a = owner ? at.model.annotations.get(owner) : undefined;
    const d = a?.pairs;
    if (!a || !d) return;
    const k = d.keys.indexOf(el.key);
    if (k < 0) return;
    const label = d.labels?.[k] || String(el.key);
    const x = `${at.words.systemName(a.target)} ${spanText(d.aStart[k], d.aEnd[k])}`;
    const y = `${at.words.systemName(a.ySystem ?? a.target)} ${spanText(d.bStart[k], d.bEnd[k])}`;
    const kind = d.kinds[d.kind[k]];
    out.title = label;
    out.row('From', x);
    out.row('To', y);
    if (kind) out.row('Kind', kind);
    out.words.push(`${label}, ${x} to ${y}${kind ? `, ${kind}` : ''}`);
    out.whole = true;
  },
};

/** An item (a range, a feature): its label, its extent, its category. */
const ItemSection: CardSection = {
  name: 'item',
  describe(at, out) {
    const el = at.place.element;
    const l = at.layer;
    if (out.whole || !el || !l?.type.items) return;
    const it = at.words.itemsOf(l).find((d) => d.key === el.key);
    if (!it) return;
    const system = at.words.systemName(it.system ?? l.layer.data?.system ?? at.section.section.axis.system);
    const category = at.words.categoryOf(l, el.key);
    out.title = it.label;
    out.row('Extent', `${system} ${spanText(it.start, it.end)}`);
    if (category) out.row('Category', category);
    out.words.push(`${itemText(system, it)}${category ? `, ${category}` : ''}`);
    out.whole = true;
  },
};

/** A rows layer's row (an alignment's sequence): its label, and the column cell a row's columns name. */
const RowSection: CardSection = {
  name: 'row',
  describe(at, out) {
    const l = at.layer;
    // by its type's name, not its class: the rows layer stays out of a bundle without it
    if (out.whole || !l || l.type.name !== 'rows') return;
    const lay = l.lay as RowsLay;
    const e = at.place.target === 'row' ? at.position : at.place.ySystem === lay.rowsSystem ? at.yPosition : Number.NaN;
    if (Number.isNaN(e)) return;
    const label = lay.label(e);
    out.title = label;
    out.row('Row', label);
    out.words.push(`row ${label}`);
    const state = at.place.state;
    const colon = state ? state.indexOf(': ') : -1;
    if (state && colon > 0) {
      out.row(state.slice(0, colon), state.slice(colon + 2));
      out.words.push(state);
    }
    // a row's header is the row alone: its place on the rows system is no position a reader counts
    if (at.place.target === 'row') out.whole = true;
  },
};

/** A grid's cell: its two positions and its value. */
const GridSection: CardSection = {
  name: 'grid',
  describe(at, out) {
    const p = at.place;
    const l = at.layer;
    if (out.whole || !l || l.type.name !== 'grid' || p.system === undefined || p.ySystem === undefined) return;
    if (Number.isNaN(at.position) || Number.isNaN(at.yPosition)) return;
    // over blocks (a grid over a complex) each axis by its block's label, the chain: `CDK2 33 × cyclin A2 263`
    const blocks = at.section.section.blocks;
    const nameOf = (system: string, k: number | undefined) =>
      (k !== undefined && blocks?.[k]?.axis.system === system ? blocks[k].label : undefined) ??
      at.words.systemName(system);
    // a residue as a structure's readout writes it: its letter before its number (`CDK2 K33`), when it is loaded
    const end = (system: string, k: number | undefined, i: number) =>
      `${nameOf(system, k)} ${at.words.letterOn(system, i)}${spanText(i, i + 1)}`;
    const x = end(p.system, p.block, at.position);
    const y = end(p.ySystem, p.yBlock ?? p.block, at.yPosition);
    const v = at.words.gridValue(l, p.system, at.position, p.ySystem, at.yPosition);
    out.title = `${x} × ${y}`;
    out.row('x', x);
    out.row('y', y);
    if (!Number.isNaN(v)) out.row('Value', valueText(v));
    out.words.push(`${at.section.id}, ${x} × ${y}${Number.isNaN(v) ? '' : `, value ${valueText(v)}`}`);
    out.whole = true;
  },
};

/**
 * A cell: the section, the system and the position in every numbering of the system (the ordinal, the index, the
 * declared ones: the section's shown numbering first), the letter shown there, a value track's value and the items
 * under it.
 */
const PlaceSection: CardSection = {
  name: 'place',
  describe(at, out) {
    const p = at.place;
    if (out.whole || p.system === undefined || Number.isNaN(at.position)) return;
    const w = at.words;
    const system = p.system;
    const position = at.position;
    const name = w.systemName(system);
    const declared = w.numberedLabels(at.section, system, position);
    const ordinal = spanText(position, position + 1);
    // the title says the system and the ordinal when nothing else named the place; else a row does
    if (!out.title) out.title = `${name} ${ordinal}`;
    else out.row(name, ordinal);
    out.row('index', String(position));
    for (const [id, label] of declared) out.row(id, label);
    const labels = declared.map(([id, label]) => `${id} ${label}`).join(', ');
    out.words.push(`${at.section.id}, ${name} ${ordinal}${labels ? ` (${labels})` : ''}`);
    const letter = w.letterAt(at.section, position);
    if (letter) {
      out.row('Letter', letter);
      out.words.push(`letter ${letter}`);
    }
    for (const [label, value] of w.valuesAt(at.section, position)) {
      out.row(label, value);
      out.words.push(`${label} ${value}`);
    }
    const items = w.itemsAt(at.section, position);
    for (const it of items) out.row(it.label, it.extent);
    if (items.length) out.words.push(`items: ${items.map((i) => `${i.label}, ${i.extent}`).join('; ')}`);
  },
};

// the block last: what is in it is said first, then its chain
const Sections: readonly CardSection[] = [
  PairSection,
  ItemSection,
  RowSection,
  GridSection,
  PlaceSection,
  BlockSection,
];

/** The words of places in a viewer: the hover's, the cursor's. Items are read once a frame. */
export class Describer {
  private readonly _src: DescribeSource;
  /** Items by layer layout, made once a frame (a description looks them up, a frame changes them). */
  private _items = new Map<LayerLayout, readonly ItemDescription[]>();
  private _itemsLayout: unknown;

  constructor(src: DescribeSource) {
    this._src = src;
  }

  /** The description of a place in the last frame; undefined when it is in none of its sections. */
  describe(place: Place): HoverCard | undefined {
    const layout = this._src.layout();
    if (!layout || place.section === undefined) return undefined;
    const section = layout.sections.find((s) => s.section.address === place.section);
    if (!section) return undefined;
    const layer = place.layer === undefined ? undefined : layerOf(section, place.layer);
    const out = new Parts();
    const at: PlaceContext = {
      words: this,
      model: this._src.spec.model.value,
      section,
      layer,
      place,
      position: Number.isNaN(place.position) ? Number.NaN : Math.floor(place.position),
      yPosition:
        place.yPosition === undefined || Number.isNaN(place.yPosition) ? Number.NaN : Math.floor(place.yPosition),
    };
    for (const s of Sections) s.describe(at, out);
    if (!out.words.length) return undefined;
    return { title: out.title || out.words[0], rows: out.rows, text: out.words.join(', ') };
  }

  /** A cell of a section: the keyboard cursor's. */
  cell(section: string, system: string, position: number): HoverCard | undefined {
    return this.describe({ section, system, position });
  }

  /** The text of a place ('' when there is none). */
  text(place: Place): string {
    return this.describe(place)?.text ?? '';
  }

  // ---- what the sections read ----

  systemName(system: string): string {
    return this._src.spec.model.value.systems.get(system)?.name ?? system;
  }

  itemsOf(l: LayerLayout): readonly ItemDescription[] {
    const layout = this._src.layout();
    if (layout !== this._itemsLayout) {
      this._itemsLayout = layout;
      this._items = new Map();
    }
    let items = this._items.get(l);
    if (!items) this._items.set(l, (items = l.type.items?.(l.lay) ?? []));
    return items;
  }

  /** An item's category: the `category` field of an inline annotation's item (a source's items are not held here). */
  categoryOf(l: LayerLayout, key: string | number): string {
    const owner = l.layer.data?.address;
    const a = owner ? this._src.spec.model.value.annotations.get(owner) : undefined;
    const f = a?.items?.find((i) => i.key === key)?.fields?.category;
    return typeof f === 'string' ? f : '';
  }

  /** The items of a section at a position: its tracks' layers' items under the position's x. */
  itemsAt(section: SectionLayout, position: number, only?: LayerLayout): { label: string; extent: string }[] {
    const axis = section.axis;
    if (!axis) return [];
    const a = axis.x(position);
    const b = axis.x(position + 1);
    if (!Number.isFinite(a) || !Number.isFinite(b)) return [];
    const x = (a + b) / 2;
    const out: { label: string; extent: string }[] = [];
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (only && l !== only) continue;
        const data = l.layer.data?.system;
        for (const it of this.itemsOf(l))
          if (it.y === undefined && it.x0 <= x && x <= it.x1)
            out.push({
              label: it.label,
              extent: `${this.systemName(it.system ?? data ?? section.section.axis.system)} ${spanText(it.start, it.end)}`,
            });
      }
    return out;
  }

  /** The letter a letters layer of the section shows at a position (none while it is not loaded). */
  letterAt(section: SectionLayout, position: number): string {
    const model = this._src.spec.model.value;
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (l.type.name !== 'letters') continue;
        const r = shownLetters(LettersType.prepare(l.layer, model), section.section.axis.dir, this._src.store.store, [
          position,
        ]);
        if (r.text.trim()) return r.text.trim();
      }
    return '';
  }

  /** A system's letter at a position: inline, else in a loaded tile of its source ('' while none holds it). */
  letterOn(system: string, position: number): string {
    const s = this._src.spec.model.value.systems.get(system);
    if (s?.letters !== undefined) return s.letters[position] ?? '';
    if (!s?.source) return '';
    const store = this._src.store.store;
    for (const k of store.keysOf({ source: s.source, system: s.name, start: position, end: position + 1, level: 1 })) {
      const e = store.get(k);
      const t = e?.status === 'loaded' ? e.tables.letters : undefined;
      if (!t || position < t.start || position >= t.start + t.letters.length) continue;
      const i = position - t.start;
      return typeof t.letters === 'string' ? t.letters[i] : (t.alphabet?.[t.letters[i]] ?? '');
    }
    return '';
  }

  /** Each value track's value at a position, by its track's label (a bin's low and high when they differ). */
  valuesAt(section: SectionLayout, position: number): [string, string][] {
    const axis = section.axis;
    if (!axis) return [];
    const x = (axis.x(position) + axis.x(position + 1)) / 2;
    if (!Number.isFinite(x)) return [];
    const out: [string, string][] = [];
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (l.type.name !== 'values') continue;
        const lay = l.lay as ValuesLay;
        for (let j = 0; j < lay.n; j++) {
          if (!(lay.x0[j] <= x && x <= lay.x1[j])) continue;
          const lo = lay.lo[j];
          const hi = lay.v[j];
          if (Number.isNaN(hi)) break;
          const label = t.track.label ?? l.key;
          out.push([label, lo === hi || Number.isNaN(lo) ? valueText(hi) : `${valueText(lo)}–${valueText(hi)}`]);
          break;
        }
      }
    return out;
  }

  /** A grid's value at a cell (the hovered layer's grid, on its two systems), NaN when none. */
  gridValue(l: LayerLayout, system: string, i: number, ySystem: string, j: number): number {
    const owner = l.layer.data?.address;
    const a = owner ? this._src.spec.model.value.annotations.get(owner) : undefined;
    const g = a?.grid;
    if (!a || !g || a.target !== system || (a.ySystem ?? a.target) !== ySystem) return Number.NaN;
    if (i < 0 || j < 0 || i >= g.nx || j >= g.ny) return Number.NaN;
    return g.values[i * g.ny + j];
  }

  /**
   * The position's labels in the system's declared numberings, with their ids: the numbering the section's layers show
   * first (as a ruler or a text's numbers do), then the others; a numbering with no label there is left out.
   */
  numberedLabels(section: SectionLayout, system: string, position: number): [string, string][] {
    const s = this._src.addresses.queryModel().system(system);
    if (!s) return [];
    const ids: string[] = [];
    for (const t of section.section.tracks)
      for (const l of t.layers) {
        const id = l.params.numbering;
        if (typeof id === 'string' && id && id !== DefaultNumbering && id !== 'index' && !ids.includes(id))
          ids.push(id);
      }
    for (const id of numberingIds(s)) if (id !== DefaultNumbering && id !== 'index' && !ids.includes(id)) ids.push(id);
    const out: [string, string][] = [];
    for (const id of ids) {
      const n = findNumbering(s, id);
      const label = n?.numbering.labelAt(position - n.shift);
      if (label !== undefined) out.push([id, label]);
    }
    return out;
  }
}

/** A layer of a section by its key: its tracks', its strips' (a section of blocks), its section layers'. */
function layerOf(section: SectionLayout, key: string): LayerLayout | undefined {
  for (const t of section.tracks) for (const l of t.layers) if (l.key === key) return l;
  for (const ls of section.lays) for (const l of ls) if (l.key === key) return l;
  for (const l of section.sectionLays) if (l.key === key) return l;
  return undefined;
}

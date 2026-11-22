/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Key, Location } from '@molstar/seqx-core/query';
import { DefaultNumbering, findNumbering, numberingIds, rangesLocation } from '@molstar/seqx-core/query';
import { extentText } from '@molstar/seqx-core/render/accessibility';
import {
  type GridLay,
  itemText,
  layerSource,
  type ItemDescription,
  type LayerLayout,
  LettersType,
  type RowsLay,
  type SectionLayout,
  shownLetters,
  spanText,
  gridSummary,
  type ValuesLay,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { AddressModel, Cursor, SpecModel, StoreModel, ViewModel } from '@molstar/seqx-viewer/models';
import { chemistryRows } from './describe-molecules.js';
import type { HoverCard, HoverCardRow } from './hover-card.js';

// What a place is, in words (_spec/README.md "Input", Announcements; "Hover card"): one description that the live
// region says and the hover card shows. A place (the hover's record, the keyboard cursor's cell) goes through the
// sections in order; each adds its part: a row of the card (a label and a value) and its words in the sentence. The
// card's `text` is the sentence the announcer says, so a reader of the canvas and a reader of the live region are told
// the same things. The sections are internal (extensibility.md: chrome parts behind internal interfaces until an outside
// user needs one): not exported. One describer a viewer (`describerFor`, by its spec model): the announcer and the card
// share it, and its items read once a frame.

/** What the words read: the models a behavior is given, and the last frame's layout. */
export interface DescribeSource {
  readonly spec: SpecModel;
  readonly store: StoreModel;
  readonly addresses: AddressModel;
  /** The resolved view (a dragged rectangle's summary reads its grids). */
  readonly view: ViewModel;
  layout(): ViewLayout | undefined;
}

/** A rectangle dragged on a grid: its summary is described. */
export interface RectPlace {
  readonly rect: Location;
}

/**
 * The keyboard cursor's place (_spec/README.md "Input", Announcements): its cell, on a section of two axes with its row
 * (the row's label, a grid cell's two positions), or its row's header (the row alone).
 */
export interface CursorPlace {
  readonly cursor: Cursor;
}

/** A place to describe: the hover's record (`HoverRecord`) or a cell (the keyboard cursor's). */
export interface Place {
  readonly section: string | undefined;
  readonly layer?: string | undefined;
  readonly system: string | undefined;
  /** NaN when none. */
  readonly position: number;
  readonly element?: Key | undefined;
  readonly ySystem?: string | undefined;
  /** NaN when none. */
  readonly yPosition?: number;
  readonly target?: string | undefined;
  /** The block (its index in the section) and a grid cell's y block, when the section has blocks. */
  readonly block?: number | undefined;
  readonly yBlock?: number | undefined;
  /** What a layer named within the position (a rows column cell: `title: value`). */
  readonly state?: string | undefined;
  /** A merged run of mixed categories (a ranges layer's hit). */
  readonly categories?: string | undefined;
  /** A ring's grouped label: the hits of the items it stands for (`target` `group`). */
  readonly members?: readonly { readonly element?: Key | undefined }[] | undefined;
  /** A run's blocks (`target` `run`). */
  readonly count?: number | undefined;
}

/** The parts a description collects. */
class Parts {
  title = '';
  readonly rows: HoverCardRow[] = [];
  readonly words: string[] = [];
  /** Set by a section that says all there is (an item: its own words): the place's sections then say nothing. */
  whole = false;

  row(label: string, value: string, item = false): void {
    this.rows.push(item ? { label, value, item } : { label, value });
  }
}

/** What a section reads. */
interface PlaceContext {
  readonly words: Words;
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

/** A value as a card says it: up to three significant decimals, no trailing zeros; its unit after it when given. */
function valueText(v: number, unit = ''): string {
  return unit ? `${numberText(v)} ${unit}` : numberText(v);
}

function numberText(v: number): string {
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
    // a run of small blocks: the block under the point, which of the run, not read at this zoom (genome coordinates)
    if (at.place.target === 'run' || at.place.count !== undefined) {
      const block = blocks[k];
      const of = `one of ${(at.place.count ?? 0).toLocaleString('en-US')}`;
      out.title = block.label;
      out.row(blockKindName(at.model, block.axis.system), block.label);
      out.row('Run', of);
      if (at.place.state) out.row('Data', at.place.state);
      out.words.push(`${block.label}, ${of}${at.place.state ? `, ${at.place.state}` : ''}`);
      out.whole = true;
      return;
    }
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
    // the row is named by the block's system: an assembly's sequence (a chromosome, a contig) is a sequence, a
    // complex's unit or a chain's sequence a chain
    const name = blockKindName(at.model, blocks[k].axis.system);
    const y = at.place.yBlock;
    if (y !== undefined && y !== k) {
      say(k, `${name} x`);
      say(y, `${blockKindName(at.model, blocks[y].axis.system)} y`);
    } else say(k, name);
  },
};

/**
 * What a block's row is called by its system: `Chain` for a complex's (a unit, a molecule, its atoms) or a protein's,
 * `Sequence` for a nucleic sequence (an assembly's chromosome, a plain DNA or RNA sequence).
 */
function blockKindName(model: StaticModel, system: string): string {
  const e = model.systems.get(system);
  return !e || e.kind !== 'sequence' || e.alphabet === 'protein' ? 'Chain' : 'Sequence';
}

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
    // a complex's link: each end on its row's unit
    const as = (a.units && d.aUnit && a.units[d.aUnit[k]]) || a.target;
    const bs = (a.units && d.bUnit && a.units[d.bUnit[k]]) || (a.ySystem ?? a.target);
    const x = `${at.words.systemName(as)} ${spanText(d.aStart[k], d.aEnd[k])}`;
    const y = `${at.words.systemName(bs)} ${spanText(d.bStart[k], d.bEnd[k])}`;
    const kind = d.kinds[d.kind[k]];
    out.title = label;
    out.row('From', x);
    out.row('To', y);
    if (kind) out.row('Kind', kind);
    out.words.push(`${label}, ${x} to ${y}${kind ? `, ${kind}` : ''}`);
    out.whole = true;
  },
};

/**
 * A ring's grouped label (a dense fan as one): every member, each with its extent. The card is the whole list: a group
 * larger than `rings.cluster.fanMax` opens no fan, so no member is left to a "more".
 */
const GroupSection: CardSection = {
  name: 'group',
  describe(at, out) {
    const l = at.layer;
    const members = at.place.members;
    if (out.whole || at.place.target !== 'group' || !members?.length || !l?.type.items) return;
    const items = at.words.itemsOf(l);
    const system = at.words.systemName(l.layer.data?.system ?? at.section.section.axis.system);
    out.title = `${members.length} items`;
    const said: string[] = [];
    for (const m of members) {
      const it = m.element ? items.find((d) => d.key === m.element!.key) : undefined;
      if (!it) continue;
      out.row(it.label, `${system} ${extentText(it)}`, true);
      said.push(it.label);
    }
    out.words.push(`${members.length} items: ${said.join(', ')}`);
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
    const strand = at.words.strandOf(l, el.key, it.start);
    out.title = it.label;
    // the track's name: a short drawn label (`CDR-H1`) is told apart by it (`CDRs (Kabat)`)
    const track = at.section.tracks.find((t) => t.layers.includes(l))?.track.label;
    if (track) out.row('Track', track);
    out.row('Extent', `${system} ${extentText(it)}`);
    // a hit in a merged run of mixed categories lists them all (the hit's `categories`, in the legend's order)
    const mixed = at.place.categories;
    if (mixed) out.row('Categories', mixed);
    else if (category) out.row('Category', category);
    if (strand) out.row('Strand', strand);
    out.words.push(
      `${itemText(system, it)}${track ? `, ${track}` : ''}${mixed ? `, ${mixed}` : category ? `, ${category}` : ''}${strand ? `, strand ${strand}` : ''}`,
    );
    out.whole = true;
  },
};

/**
 * A rows layer's hidden-columns chip (core _spec "Render", The rows layer, Row columns, Narrow): each column it hides
 * and why.
 */
const ColumnsHiddenSection: CardSection = {
  name: 'columns-hidden',
  describe(at, out) {
    const l = at.layer;
    if (out.whole || at.place.target !== 'columns-hidden' || !l || l.type.name !== 'rows') return;
    const cols = (l.lay as RowsLay).cols;
    if (!cols.chipCols.length) return;
    const titles = cols.chipCols.map((k) => cols.columns[k].title);
    out.title = 'Hidden columns';
    for (const t of titles) out.row(t, HiddenWhy);
    out.words.push(`${titles.join(', ')} hidden, ${HiddenWhy}`);
    out.whole = true;
  },
};

/** Why a row column is hidden: the only reason there is (the view's width, EC-107). */
const HiddenWhy = 'no room at this width';

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
    // the row the keyboard cursor stands on: a folded row's summary row, labeled as drawn (`first +n`), a hidden one's
    // nearest shown (`cursorRow`)
    const i = lay.cursorRow(e);
    const label = i >= 0 ? lay.displayLabel(i) : lay.label(e);
    const shown = i >= 0 ? lay.display.rowAt(i) : e;
    // the title names the row: no row repeats it
    out.title = label;
    // a lazy row whose tile has not come yet: said pending (its label is its number until it does)
    const pending = !lay.held && shown >= 0 && !lay.lazyLabels.has(shown);
    out.words.push(`row ${label}${pending ? ', pending' : ''}`);
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
    const value = valueText(v, unitOf(l));
    out.title = `${x} × ${y}`;
    out.row('x', x);
    out.row('y', y);
    if (!Number.isNaN(v)) out.row('Value', value);
    out.words.push(`${at.section.id}, ${x} × ${y}${Number.isNaN(v) ? '' : `, value ${value}`}`);
    out.whole = true;
  },
};

/**
 * A letters layer on another system than its section's through a codon step (a protein's letters over its gene, a residue a codon wide): the
 * residue under the pointer, its letter, and the codon's positions on the section's system (a codon split by an intron
 * names each part). The place's own rows would say the nucleotide and the first letters layer's letter.
 */
const ResidueSection: CardSection = {
  name: 'residue',
  describe(at, out) {
    const l = at.layer;
    const data = l?.layer.data?.system;
    const system = at.place.system;
    if (out.whole || !l || l.type.name !== 'letters' || !data || !system || data === system) return;
    if (Number.isNaN(at.position)) return;
    // only a codon step reads a residue a codon wide: an alignment row, an offset or a pairwise map keeps the place's card
    const step = at.model.path(data, system).path?.step;
    if (!step) return;
    const residue = runsOf(
      at.model.project(rangesLocation(system, [at.position, at.position + 1]), data).location,
      data,
    );
    if (residue.length !== 2) return;
    const r = residue[0];
    const codon = runsOf(at.model.project(rangesLocation(data, [r, r + 1]), system).location, system);
    const name = at.words.systemName(data);
    const where = at.words.systemName(system);
    const letter = at.words.letterOn(data, r);
    const parts = Array.from({ length: codon.length / 2 }, (_, k) => spanText(codon[2 * k], codon[2 * k + 1])).join(
      ', ',
    );
    out.title = `${name} ${spanText(r, r + 1)}`;
    out.row('Index', String(r));
    if (letter) out.row('Letter', letter);
    if (parts) out.row('Codon', `${where} ${parts}`);
    out.words.push(
      `${at.section.id}, ${name} ${spanText(r, r + 1)}${letter ? `, letter ${letter}` : ''}${parts ? `, codon on ${where} ${parts}` : ''}`,
    );
    const items = at.words.itemsAt(at.section, at.position);
    for (const it of items) out.row(it.label, it.extent, true);
    if (items.length) out.words.push(`items: ${items.map((i) => `${i.label}, ${i.extent}`).join('; ')}`);
    out.whole = true;
  },
};

/** The ranges a location holds on a system, `[start, end, …]` (empty when none). */
function runsOf(loc: Location, system: string): ArrayLike<number> {
  return loc.parts.find((p) => p.system === system)?.ranges ?? [];
}

/** The row's value for a position a numbering has no label for (the ruler's note, `unnumbered`). */
const NotNumbered = 'not numbered';

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
    // a numbering with no label there is said as not numbering it (the ruler marks the same run), apart from the labels
    const labels = declared
      .filter(([, label]) => label !== NotNumbered)
      .map(([id, label]) => `${id} ${label}`)
      .join(', ');
    out.words.push(`${at.section.id}, ${name} ${ordinal}${labels ? ` (${labels})` : ''}`);
    for (const [id, label] of declared) if (label === NotNumbered) out.words.push(`not numbered by ${id}`);
    const letter = w.letterAt(at.section, position);
    if (letter) {
      out.row('Letter', letter);
      out.words.push(`letter ${letter}`);
    }
    // a modified residue's component, its alternatives; a molecule's residue; an atom and its residue
    for (const r of chemistryRows(at.model, system, position)) {
      out.row(r.label, r.value);
      out.words.push(r.words);
    }
    for (const [label, value] of w.valuesAt(at.section, position)) {
      out.row(label, value);
      out.words.push(`${label} ${value}`);
    }
    const items = w.itemsAt(at.section, position);
    for (const it of items) out.row(it.label, it.extent, true);
    if (items.length) out.words.push(`items: ${items.map((i) => `${i.label}, ${i.extent}`).join('; ')}`);
  },
};

// the block last: what is in it is said first, then its chain
const Sections: readonly CardSection[] = [
  PairSection,
  GroupSection,
  ItemSection,
  ColumnsHiddenSection,
  RowSection,
  GridSection,
  ResidueSection,
  PlaceSection,
  BlockSection,
];

/**
 * A unit the layer's theme declares for its values (a `field` form's or a scale's `unit`, `Å`): the first found in its
 * channels' forms; '' when none does.
 */
function unitOf(l: LayerLayout): string {
  const find = (v: unknown, depth: number): string => {
    if (depth > 4 || typeof v !== 'object' || v === null) return '';
    const node = v as { params?: Record<string, unknown>; args?: readonly unknown[] };
    if (typeof node.params?.unit === 'string') return node.params.unit;
    for (const x of [...Object.values(node.params ?? {}), ...(node.args ?? [])]) {
      const u = find(x, depth + 1);
      if (u) return u;
    }
    return '';
  };
  for (const v of Object.values(l.layer.theme ?? {})) {
    const u = find(v, 0);
    if (u) return u;
  }
  return '';
}

/** What the sections read, and the description itself (internal: the `Describer` is its face). */
class Words {
  private readonly _src: DescribeSource;
  /** Items by layer layout, made once a frame (a description looks them up, a frame changes them). */
  private _items = new Map<LayerLayout, readonly ItemDescription[]>();
  private _itemsLayout: unknown;

  constructor(src: DescribeSource) {
    this._src = src;
  }

  /** The header words (the `headers` behavior's), when it is listed. */
  headers: ((place: Place) => { readonly title: string; readonly text: string } | undefined) | undefined = undefined;

  /** The description of a place in the last frame; undefined when it is in none of its sections. */
  describe(place: Place): HoverCard | undefined {
    // a header's part: its entry's words alone (the card a title alone; _spec/README.md "Behaviors", Headers)
    const t = place.target;
    if (t === 'header-slot' || t === 'header-state' || t === 'header-label' || t === 'section-header') {
      const h = this.headers?.(place);
      return h ? { title: h.title, rows: [], text: h.text } : undefined;
    }
    const layout = this._src.layout();
    if (!layout || place.section === undefined) return undefined;
    const section = layout.sections.find((s) => s.section.address === place.section);
    if (!section) return undefined;
    const layer = place.layer === undefined ? undefined : layerOf(section, place.layer);
    const out = new Parts();
    const at: PlaceContext = {
      words: this,
      model: this._src.view.model.value,
      section,
      layer,
      place,
      position: Number.isNaN(place.position) ? Number.NaN : Math.floor(place.position),
      yPosition:
        place.yPosition === undefined || Number.isNaN(place.yPosition) ? Number.NaN : Math.floor(place.yPosition),
    };
    // a run of small blocks says its block and the run alone
    for (const s of Sections) if ((t !== 'run' && place.count === undefined) || s === BlockSection) s.describe(at, out);
    if (!out.words.length) return undefined;
    return { title: out.title || out.words[0], rows: out.rows, text: out.words.join(', ') };
  }

  /**
   * The cursor as a place: on a section of two axes its row's layer and row, its header the row alone; on a track's
   * header the header (`state: 'cursor'`: its words say the key that opens its actions).
   */
  cursorPlace(c: Cursor): Place {
    const cell: Place = { section: c.section, system: c.system, position: c.position };
    if (c.track !== undefined)
      return {
        ...cell,
        layer: c.track,
        target: c.track === c.section ? 'section-header' : 'header-label',
        state: 'cursor',
      };
    const layout = this._src.layout();
    if (c.ySystem === undefined || c.yPosition === undefined || !layout) return cell;
    const s = layout.sections.find((x) => x.section.address === c.section);
    const l = s && twoAxesOf(s, c.ySystem);
    if (!l) return cell;
    const y = { layer: l.key, ySystem: c.ySystem, yPosition: c.yPosition };
    if (c.part === 'header')
      return { section: c.section, system: c.ySystem, position: c.yPosition, ...y, target: 'row' };
    return { ...cell, ...y, target: 'cell' };
  }

  /** A rectangle on a grid: its cells counted, their mean, least and greatest value (`gridSummary`). */
  rect(loc: Location): HoverCard | undefined {
    const view = this._src.view.resolved.value;
    if (!loc.rects.length || !view) return undefined;
    const s = gridSummary({ model: this._src.view.model.value, view }, loc);
    if (!s.count) return undefined;
    const rows = [
      { label: 'Mean', value: numberText(s.mean) },
      { label: 'Min', value: numberText(s.min) },
      { label: 'Max', value: numberText(s.max) },
    ];
    const cells = `${s.count.toLocaleString('en-US')} cell${s.count === 1 ? '' : 's'}`;
    const text = `Rectangle: ${cells}, ${rows.map((r) => `${r.label.toLowerCase()} ${r.value}`).join(', ')}`;
    return { title: cells, rows, text };
  }

  // ---- what the sections read ----

  systemName(system: string): string {
    return this._src.view.model.value.systems.get(system)?.name ?? system;
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

  /**
   * An item's strand (`+` or `−`) when its data has one: the `strand` field of an inline annotation's item, else the
   * strand column of the loaded list tile of its source holding the item's first position; '' when it has none or its
   * tile is not loaded.
   */
  strandOf(l: LayerLayout, key: string | number, start: number): string {
    const model = this._src.view.model.value;
    const owner = l.layer.data?.address;
    const a = owner ? model.annotations.get(owner) : undefined;
    const f = a?.items?.find((i) => i.key === key)?.fields?.strand;
    if (typeof f === 'string') return f === '-' ? '−' : f === '+' ? '+' : '';
    if (typeof f === 'number') return f < 0 ? '−' : f > 0 ? '+' : '';
    const src = layerSource(model, l.layer);
    if (!src) return '';
    const store = this._src.store.store;
    // the layer reads the level its pixels need: the list tile holding the item is at one of the source's levels
    for (const level of store.coverage(src.source)?.levels ?? [1])
      for (const k of store.keysOf({ source: src.source, system: src.system, start, end: start + 1, level })) {
        const e = store.get(k);
        const t = e?.status === 'loaded' ? e.tables.ranges : undefined;
        if (!t || t.form !== 'list' || !t.strand) continue;
        for (let i = 0; i < t.start.length; i++) {
          if (t.start[i] !== start || (t.id ? t.id.get(i) : t.row?.[i]) !== key) continue;
          return t.strand[i] < 0 ? '−' : t.strand[i] > 0 ? '+' : '';
        }
      }
    return '';
  }

  /** An item's category: the `category` field of an inline annotation's item (a source's items are not held here). */
  categoryOf(l: LayerLayout, key: string | number): string {
    const owner = l.layer.data?.address;
    const a = owner ? this._src.view.model.value.annotations.get(owner) : undefined;
    const f = a?.items?.find((i) => i.key === key)?.fields?.category;
    return typeof f === 'string' ? f : '';
  }

  /**
   * The items of a section at a position: its tracks' layers' items under the position's x; an item across a cyclic
   * origin (its two `runs`) by its runs, never by its drawn x span (the complement of what it covers).
   */
  itemsAt(section: SectionLayout, position: number, only?: LayerLayout): { label: string; extent: string }[] {
    const axis = section.axis;
    if (!axis) return [];
    const a = axis.x(position);
    const b = axis.x(position + 1);
    if (!Number.isFinite(a) || !Number.isFinite(b)) return [];
    const x = (a + b) / 2;
    const axisSystem = section.section.axis.system;
    const out: { label: string; extent: string }[] = [];
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (only && l !== only) continue;
        const data = l.layer.data?.system;
        const own = (it: ItemDescription) => (it.system ?? data ?? section.section.axis.system) === axisSystem;
        for (const it of this.itemsOf(l))
          if (it.y === undefined && (it.runs ? own(it) && inRuns(it.runs, position) : it.x0 <= x && x <= it.x1))
            out.push({
              label: it.label,
              extent: `${this.systemName(it.system ?? data ?? section.section.axis.system)} ${extentText(it)}`,
            });
      }
    return out;
  }

  /** The letter a letters layer of the section shows at a position (none while it is not loaded). */
  letterAt(section: SectionLayout, position: number): string {
    const model = this._src.view.model.value;
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
    const s = this._src.view.model.value.systems.get(system);
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
          const label = t.track.label || l.key;
          const unit = unitOf(l);
          out.push([
            label,
            lo === hi || Number.isNaN(lo) ? valueText(hi, unit) : `${numberText(lo)}–${valueText(hi, unit)}`,
          ]);
          break;
        }
      }
    return out;
  }

  /** A grid's value at a cell (the hovered layer's grid, on its two systems), NaN when none. */
  gridValue(l: LayerLayout, system: string, i: number, ySystem: string, j: number): number {
    const owner = l.layer.data?.address;
    const a = owner ? this._src.view.model.value.annotations.get(owner) : undefined;
    const g = a?.grid;
    if (!a || !g || a.target !== system || (a.ySystem ?? a.target) !== ySystem) return Number.NaN;
    if (i < 0 || j < 0 || i >= g.nx || j >= g.ny) return Number.NaN;
    return g.values[i * g.ny + j];
  }

  /**
   * The position's labels in the system's declared numberings, with their ids: the numbering the section's layers show
   * first (as a ruler or a text's numbers do), then the others; one with no label there (absent, or not arrived) says
   * `–`, so a reader sees it is missing.
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
      // absent: `–` for a numbering not arrived; a numbering that is there and has no label here does not number it
      out.push([id, label ?? (n ? NotNumbered : '–')]);
    }
    return out;
  }
}

/** The words of places in a viewer: the hover's, the cursor's, a dragged rectangle's. */
export class Describer {
  private readonly _words: Words;

  constructor(src: DescribeSource) {
    this._words = new Words(src);
  }

  /**
   * The description of a place, a rectangle's summary or the keyboard cursor's place in the last frame; undefined when
   * there is none.
   */
  describe(place: Place | RectPlace | CursorPlace): HoverCard | undefined {
    if ('rect' in place) return this._words.rect(place.rect);
    return this._words.describe('cursor' in place ? this._words.cursorPlace(place.cursor) : place);
  }

  /** The header words of a place (a slot, a chip, a label): the `headers` behavior sets them when it is listed. */
  get headers(): Words['headers'] {
    return this._words.headers;
  }
  set headers(fn: Words['headers']) {
    this._words.headers = fn;
  }

  /** The text of a place ('' when there is none). */
  text(place: Place | RectPlace | CursorPlace): string {
    return this.describe(place)?.text ?? '';
  }
}

const Describers = new WeakMap<SpecModel, Describer>();

/** The viewer's one describer (by its spec model): the announcer's and the hover card's. */
export function describerFor(src: DescribeSource): Describer {
  let d = Describers.get(src.spec);
  if (!d) Describers.set(src.spec, (d = new Describer(src)));
  return d;
}

/** The layer of two axes a section shows `ySystem` with (its rows layer, its grid): the cursor's row's. */
function twoAxesOf(section: SectionLayout, ySystem: string): LayerLayout | undefined {
  const is = (l: LayerLayout) =>
    (l.type.name === 'rows' && (l.lay as RowsLay).rowsSystem === ySystem) ||
    (l.type.name === 'grid' && (l.lay as GridLay).ySystem === ySystem);
  for (const t of section.tracks) for (const l of t.layers) if (is(l)) return l;
  for (const ls of section.lays) for (const l of ls) if (is(l)) return l;
  for (const l of section.sectionLays) if (is(l)) return l;
  return undefined;
}

/** A layer of a section by its key: its tracks', its strips' (a section of blocks), its section layers'. */
function layerOf(section: SectionLayout, key: string): LayerLayout | undefined {
  for (const t of section.tracks) for (const l of t.layers) if (l.key === key) return l;
  for (const ls of section.lays) for (const l of ls) if (l.key === key) return l;
  for (const l of section.sectionLays) if (l.key === key) return l;
  return undefined;
}

/** Whether a position is in one of `[start, end)` runs (an item's across a cyclic origin). */
function inRuns(runs: readonly number[], position: number): boolean {
  for (let k = 0; k + 1 < runs.length; k += 2) if (position >= runs[k] && position < runs[k + 1]) return true;
  return false;
}

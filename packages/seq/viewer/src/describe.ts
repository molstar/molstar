/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesList } from '@molstar/seqx-core/data';
import type { Key, Location } from '@molstar/seqx-core/query';
import { DefaultNumbering, findNumbering, numberingIds } from '@molstar/seqx-core/query';
import { extentText } from '@molstar/seqx-core/render/accessibility';
import {
  type GridLay,
  layerSource,
  type ItemDescription,
  type LayerLayout,
  LettersType,
  type RowsLay,
  type SectionLayout,
  shownLetters,
  gridSummary,
  type ItemTitle,
  type ValuesLay,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import type { AddressModel, Cursor, SpecModel, StoreModel, ViewModel } from '@molstar/seqx-viewer/models';
import { issueOf, type Issue } from '@molstar/core/repr/tree';
import { BridgeView } from '@molstar/seqx-viewer/models/bridge-hooks';
import {
  type CardAt,
  type CardParts,
  type CardProviderEntry,
  type CardReads,
  composeCardProviders,
} from './card-provider.js';
import {
  BuiltInCardProviders,
  CardProviderNames,
  NotNumbered,
  numberText,
  unitOf,
  valueText,
} from './card-providers.js';
import type { HoverCard, HoverCardRow } from './hover-card.js';
import { ViewerIssues } from './issues.js';

// What a place is, in words (_spec/README.md "Input", Announcements; "Hover card"): one description that the live
// region says and the hover card shows. A place (the hover's record, the keyboard cursor's cell) goes through the
// card's providers in order (_spec/README.md "Hover card", Providers: the registry's `card-provider` entries, the
// built-ins `BuiltInCardProviders`); each adds its part: rows of the card (a label and a value) and its words in the
// sentence. The card's `text` is the sentence the announcer says, so a reader of the canvas and a reader of the live
// region are told the same things. One describer a viewer (`describerFor`, by its spec model): the announcer and the
// card share it, and its items read once a frame.

/** The built-in providers composed: what a describer runs until its viewer sets its registry's. */
const DefaultProviders = composeCardProviders(BuiltInCardProviders);

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
class Parts implements CardParts {
  title = '';
  readonly rows: HoverCardRow[] = [];
  readonly words: string[] = [];
  /** Set by a section that says all there is (an item: its own words): the place's sections then say nothing. */
  whole = false;
  moreCount = 0;

  row(label: string, value: string, item = false): void {
    this.rows.push(item ? { label, value, item } : { label, value });
  }
}

/** What the providers read, and the description itself (internal: the `Describer` is its face). */
class Words implements CardReads {
  private readonly _src: DescribeSource;
  /** Items by layer layout, made once a frame (a description looks them up, a frame changes them). */
  private _items = new Map<LayerLayout, readonly ItemDescription[]>();
  private _itemsLayout: unknown;
  /** The providers composed (`composeCardProviders`): the built-ins until the viewer sets its registry's. */
  providers: readonly CardProviderEntry[] = DefaultProviders;
  /** The `hover.cardItems` cap the providers listing matches read. */
  cardItems = 3;
  /** Says a provider that threw (its rows dropped), once a provider. */
  report: ((name: string, error: unknown) => void) | undefined = undefined;
  private readonly _threw = new Set<string>();

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
    // a ribbon in a bridge's band (its `section` the bridge's id): the bridge chunk's card
    if (!section) return BridgeView.parts?.card(this._src.view, layout, place);
    const layer = place.layer === undefined ? undefined : layerOf(section, place.layer);
    const out = new Parts();
    const at: CardAt = {
      words: this,
      model: this._src.view.model.value,
      view: this._src.view.resolved.value,
      store: this._src.store,
      section,
      layer,
      place,
      position: Number.isNaN(place.position) ? Number.NaN : Math.floor(place.position),
      yPosition:
        place.yPosition === undefined || Number.isNaN(place.yPosition) ? Number.NaN : Math.floor(place.yPosition),
      cardItems: this.cardItems,
    };
    // a run of small blocks says its block and the run alone
    const run = t === 'run' || place.count !== undefined;
    for (const p of this.providers) if (!run || p.name === CardProviderNames.block) this._run(p, at, out);
    if (!out.words.length) return undefined;
    const card = { title: out.title || out.words[0], rows: out.rows, text: out.words.join(', ') };
    return out.moreCount ? { ...card, moreCount: out.moreCount } : card;
  }

  /** One provider's part: a throw drops what it added (its rows, words and title), said once a provider. */
  private _run(p: CardProviderEntry, at: CardAt, out: Parts): void {
    const rows = out.rows.length;
    const words = out.words.length;
    const { title, whole, moreCount } = out;
    try {
      p.describe(at, out);
    } catch (e) {
      out.rows.length = rows;
      out.words.length = words;
      out.title = title;
      out.whole = whole;
      out.moreCount = moreCount;
      if (!this._threw.has(p.name)) {
        this._threw.add(p.name);
        this.report?.(p.name, e);
      }
    }
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
  strandOf(
    l: LayerLayout,
    key: string | number,
    start: number,
    /** Its source's row when already looked up (`sourceRowOf`; null: none). */
    row?: { t: RangesList; i: number } | null,
  ): string {
    const model = this._src.view.model.value;
    const owner = l.layer.data?.address;
    const a = owner ? model.annotations.get(owner) : undefined;
    const f = a?.items?.find((i) => i.key === key)?.fields?.strand;
    if (typeof f === 'string') return f === '-' ? '−' : f === '+' ? '+' : '';
    if (typeof f === 'number') return f < 0 ? '−' : f > 0 ? '+' : '';
    const r = row === undefined ? this.sourceRowOf(l, key, start) : row;
    const s = r?.t.strand?.[r.i] ?? 0;
    return s < 0 ? '−' : s > 0 ? '+' : '';
  }

  /**
   * An item's row in the loaded list tile of its layer's source holding its first position (its source's fields: a VCF
   * allele's); undefined for an inline annotation's item, or while no tile holding it is loaded.
   */
  sourceRowOf(l: LayerLayout, key: string | number, start: number): { t: RangesList; i: number } | undefined {
    const src = layerSource(this._src.view.model.value, l.layer);
    if (!src) return undefined;
    const store = this._src.store.store;
    // the layer reads the level its pixels need: the list tile holding the item is at one of the source's levels
    for (const level of store.coverage(src.source)?.levels ?? [1])
      for (const k of store.keysOf({ source: src.source, system: src.system, start, end: start + 1, level })) {
        const e = store.get(k);
        const t = e?.status === 'loaded' ? e.tables.ranges : undefined;
        if (!t || t.form !== 'list') continue;
        for (let i = 0; i < t.start.length; i++)
          if (t.start[i] === start && (t.id ? t.id.get(i) : t.row?.[i]) === key) return { t, i };
      }
    return undefined;
  }

  /** An address as positions in the spec shown, once it resolves (the address model's). */
  locate(address: string): Location | undefined {
    const r = this._src.addresses.locate(address);
    return r.resolution.status === 'resolved' ? r.location : undefined;
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

  /** Sets the providers the card is composed by, the cap they read, and where a throw is said (`composeCard`). */
  compose(
    providers: readonly CardProviderEntry[],
    cardItems: number,
    report: (name: string, error: unknown) => void,
  ): void {
    this._words.providers = providers;
    this._words.cardItems = cardItems;
    this._words.report = report;
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

/**
 * An export's `card` titles from a describer (_spec/README.md "Hover card", Export titles): each item's card, its title
 * and its rows as plain text, a line each, as the composed providers give it at the item's middle.
 */
export function cardTitles(describer: Describer): ItemTitle {
  return (s, l, it) => {
    const owner = l.layer.data?.address;
    const x = (it.x0 + it.x1) / 2;
    const position = s.axis ? s.axis.positionAt(x) : Number.NaN;
    const card = describer.describe({
      section: s.section.address,
      layer: l.key,
      system: s.section.axis.system,
      position,
      ...(owner !== undefined ? { element: { owner, key: it.key } } : {}),
    });
    if (!card) return undefined;
    return [card.title, ...card.rows.map((r) => (r.value ? `${r.label}: ${r.value}` : r.label))].join('\n');
  };
}

/**
 * A viewer's card composed from its registry (_spec/README.md "Hover card", Providers): its `card-provider` entries by
 * order, then name, a later one under a name replacing; a registry listing none keeps the built-ins (a spec written
 * before providers keeps its card). A composition without `seq:item` or `seq:place`, the hover contract's floor, is
 * said (`seq.viewer.card.floor`); a provider that throws has its rows dropped, said once (`seq.viewer.card.provider`).
 */
export function composeCard(
  describer: Describer,
  registry: readonly { readonly kind: string }[],
  cardItems: number,
  report: (issue: Issue) => void,
): void {
  const listed = registry.filter((e): e is CardProviderEntry => e.kind === 'card-provider');
  const providers = listed.length ? composeCardProviders(listed) : DefaultProviders;
  const missing = [CardProviderNames.item, CardProviderNames.place].filter((n) => !providers.some((p) => p.name === n));
  if (missing.length)
    report(
      issueOf(ViewerIssues, 'seq.viewer.card.floor', `The card's providers leave out ${missing.join(' and ')}`, {
        hint: ViewerIssues.codes['seq.viewer.card.floor'].hint,
      }),
    );
  describer.compose(providers, cardItems, (name, e) =>
    report(
      issueOf(
        ViewerIssues,
        'seq.viewer.card.provider',
        `The card provider ${name} threw: its rows are left out (${e instanceof Error ? e.message : String(e)})`,
        { hint: ViewerIssues.codes['seq.viewer.card.provider'].hint },
      ),
    ),
  );
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

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/seqx-core/query';
import { rangesLocation } from '@molstar/seqx-core/query';
import { extentText } from '@molstar/seqx-core/render/accessibility';
import { itemText, type LayerLayout, type RowsLay, spanText } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { CardProviderEntry } from './card-provider.js';
import { alleleOf } from './describe-allele.js';
import { chemistryRows } from './describe-molecules.js';
import { TooltipsCardProvider } from './card-tooltips.js';
import { regionHoverOf } from './describe-region.js';

// The card's built-in providers (_spec/README.md "Hover card", Providers): the viewer's own sections, registry entries
// of kind `card-provider` like a host's, ordered by `order` then name. Each adds its part of a place's description: rows
// of the card and its words in the sentence. Their orders keep the card's order (what is under the point first, the
// place, the spec's tooltips, the block last); a host's entry sits between them by its `order`.

/** The built-in providers' names. */
export const CardProviderNames = {
  region: 'seq:region',
  pair: 'seq:pair',
  group: 'seq:group',
  item: 'seq:item',
  columnsHidden: 'seq:columns-hidden',
  row: 'seq:row',
  grid: 'seq:grid',
  residue: 'seq:residue',
  place: 'seq:place',
  tooltips: 'seq:tooltips',
  block: 'seq:block',
} as const;

/** A value as a card says it: up to three significant decimals, no trailing zeros; its unit after it when given. */
export function valueText(v: number, unit = ''): string {
  return unit ? `${numberText(v)} ${unit}` : numberText(v);
}

export function numberText(v: number): string {
  if (!Number.isFinite(v)) return String(v);
  if (Number.isInteger(v)) return v.toLocaleString('en-US');
  const a = Math.abs(v);
  const digits = a >= 100 ? 1 : a >= 1 ? 2 : 3;
  return Number(v.toFixed(digits)).toString();
}

/**
 * A unit the layer's theme declares for its values (a `field` form's or a scale's `unit`, `Å`): the first found in its
 * channels' forms; '' when none does.
 */
export function unitOf(l: LayerLayout): string {
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

// ---- the built-ins ----

/** A block of a section of blocks (a chain of a complex): its label, and which copy when its system has several. */
const BlockSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.block,
  order: 110,
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
    // a grid over two sections of blocks: its y blocks are the second section's (EC-209)
    const yBlocks = at.section.section.y?.blocks ?? blocks;
    const say = (b: number, which: string, of = blocks) => {
      const block = of[b];
      if (!block) return;
      out.row(which, block.label);
      let copy = '';
      if (block.members.length > 1) copy = `${block.members.length} copies merged`;
      else {
        const twins = of.filter((o) => o.axis.system === block.axis.system);
        if (twins.length > 1) copy = `copy ${twins.indexOf(block) + 1} of ${twins.length}`;
      }
      if (copy) out.row('Copy', copy);
      out.words.push(`${which.toLowerCase()} ${block.label}${copy ? `, ${copy}` : ''}`);
    };
    // the row is named by the block's system: an assembly's sequence (a chromosome, a contig) is a sequence, a
    // complex's unit or a chain's sequence a chain
    const name = blockKindName(at.model, blocks[k].axis.system);
    const y = at.place.yBlock;
    if (y !== undefined && yBlocks[y] && (y !== k || yBlocks !== blocks)) {
      say(k, `${name} x`);
      say(y, `${blockKindName(at.model, yBlocks[y].axis.system)} y`, yBlocks);
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

/** Data read by region (a pairs bin, a region grid's cell, a matrix cell): the core's words (`describe-region.ts`). */
const RegionSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.region,
  order: 10,
  describe(at, out) {
    const h = regionHoverOf(at.layer, at.place.element, at.position, at.yPosition, (s) => at.words.systemName(s));
    if (!h) return;
    out.title = h.title;
    for (const [label, value] of h.rows) out.row(label, value);
    out.words.push(h.text);
    out.whole = true;
  },
};

/** A pairs layer's pair: its label, its two ends, its kind. */
const PairSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.pair,
  order: 20,
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
const GroupSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.group,
  order: 30,
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
const ItemSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.item,
  order: 40,
  describe(at, out) {
    const el = at.place.element;
    const l = at.layer;
    if (out.whole || !el || !l?.type.items) return;
    const it = at.words.itemsOf(l).find((d) => d.key === el.key);
    if (!it) return;
    const system = at.words.systemName(it.system ?? l.layer.data?.system ?? at.section.section.axis.system);
    // a source's allele (a VCF's: its row holds `ref` and `alt`) is named by its id, `pos:ref>alt`, never its record's ID
    // (its source's row looked up once: the strand reads it too)
    const row = at.words.sourceRowOf(l, el.key, it.start);
    const allele = row && it.system === undefined ? alleleOf(row.t, row.i) : undefined;
    const category = allele?.category || at.words.categoryOf(l, el.key);
    const strand = at.words.strandOf(l, el.key, it.start, row ?? null);
    out.title = allele ? String(el.key) : it.label;
    // the track's name: a short drawn label (`CDR-H1`) is told apart by it (`CDRs (Kabat)`)
    const track = at.section.tracks.find((t) => t.layers.includes(l))?.track.label;
    if (track) out.row('Track', track);
    out.row('Extent', `${system} ${extentText(it)}`);
    if (allele) {
      if (category) out.row('Category', category);
      for (const [name, value] of allele.rows) out.row(name, value);
    }
    // a hit in a merged run of mixed categories lists them all (the hit's `categories`, in the legend's order)
    const mixed = at.place.categories;
    if (mixed) out.row('Categories', mixed);
    else if (category && !allele) out.row('Category', category);
    if (strand) out.row('Strand', strand);
    out.words.push(
      `${allele ? `${String(el.key)}, ` : ''}${itemText(system, it)}${track ? `, ${track}` : ''}${mixed ? `, ${mixed}` : category ? `, ${category}` : ''}${strand ? `, strand ${strand}` : ''}`,
    );
    out.whole = true;
  },
};

/**
 * A rows layer's hidden-columns chip (core _spec "Render", The rows layer, Row columns, Narrow): each column it hides
 * and why.
 */
const ColumnsHiddenSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.columnsHidden,
  order: 50,
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
const RowSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.row,
  order: 60,
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
const GridSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.grid,
  order: 70,
  describe(at, out) {
    const p = at.place;
    const l = at.layer;
    if (out.whole || !l || l.type.name !== 'grid' || p.system === undefined || p.ySystem === undefined) return;
    if (Number.isNaN(at.position) || Number.isNaN(at.yPosition)) return;
    // over blocks (a grid over a complex) each axis by its block's label, the chain: `CDK2 33 × cyclin A2 263`
    const blocks = at.section.section.blocks;
    const yBlocks = at.section.section.y?.blocks ?? blocks;
    const nameOf = (system: string, k: number | undefined, of = blocks) =>
      (k !== undefined && of?.[k]?.axis.system === system ? of[k].label : undefined) ?? at.words.systemName(system);
    // a residue as a structure's readout writes it: its letter before its number (`CDK2 K33`), when it is loaded
    const end = (system: string, k: number | undefined, i: number, of = blocks) =>
      `${nameOf(system, k, of)} ${at.words.letterOn(system, i)}${spanText(i, i + 1)}`;
    // a categorical axis (EC-35): its member by label (`NA12878 × NA12891`)
    const st = at.section.placement.strips;
    const xm = st.n ? st.axes[0].members : undefined;
    const ym = st.n ? st.ys[0]?.members : undefined;
    const x = xm?.labels[at.position] ?? end(p.system, p.block, at.position);
    const y = ym?.labels[at.yPosition] ?? end(p.ySystem, p.yBlock ?? p.block, at.yPosition, yBlocks);
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
const ResidueSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.residue,
  order: 80,
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
export const NotNumbered = 'not numbered';

/**
 * A cell: the section, the system and the position in every numbering of the system (the ordinal, the index, the
 * declared ones: the section's shown numbering first), the letter shown there, a value track's value and the items
 * under it.
 */
const PlaceSection: CardProviderEntry = {
  kind: 'card-provider',
  name: CardProviderNames.place,
  order: 90,
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

/**
 * The built-in providers (_spec/README.md "Hover card", Providers), in the default spec's registry; a host leaves one
 * out of the registry it spreads to remove it, or lists its own under the same name to replace it.
 */
export const BuiltInCardProviders: readonly CardProviderEntry[] = [
  RegionSection,
  PairSection,
  GroupSection,
  ItemSection,
  ColumnsHiddenSection,
  RowSection,
  GridSection,
  ResidueSection,
  PlaceSection,
  TooltipsCardProvider,
  BlockSection,
];

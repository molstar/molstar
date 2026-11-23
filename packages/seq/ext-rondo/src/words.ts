/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayerPlace, PairsSource } from '@molstar/seqx-core/extension';
import { foldedOf, Kind, type RondoLay, stemAt } from './lay.js';

// What the ribbons say (_spec/README.md, Hover, the keyboard): the card's lines for a hit (core's `LayerType.describe`:
// the viewer's own text after), the layer's ARIA summary, and what the keyboard cursor adds at a stem.

/** A length in bases as a reader says it: `bp`, `kb`, `Mb` (two significant decimals at most). */
export function basesText(n: number): string {
  const unit = n >= 1e6 ? ['Mb', 1e6] : n >= 1e3 ? ['kb', 1e3] : ['bp', 1];
  const v = n / (unit[1] as number);
  return `${Number(v.toFixed(v >= 100 || unit[1] === 1 ? 0 : v >= 10 ? 1 : 2))} ${unit[0]}`;
}

/** A position as a coordinate (`50` Mb at most two decimals, a bare number below a kb). */
function at(n: number): string {
  if (n >= 1e6) return `${Number((n / 1e6).toFixed(2))}`;
  if (n >= 1e3) return `${Number((n / 1e3).toFixed(2))}`;
  return `${n}`;
}

/** A stem's range on its system: `chr4:50–58 Mb`. */
function rangeText(lay: RondoLay, s: number): string {
  const start = lay.stemStart[s];
  const end = lay.stemEnd[s];
  const unit = end >= 1e6 ? ' Mb' : end >= 1e3 ? ' kb' : '';
  return `${lay.names[lay.stemSystem[s]]}:${at(start)}–${at(end)}${unit}`;
}

const links = (n: number) => `${Math.round(n).toLocaleString('en-US')} link${Math.round(n) === 1 ? '' : 's'}`;

/** A field's value as the card says it: a number to 4 significant digits, a strand's sign as a minus. */
const valueText = (v: string | number) => (typeof v === 'number' ? String(Number(v.toPrecision(4))) : v);

/**
 * A link's fields (_spec/README.md, Kinds of links, The card): its kind, its ends' strands, then every other field the
 * data carries for it (`q`, gene names where the file has them), absent values left out.
 */
function fieldLines(lay: RondoLay, k: number): string[] {
  const read = lay.read as PairsSource | undefined;
  const d = read?.data;
  const datum = read?.datum;
  if (!d || !datum) return [];
  datum.set(k);
  const out: string[] = [];
  const type = datum.field('type');
  if (type !== undefined) out.push(`Type: ${type}`);
  out.push(`Strands: ${datum.field('a_strand') === '-' ? '−' : '+'} ${datum.field('b_strand') === '-' ? '−' : '+'}`);
  for (const name of Object.keys(d.extra ?? {})) {
    if (name === 'count') continue;
    const v = datum.field(name);
    if (v !== undefined && v !== null && v !== '' && !(typeof v === 'number' && Number.isNaN(v)))
      out.push(`${name}: ${valueText(v)}`);
  }
  return out;
}

/**
 * Each layer's links on the section between two ranges (`b` absent: with an end in `a`), one line a layer
 * (`Contacts: 140 links`) when the section has more than one: a bin's or a stem's count per kind. Counted from each
 * layer's rows as read (bins' starts, items' ends), each row once.
 */
function perLayer(
  lay: RondoLay,
  sysA: string,
  a0: number,
  a1: number,
  sysB?: string,
  b0 = 0,
  b1 = 0,
): string[] {
  if (lay.peers.length < 2) return [];
  const out: string[] = [];
  for (const p of lay.peers) {
    const read = p.read;
    const d = read?.data;
    if (!read || !d || p.empty) continue;
    const sysOf = (unit: Int32Array | undefined, k: number, side: 0 | 1) =>
      unit ? p.systems[unit[k]] : p.systems[Math.min(side, p.systems.length - 1)];
    let n = 0;
    for (let k = 0; k < d.n; k++) {
      if (!read.shown[k]) continue;
      const sa = sysOf(d.aUnit, k, 0);
      const sb = sysOf(d.bUnit, k, 1);
      const inA = (s: string, x: number) => s === sysA && x >= a0 && x < a1;
      const inB = (s: string, x: number) => s === sysB && x >= b0 && x < b1;
      const hit =
        sysB === undefined
          ? inA(sa, d.aStart[k]) || inA(sb, d.bStart[k])
          : (inA(sa, d.aStart[k]) && inB(sb, d.bStart[k])) || (inA(sb, d.bStart[k]) && inB(sa, d.aStart[k]));
      if (hit) n += read.weight ? read.weight[k] : 1;
    }
    out.push(`${p.label || 'Links'}: ${links(n)}`);
  }
  return out;
}

/** The card's lines for a hit of the layer: a ribbon's two stems and count, a stem's, a stub's partner. */
export function describeRondo(hit: LayerPlace, lay: RondoLay): readonly string[] | undefined {
  const key = hit.element?.key;
  if (typeof key !== 'string') return undefined;
  const i = lay.byKey.get(key);
  if (i === undefined) return undefined;
  const kind = lay.polyKind[i];
  const ref = lay.polyRef[i];
  const bins = lay.items ? undefined : basesText(lay.level);
  if (kind === Kind.ribbon) {
    const a = lay.ribA[ref];
    const b = lay.ribB[ref];
    const kindName = lay.kinds[lay.ribKind[ref]];
    const lines = [
      `${rangeText(lay, a)} ↔ ${rangeText(lay, b)}`,
      `Links: ${links(lay.ribLinks[ref])}${kindName ? ` (${kindName})` : ''}`,
    ];
    if (lay.width !== 'count') lines.push(`Width: ${lay.width} ${valueText(lay.ribCount[ref])}`);
    if (lay.ribStrand[ref] < 0) lines.push('Strand: reversed');
    // one link: its own fields (its kind, strands, q, gene names); a bin: its count on each layer of the section
    if (lay.items && lay.ribLinks[ref] === 1) lines.push(...fieldLines(lay, lay.ribRow[ref]));
    if (bins) {
      lines.push(`Bins: ${bins} bins, zoom in to see them`);
      const sa = lay.systems[lay.stemSystem[a]];
      const sb = lay.systems[lay.stemSystem[b]];
      lines.push(...perLayer(lay, sa, lay.stemStart[a], lay.stemEnd[a], sb, lay.stemStart[b], lay.stemEnd[b]));
    }
    return lines;
  }
  if (kind === Kind.stub) {
    const s = lay.stubStem[ref];
    const partner = lay.stubPartner[ref];
    if (partner < 0) {
      // a folded stub: its partners listed, each with its links
      const folded = foldedOf(lay, ref);
      return [
        `${rangeText(lay, s)} → other chromosomes`,
        `Links: ${links(lay.stubCount[ref])} to ${folded.length} chromosomes, not shown`,
        `Partners: ${folded.map((f) => `${lay.names[f.partner]} (${Math.round(f.count).toLocaleString('en-US')})`).join(', ')}`,
      ];
    }
    return [`${rangeText(lay, s)} → ${lay.names[partner]}`, `Links: ${links(lay.stubCount[ref])}, not shown`];
  }
  if (kind === Kind.trumpet || kind === Kind.fan) {
    const partners = new Set<number>();
    for (let r = 0; r < lay.nr; r++) {
      if (lay.ribA[r] === ref) partners.add(lay.stemSystem[lay.ribB[r]]);
      else if (lay.ribB[r] === ref) partners.add(lay.stemSystem[lay.ribA[r]]);
    }
    for (let b = 0; b < lay.nb; b++)
      if (lay.stubStem[b] === ref && lay.stubPartner[b] >= 0) partners.add(lay.stubPartner[b]);
    for (let e = 0; e < lay.nf; e++) if (lay.stubStem[lay.foldOwner[e]] === ref) partners.add(lay.foldPartner[e]);
    const lines = [
      rangeText(lay, ref),
      `${lay.width === 'count' ? `Links: ${links(lay.stemCount[ref])}` : `Width: ${lay.width} ${valueText(lay.stemCount[ref])}`} to ${partners.size} chromosome${partners.size === 1 ? '' : 's'}`,
    ];
    if (lay.stemWithin[ref] > 0) lines.push(`Within: ${links(lay.stemWithin[ref])}`);
    if (lay.stemFan[ref] > 0)
      lines.push(
        `Fan: +${Math.round(lay.stemFan[ref]).toLocaleString('en-US')} weaker links to ${lay.stemFanTo[ref]} stems`,
      );
    lines.push(...perLayer(lay, lay.systems[lay.stemSystem[ref]], lay.stemStart[ref], lay.stemEnd[ref]));
    return lines;
  }
  return undefined;
}

/** The layer's ARIA description: its ribbons and chromosomes in view, the strongest pair; "no links" when none. */
export function summaryRondo(lay: RondoLay): string | undefined {
  if (lay.empty === 'layout') return undefined;
  if (lay.empty === 'pending') return 'links loading';
  if (lay.empty === 'failed') return 'links failed';
  let drawn = 0;
  for (let r = 0; r < lay.nr; r++) drawn += lay.ribDrawn[r];
  if (!drawn && !lay.nb) return 'no links';
  const blocks = new Set<number>();
  for (let s = 0; s < lay.ns; s++) if (lay.stemPlaced[s] && lay.stemCount[s] > 0) blocks.add(lay.stemSystem[s]);
  const r = lay.strongest;
  const strongest =
    r >= 0
      ? `; strongest ${lay.names[lay.stemSystem[lay.ribA[r]]]} ↔ ${lay.names[lay.stemSystem[lay.ribB[r]]]}, ${Math.round(lay.ribCount[r]).toLocaleString('en-US')}`
      : '';
  return `${drawn.toLocaleString('en-US')} ribbon${drawn === 1 ? '' : 's'}, ${blocks.size} chromosome${blocks.size === 1 ? '' : 's'}${strongest}`;
}

/** What the keyboard cursor adds at a position: the stem under it and its links. */
export function announceRondo(lay: RondoLay, system: string, position: number): string | undefined {
  const s = stemAt(lay, system, position);
  if (s < 0 || lay.stemCount[s] <= 0) return undefined;
  return `stem of ${links(lay.stemCount[s])}`;
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayerPlace, Legend } from '@molstar/seqx-core/extension';
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

/** A range as a span: `50–58 Mb`. */
export function spanText(start: number, end: number): string {
  const unit = end >= 1e6 ? ' Mb' : end >= 1e3 ? ' kb' : '';
  return `${at(start)}–${at(end)}${unit}`;
}

/** A stem's range on its system: `chr4:50–58 Mb`. */
function rangeText(lay: RondoLay, s: number): string {
  return `${lay.names[lay.stemSystem[s]]}:${spanText(lay.stemStart[s], lay.stemEnd[s])}`;
}

export const links = (n: number) => `${Math.round(n).toLocaleString('en-US')} link${Math.round(n) === 1 ? '' : 's'}`;

/** A field's value as the card says it: a number to 4 significant digits, a strand's sign as a minus. */
export const valueText = (v: string | number) => (typeof v === 'number' ? String(Number(v.toPrecision(4))) : v);

/**
 * The card's lines on kinds of links (_spec/README.md, Kinds of links, The card): a ribbon's links' shared fields and a
 * bin's or a stem's count on each layer, from `card.ts`, a chunk `rondoViewer` loads once composed (`install`); none
 * before it lands, and none headless (no card there).
 */
export const CardLines: {
  fields?: (lay: RondoLay, r: number) => string[];
  perLayer?: (lay: RondoLay, sysA: string, a0: number, a1: number, sysB?: string, b0?: number, b1?: number) => string[];
  /** The layer's legend (Colour, The legend). */
  legend?: (lay: RondoLay, options: { readonly title: string }) => readonly Legend[];
  /** A range's cytobands as the viewer whose frame holds the lay draws them (`q21.3`, `q21.1–q22.3`); none: none. */
  band?: (lay: RondoLay, system: string, start: number, end: number) => string | undefined;
} = {};

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
    // both ends, each its chromosome, band and span (novice 25: the card named one end)
    const end = (s: number) => {
      const band = CardLines.band?.(lay, lay.systems[lay.stemSystem[s]], lay.stemStart[s], lay.stemEnd[s]);
      return `${lay.names[lay.stemSystem[s]]}: ${band ? `${band} · ` : ''}${spanText(lay.stemStart[s], lay.stemEnd[s])}`;
    };
    const lines = [
      `${lay.names[lay.stemSystem[a]]} ↔ ${lay.names[lay.stemSystem[b]]}`,
      end(a),
      end(b),
      `Links: ${links(lay.ribLinks[ref])}${kindName ? ` (${kindName})` : ''}`,
    ];
    if (lay.width !== 'count') lines.push(`Width: ${lay.width} ${valueText(lay.ribCount[ref])}`);
    if (lay.ribStrand[ref] < 0) lines.push('Strand: reversed');
    // its links' fields (one link's own: its kind, strands, q, gene names); a bin: its count on each layer of the section
    if (lay.items && CardLines.fields) lines.push(...CardLines.fields(lay, ref));
    if (bins) {
      lines.push(`${bins} bins: zoom in to see them`);
      const sa = lay.systems[lay.stemSystem[a]];
      const sb = lay.systems[lay.stemSystem[b]];
      if (CardLines.perLayer)
        lines.push(
          ...CardLines.perLayer(lay, sa, lay.stemStart[a], lay.stemEnd[a], sb, lay.stemStart[b], lay.stemEnd[b]),
        );
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
        `Links: ${links(lay.stubLinks[ref])} to ${folded.length} chromosomes, not shown`,
        `Partners: ${folded.map((f) => `${lay.names[f.partner]} (${Math.round(f.count).toLocaleString('en-US')})`).join(', ')}`,
      ];
    }
    return [`${rangeText(lay, s)} → ${lay.names[partner]}`, `Links: ${links(lay.stubLinks[ref])}, not shown`];
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
    // the first row names what it counts: this stem's links, to its one partner by name (novice 21)
    const to = partners.size === 1 ? lay.names[partners.values().next().value!] : `${partners.size} chromosomes`;
    const lines = [rangeText(lay, ref), `This stem: ${links(lay.stemLinks[ref])} to ${to}`];
    if (lay.width !== 'count') lines.push(`Width: ${lay.width} ${valueText(lay.stemCount[ref])}`);
    if (lay.withinLinks[ref] > 0) lines.push(`Within: ${links(lay.withinLinks[ref])}`);
    if (lay.stemFan[ref] > 0)
      lines.push(
        `Fan: +${Math.round(lay.fanLinks[ref]).toLocaleString('en-US')} weaker links to ${lay.stemFanTo[ref]} stems`,
      );
    // each layer's links on the stem's chromosome, under a heading of their own (a row with no value)
    const all = CardLines.perLayer?.(lay, lay.systems[lay.stemSystem[ref]], 0, Number.POSITIVE_INFINITY) ?? [];
    if (all.length) lines.push(`${lay.names[lay.stemSystem[ref]]}, all links`, ...all);
    return lines;
  }
  return undefined;
}

/** The layer's ARIA description: its ribbons and chromosomes in view, the strongest pair; "no links" when none. */
export function summaryRondo(lay: RondoLay): string | undefined {
  if (lay.empty === 'layout') return undefined;
  if (lay.empty === 'pending') return 'links loading';
  if (lay.empty === 'shared') return 'links drawn by another layer';
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
  return `stem of ${links(lay.stemLinks[s])}`;
}

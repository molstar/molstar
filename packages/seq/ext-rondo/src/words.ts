/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayerPlace } from '@molstar/seqx-core/extension';
import { Kind, type RondoLay, stemAt } from './lay.js';

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
      `Links: ${links(lay.ribCount[ref])}${kindName ? ` (${kindName})` : ''}`,
    ];
    if (lay.ribStrand[ref] < 0) lines.push('Strand: reversed');
    if (bins) lines.push(`Bins: ${bins} bins, zoom in to see them`);
    return lines;
  }
  if (kind === Kind.stub) {
    const s = lay.stubStem[ref];
    return [
      `${rangeText(lay, s)} → ${lay.names[lay.stubPartner[ref]]}`,
      `Links: ${links(lay.stubCount[ref])}, not shown`,
    ];
  }
  if (kind === Kind.trumpet || kind === Kind.fan) {
    const partners = new Set<number>();
    for (let r = 0; r < lay.nr; r++) {
      if (lay.ribA[r] === ref) partners.add(lay.stemSystem[lay.ribB[r]]);
      else if (lay.ribB[r] === ref) partners.add(lay.stemSystem[lay.ribA[r]]);
    }
    for (let b = 0; b < lay.nb; b++) if (lay.stubStem[b] === ref) partners.add(lay.stubPartner[b]);
    const lines = [
      rangeText(lay, ref),
      `Links: ${links(lay.stemCount[ref])} to ${partners.size} chromosome${partners.size === 1 ? '' : 's'}`,
    ];
    if (lay.stemWithin[ref] > 0) lines.push(`Within: ${links(lay.stemWithin[ref])}`);
    if (lay.stemFan[ref] > 0)
      lines.push(
        `Fan: +${Math.round(lay.stemFan[ref]).toLocaleString('en-US')} weaker links to ${lay.stemFanTo[ref]} stems`,
      );
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

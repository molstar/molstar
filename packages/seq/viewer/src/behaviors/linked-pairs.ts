/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseAddress, printAddress, type Segment } from '@molstar/core/repr/tree';
import {
  type ElementKey,
  type Elements,
  EmptyElements,
  EmptyLocation,
  type Location,
  unionElements,
} from '@molstar/seqx-core/query';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import type { Emphasized } from '@molstar/seqx-viewer/models';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';

// The `linkedPairs` behavior (_spec/README.md "Behaviors", Linked pairs; the design's viewer-models.md, Behaviors): on
// a hovered or selected range element (a feature, a CDR), the pairs with an end in its extent, on every pairs layer of
// the view shown, written to a highlight channel of its own (`linked-pairs`: the highlight pass, in the channel's
// color, never the selection's: they are related to the selection, not in it), cleared when the hover or selection
// goes. Built on the public API alone: the query language (`<pairs>/(items :overlapping (extent :of <element>))`,
// pairs answering `items`), the highlight channel (its element keys, `EmphasisHighlight.elements`), the hover and
// selection models. The pairs layer's own rule (a hovered position lights its partners) is the layer's, unchanged.
// The hover's query is cached per hovered element key (the view and the model it was made on): a hover allocates
// nothing a frame past its first.

/** What lights the pairs: the hover, the selection. */
export type LinkedPairsSource = 'hover' | 'selection';

export interface LinkedPairsParams {
  /** What lights the pairs (both by default). */
  readonly on?: readonly LinkedPairsSource[];
  /** The pairs layers it reads, by address (`/views/v/s/contacts/arcs`); every pairs layer of the view when absent. */
  readonly layers?: readonly string[];
}

export interface LinkedPairsBehavior {
  /** The highlight channel it writes (`linked-pairs`). */
  readonly channel: string;
  /** The pairs lit now (the channel's keys): the hover's and the selection's. */
  lit(): Elements;
  dispose(): void;
}

/** The channel the behavior writes: its color is `highlight.colors['linked-pairs']`, else `highlight.color`. */
export const LinkedPairsChannel = 'linked-pairs';

const BothSources: readonly LinkedPairsSource[] = ['hover', 'selection'];

/** A pairs layer's annotation: its data's address up to its first form (`<pairs>/(except …)` reads `<pairs>`). */
function annotationOf(layer: ResolvedLayer): string | undefined {
  const address = layer.data?.address;
  const segments = address ? parseAddress(address).pattern?.segments : undefined;
  if (!segments) return undefined;
  const ids: string[] = [];
  for (const s of segments) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  return ids.length ? printAddress(ids) : undefined;
}

/** The pairs annotations the view's pairs layers draw (those `layers` names, when given), each once. */
function pairsOf(view: ResolvedView, layers: readonly string[] | undefined): string[] {
  const out = new Set<string>();
  const add = (l: ResolvedLayer) => {
    if (l.type !== 'pairs' || (layers && (l.address === undefined || !layers.includes(l.address)))) return;
    const a = annotationOf(l);
    if (a !== undefined) out.add(a);
  };
  for (const s of view.sections) {
    for (const t of s.tracks) for (const l of t.layers) add(l);
    for (const b of s.blocks ?? []) for (const ls of b.layers) for (const l of ls) add(l);
  }
  return [...out];
}

/**
 * Where an element is, as a location form: its extent's ranges as index ranges on their systems (a hit's extent: what a
 * layer read, a source's item too), else the extent of its item form (`(extent :of <owner>/(item :id <key>))`).
 */
function extentForm(k: ElementKey, extent: Location | undefined): string | undefined {
  if (extent?.parts.length) {
    const ranges: string[] = [];
    // a zero-length range (a boundary) has no index range: left to the item form, whose extent the forms read
    for (const p of extent.parts)
      for (let i = 0; i + 1 < p.ranges.length; i += 2)
        if (p.ranges[i + 1] > p.ranges[i])
          ranges.push(`${p.system}/(range :from ${p.ranges[i]} :to ${p.ranges[i + 1] - 1} :numbering index)`);
    if (ranges.length) return ranges.length === 1 ? ranges[0] : `(merge :of [${ranges.join(' ')}])`;
  }
  const owner = parseAddress(k.owner).pattern?.segments;
  if (!owner) return undefined;
  const item: Segment = { kind: 'item', params: { id: k.key } };
  return `(extent :of ${printAddress([...owner, item])})`;
}

/** Whether a location has a range of at least a position (a boundary's has none). */
function hasRanges(l: Location): boolean {
  for (const p of l.parts)
    for (let i = 0; i + 1 < p.ranges.length; i += 2) if (p.ranges[i + 1] > p.ranges[i]) return true;
  return false;
}

export function linkedPairsBehavior(ctx: BehaviorContext, params: LinkedPairsParams = {}): LinkedPairsBehavior {
  const { view, hover, selection, highlight, addresses, spec } = ctx;
  const on = params.on ?? BothSources;
  const onHover = on.includes('hover');
  const onSelection = on.includes('selection');

  // what the caches were made on: the view shown and its model
  let madeFor: [unknown, unknown] = [undefined, undefined];
  let pairs: readonly string[] = [];
  const owners = new Set<string>();
  /** The pairs lit by an element, by its owner and key. */
  const byElement = new Map<string, Map<string | number, Elements>>();
  /** The channel's value for each set of keys written (by identity: a hover back on an element writes the same). */
  let values = new WeakMap<Elements, Emphasized>();

  const fresh = () => {
    const v = view.resolved.value;
    const m = spec.model.value;
    if (madeFor[0] === v && madeFor[1] === m) return;
    madeFor = [v, m];
    pairs = v ? pairsOf(v, params.layers) : [];
    owners.clear();
    for (const p of pairs) owners.add(p);
    byElement.clear();
    values = new WeakMap();
    selectionFor = [undefined, EmptyElements];
  };

  /** The pairs with an end in an element's extent, on every pairs annotation drawn; cached by its owner and key. */
  const litBy = (k: ElementKey, extent?: Location): { lit: Elements; waiting: boolean } => {
    let byKey = byElement.get(k.owner);
    const kept = byKey?.get(k.key);
    if (kept) return { lit: kept, waiting: false };
    const r = query(extentForm(k, extent));
    // what waits for its data (an item read from a source, named by its address) is asked again next time
    if (r.waiting) return r;
    if (!byKey) byElement.set(k.owner, (byKey = new Map()));
    byKey.set(k.key, r.lit);
    return r;
  };

  /** The pairs with an end in a location form, on every pairs annotation drawn, and whether a query waits. */
  const query = (where: string | undefined): { lit: Elements; waiting: boolean } => {
    let lit = EmptyElements;
    let waiting = false;
    if (where !== undefined)
      for (const p of pairs) {
        const r = addresses.resolve(`${p}/(items :overlapping ${where})`);
        if (r.result?.type === 'keys') lit = unionElements(lit, r.result.keys);
        else if (r.status === 'waiting') waiting = true;
      }
    return { lit, waiting };
  };

  /** The selection's last query: what it was made of (its value) and what it lit. */
  let selectionFor: [unknown, Elements] = [undefined, EmptyElements];

  let hoverLit = EmptyElements;
  let selectionLit = EmptyElements;
  /** The union last made, and of what. */
  let both: [Elements, Elements, Elements] = [EmptyElements, EmptyElements, EmptyElements];

  const write = () => {
    let keys: Elements;
    if (!hoverLit.keys.length) keys = selectionLit;
    else if (!selectionLit.keys.length) keys = hoverLit;
    else {
      if (both[0] !== hoverLit || both[1] !== selectionLit)
        both = [hoverLit, selectionLit, unionElements(hoverLit, selectionLit)];
      keys = both[2];
    }
    if (!keys.keys.length) {
      highlight.clear(LinkedPairsChannel);
      return;
    }
    let value = values.get(keys);
    if (!value) values.set(keys, (value = { elements: keys, location: EmptyLocation, shown: EmptyLocation }));
    highlight.set(LinkedPairsChannel, value);
  };

  const readHover = () => {
    if (!onHover) return;
    fresh();
    const r = hover.current;
    const el = r.owner === undefined ? undefined : r.element;
    // a pair hovered is the pairs layer's own (its arc lights); a position is the layer's rule (its partners)
    hoverLit = el && pairs.length && !owners.has(el.owner) ? litBy(el, r.extent).lit : EmptyElements;
  };

  const readSelection = () => {
    if (!onSelection) return;
    fresh();
    const sel = selection.value.value;
    const keys = sel.elements.keys;
    let out = EmptyElements;
    if (!pairs.length || !keys.some((k) => !owners.has(k.owner))) out = EmptyElements;
    else if (selectionFor[0] === sel) out = selectionFor[1];
    else {
      let waiting = false;
      if (!sel.location.parts.length && hasRanges(sel.shown) && keys.every((k) => !owners.has(k.owner))) {
        // range elements alone, drawn: what is drawn of them is their extents (a source's items too), one query
        const r = query(extentForm(keys[0], sel.shown));
        out = r.lit;
        waiting = r.waiting;
      }
      // elements beside positions or pairs, or none of their extents known: each range element by its own
      else
        for (const k of keys) {
          if (owners.has(k.owner)) continue;
          const r = litBy(k);
          out = unionElements(out, r.lit);
          waiting ||= r.waiting;
        }
      // what waited is asked again on the next change
      selectionFor = waiting ? [undefined, EmptyElements] : [sel, out];
    }
    selectionLit = out;
  };

  const stops: (() => void)[] = [];
  if (onHover)
    stops.push(
      hover.version.listen(() => {
        readHover();
        write();
      }),
    );
  if (onSelection)
    stops.push(
      selection.value.listen(() => {
        readSelection();
        write();
      }),
    );
  // a new view or model: the pairs drawn and the queries made again
  const again = () => {
    readHover();
    readSelection();
    write();
  };
  stops.push(view.resolved.listen(again), spec.model.listen(again));
  again();

  return {
    channel: LinkedPairsChannel,
    lit: () => highlight.get(LinkedPairsChannel)?.elements ?? EmptyElements,
    dispose() {
      for (const s of stops) s();
      stops.length = 0;
      highlight.clear(LinkedPairsChannel);
    },
  };
}

/** The `linkedPairs` behavior's provider: listed in `DefaultSeqViewerSpec()`, removed to opt out. */
export const LinkedPairsProvider: BehaviorProvider<LinkedPairsParams | undefined, LinkedPairsBehavior> = {
  name: 'linkedPairs',
  create: (ctx, params) => linkedPairsBehavior(ctx, params),
};

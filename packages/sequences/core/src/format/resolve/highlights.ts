/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's `highlights` (seqspec:format#sections): its location sets shown persistently, each entry its data's items
// as item loci on the set, for a viewer's highlight channels (the Seq* viewer's `set:<section>/<n>`) and Nightingale's
// `highlight`. A query node is an annotation, so it is data here as any other. Edge rate: allocates.

import { Annotations, type ItemLoci, KeyKind, type KeyTarget, OrderedSet, type SeqSet } from "@molstar/seq-core/model";
import type { ResolvedView } from "./view.js";

export type SectionHighlightStyle = "band" | "highlight" | "accent";
export interface SectionHighlight {
  /** The entry's `data`, `items` and `style` as written (defaults applied). */
  readonly data: string;
  readonly items: readonly string[] | null;
  readonly style: SectionHighlightStyle;
  /** Its items on the set; null: its data is not in the set (an unavailable computed node) or names none of `items`. */
  readonly loci: ItemLoci | null;
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };

/** Section `s`'s `highlights` (its index in the view), in order. */
export function sectionHighlights(view: ResolvedView, set: SeqSet, s: number): SectionHighlight[] {
  const list = view.sections.params[s]?.highlights;
  if (!Array.isArray(list)) return [];
  return list.map((e: { data: string; items?: string[] | null; style?: SectionHighlightStyle }) => {
    const items = e.items ?? null;
    return { data: e.data, items, style: e.style ?? "band", loci: itemLoci(set, e.data, items) };
  });
}

function itemLoci(set: SeqSet, ref: string, ids: readonly string[] | null): ItemLoci | null {
  if (!set.byName(ref, t) || t.kind !== KeyKind.annotation) return null;
  const ann = set.annotations[t.id];
  if (!ann || ann.itemCount === 0) return null;
  let items = OrderedSet.ofRange(0, ann.itemCount);
  if (ids) {
    const picked: number[] = [];
    for (const id of ids) {
      const key = set.strings.find(id);
      const item = key >= 0 ? Annotations.itemOf(ann, key as never) : -1;
      if (item >= 0 && !picked.includes(item)) picked.push(item);
    }
    if (picked.length === 0) return null;
    items = OrderedSet.ofSortedArray(Int32Array.from(picked).sort());
  }
  return { kind: "item-loci", set, elements: [{ annotation: t.id as never, items: items as never }] };
}

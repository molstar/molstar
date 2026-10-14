/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vocabulary (seqspec:libraries#interactions): the kinds of an interaction, a bit per kind (`KindMask`), and how
// each is labeled, colored and drawn. One set of opinions; another vocabulary is another extension.

/** The kinds, in bit order. `other`: any category this vocabulary does not know. */
export const Kinds = [
  "contact",
  "hbond",
  "salt_bridge",
  "pi_stacking",
  "cation_pi",
  "hydrophobic",
  "halogen",
  "metal",
  "water_bridge",
  "covalent",
  "other",
] as const;
export type KindName = (typeof Kinds)[number];

/** A kind's index (its bit is `1 << index`). */
export const Kind = Object.fromEntries(Kinds.map((k, i) => [k, i])) as { readonly [K in KindName]: number };

/** Label, color and glyph shape per kind: written out into the kernel fallbacks, so any reader draws the same marks. */
export const KindStyle: { readonly [K in KindName]: readonly [label: string, color: string, shape: string] } = {
  contact: ["contact", "#bab0ac", "pin"],
  hbond: ["hydrogen bond", "#4e79a7", "diamond"],
  salt_bridge: ["salt bridge", "#e15759", "cross"],
  pi_stacking: ["π stacking", "#59a14f", "square"],
  cation_pi: ["cation–π", "#f28e2b", "star"],
  hydrophobic: ["hydrophobic contact", "#edc948", "circle"],
  halogen: ["halogen bond", "#76b7b2", "triangle"],
  metal: ["metal coordination", "#b07aa1", "star"],
  water_bridge: ["water-mediated", "#9c9ede", "triangle"],
  covalent: ["covalent bond", "#9c755f", "square"],
  other: ["other", "#79706e", "circle"],
};

const byName = new Map<string, number>(Kinds.map((k, i) => [k, i]));

/** One bit per kind: this extension's only. */
export type KindMask = number;
export const KindMask = {
  all: (1 << Kinds.length) - 1,
  of: (...kinds: readonly KindName[]): KindMask => kinds.reduce((m, k) => m | (1 << Kind[k]), 0),
  has: (mask: KindMask, kind: number) => (mask & (1 << kind)) !== 0,
  /** The number of kinds set. */
  count(mask: KindMask): number {
    let n = 0;
    for (let m = mask; m; m &= m - 1) n++;
    return n;
  },
  names: (mask: KindMask): KindName[] => Kinds.filter((_, i) => mask & (1 << i)),
};

/** A category's kind: an unknown or a missing one is `other`; `contact` is a stated kind. */
export const kindByName = (category: string | null): number =>
  category == null ? Kind.other : (byName.get(category) ?? Kind.other);

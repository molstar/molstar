/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's color schemes: `b.colors`, the kernel's schemes as typed `scheme` theme refs (a value scheme takes its
// domain, field and direction), a registered scheme by name, a custom one, one derived from another, a sensible default
// per subject, and the schemes for a subject. No runtime dependency (seqspec:format#builders): the table of names and
// what each fits is a twin of the format's `KernelColorSchemes`, kept equal by its schemes test.

import type { Alphabet, ColorStop, LetterColors } from "@molstar/seq-core/format";

/**
 * An inline scheme (seqspec:format#themes, `scheme`'s `custom`); `dark` absent: derived (seqspec:format#color-schemes).
 */
export type CustomScheme =
  | (LetterColors & {
      readonly kind: "letters";
      readonly alphabets: readonly Alphabet[];
      readonly dark?: "same" | LetterColors;
    })
  | {
      readonly kind: "sequential" | "diverging";
      readonly colors: readonly ColorStop[];
      readonly dark?: "same" | readonly ColorStop[];
    }
  | {
      readonly kind: "categorical";
      readonly colors: readonly string[] | Readonly<Record<string, string>>;
      readonly missing?: string | null;
      readonly dark?: "same" | readonly string[] | Readonly<Record<string, string>>;
    };
export type SchemeOptions = {
  /** A letters scheme's letters recolored. */
  readonly overrides?: Readonly<Record<string, string>>;
  /** Values and categories: the field colored (`field_palette`'s). */
  readonly field?: string | null;
  /** Values: the ramp's ends; null: the layer's. */
  readonly domain?: readonly [number | null, number | null] | null;
  readonly reverse?: boolean;
};
/** A `scheme` color theme reference (a type alias, so it fits the format's `ThemeRef`). */
export type SchemeRef = {
  readonly theme: "scheme";
  readonly params: SchemeOptions & { readonly name?: string; readonly custom?: CustomScheme };
};

type Fits = { readonly letters: readonly Alphabet[]; readonly columns?: true } | "values" | "categories";
const P: Fits = { letters: ["protein"] };
const NT: Fits = { letters: ["dna", "rna"] };
/** The kernel's schemes and what each fits. Twin: `KernelColorSchemes` in the format's `palettes/schemes.ts`. */
export const KernelSchemeFits = {
  clustal_x: P,
  zappo: P,
  taylor: P,
  hydrophobicity: P,
  percent_identity: { letters: ["protein", "dna", "rna", "generic"], columns: true },
  nucleotide: NT,
  purine_pyrimidine: NT,
  nucleotide_cb: NT,
  viridis: "values",
  cividis: "values",
  magma: "values",
  blues: "values",
  rd_bu: "values",
  category10: "categories",
  okabe_ito: "categories",
  tableau10: "categories",
} as const satisfies Readonly<Record<string, Fits>>;
export type KernelSchemeName = keyof typeof KernelSchemeFits;

/** What a layer's color is about: letters of an alphabet (on an alignment: column tables too), values, categories. */
export type SchemeFor =
  | { readonly alphabet: Alphabet; readonly alignment?: boolean }
  | { readonly values: true | "sequential" | "diverging" }
  | { readonly categories: true };

const ref = (params: SchemeRef["params"]): SchemeRef => ({ theme: "scheme", params });
const named = (name: KernelSchemeName) => ref({ name });
type ValueOptions = Pick<SchemeOptions, "field" | "domain" | "reverse">;
const ramp =
  (name: KernelSchemeName) =>
  (o: ValueOptions = {}): SchemeRef =>
    ref({ name, ...o });

export const Colors = {
  clustalX: named("clustal_x"),
  zappo: named("zappo"),
  taylor: named("taylor"),
  hydrophobicity: named("hydrophobicity"),
  percentIdentity: named("percent_identity"),
  nucleotide: named("nucleotide"),
  purinePyrimidine: named("purine_pyrimidine"),
  nucleotideCb: named("nucleotide_cb"),
  viridis: ramp("viridis"),
  cividis: ramp("cividis"),
  magma: ramp("magma"),
  blues: ramp("blues"),
  rdBu: ramp("rd_bu"),
  category10: named("category10"),
  okabeIto: named("okabe_ito"),
  tableau10: named("tableau10"),
  /** A registered scheme by name (an extension's `<ns>:<name>`). */
  scheme: (name: string, o: SchemeOptions = {}): SchemeRef => ref({ name, ...o }),
  /** An inline scheme. */
  custom: (scheme: CustomScheme, o: SchemeOptions = {}): SchemeRef => ref({ custom: scheme, ...o }),
  /** `base` with options over its own; `overrides` merged, the given letters last (none: none written). */
  from: (base: SchemeRef, o: SchemeOptions): SchemeRef => {
    const merged = { ...base.params, ...o };
    if (base.params.overrides || o.overrides) merged.overrides = { ...base.params.overrides, ...o.overrides };
    return ref(merged);
  },
  /**
   * A sensible default: Clustal X on a protein alignment, Zappo on a protein, Nucleotide on DNA and RNA, percent
   * identity on another alphabet's alignment, Viridis on values (RdBu when diverging), Okabe–Ito on categories
   * (colorblind-safe); null: none fits (a generic alphabet's lone sequence).
   */
  auto(of: SchemeFor): SchemeRef | null {
    if ("values" in of) return of.values === "diverging" ? Colors.rdBu() : Colors.viridis();
    if ("categories" in of) return Colors.okabeIto;
    if (of.alphabet === "protein") return of.alignment ? Colors.clustalX : Colors.zappo;
    if (of.alphabet !== "generic") return Colors.nucleotide;
    return of.alignment ? Colors.percentIdentity : null;
  },
  /** The kernel's schemes that fit, in the picker's order. */
  list(of: SchemeFor): KernelSchemeName[] {
    return (Object.entries(KernelSchemeFits) as [KernelSchemeName, Fits][]).flatMap(([name, fits]) => {
      if (fits === "values" || fits === "categories") return fits in of ? [name] : [];
      if (!("alphabet" in of) || !fits.letters.includes(of.alphabet)) return [];
      return !("columns" in fits) || of.alignment ? [name] : [];
    });
  },
};

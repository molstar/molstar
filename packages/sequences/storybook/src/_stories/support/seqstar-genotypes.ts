/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 17's choices (the genotype matrix, `core/_test/format/cases.ts`, in the Viewer cases story), shared with the smoke
// test that applies each: the matrix layer's `order` and `filter` over its sample table, and its `density`, as patches;
// generic params only (there is no genotypes preset).

import { GENOTYPE_SAMPLES, GENOTYPE_SITES, genotypeSamples } from "../../../../core/_test/format/cases.ts";

export { GENOTYPE_SAMPLES, GENOTYPE_SITES };
/** The matrix layer a mode's patch sets. */
export const GENOTYPES_LAYER = "genotypes.matrix";

export interface GenotypeMode {
  readonly key: string;
  readonly label: string;
  readonly order: unknown;
  readonly filter: readonly unknown[];
}
const BY_POPULATION = { data: "samples", field: "population" };
const BY_ALLELES = { data: "samples", field: "alt_alleles", direction: "desc" };
/** The layer's `order` and `filter` for each choice. */
export const GENOTYPE_MODES: readonly GenotypeMode[] = [
  { key: "data", label: "Data order", order: "data", filter: [] },
  { key: "population", label: "By population", order: BY_POPULATION, filter: [] },
  { key: "alleles", label: "By population, then alternative alleles", order: [BY_POPULATION, BY_ALLELES], filter: [] },
  {
    key: "eas",
    label: "EAS and EUR only",
    order: [BY_POPULATION, BY_ALLELES],
    filter: [{ data: "samples", field: "population", values: ["EAS", "EUR"] }],
  },
  {
    key: "carriers",
    label: "Samples with the most alternative alleles (cells)",
    order: [BY_POPULATION, BY_ALLELES],
    filter: [{ data: "samples", field: "alt_alleles", min: carrierMin() }],
  },
];

/** The fewest alternative alleles that keeps about 24 samples: the matrix then draws one band each, as cells. */
function carrierMin(): number {
  const sorted = [...genotypeSamples().alleles].sort((a, b) => b - a);
  return sorted[23];
}

/** How a device row of many samples picks its color: the layer's `density`. */
export const GENOTYPE_DENSITIES = [
  { key: "majority", label: "The most common color" },
  { key: "maximum", label: "The strongest genotype" },
] as const;
/** The patch that sets the matrix layer's `density`. */
export const genotypeDensityPatch = (density: string) => ({
  format: "seq_view_spec_patch" as const,
  ops: [{ op: "set" as const, ref: GENOTYPES_LAYER, params: { density } }],
});

/** The patch that sets a mode on the matrix layer. */
export const genotypeModePatch = (m: GenotypeMode) => ({
  format: "seq_view_spec_patch" as const,
  ops: [{ op: "set" as const, ref: GENOTYPES_LAYER, params: { order: m.order, filter: m.filter } }],
});

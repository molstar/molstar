/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the record script keeps of the two responses it fetches (`scripts/record.ts`): members and records removed,
// none renamed or moved, so a spec's paths read the live response too. Kept here, in the package's sources, so a unit
// test holds the rules; the script imports them.

/** The licences a fixture may carry (examples.md, "Licences allowed"); the script refuses any other. */
export const AllowedLicences = [
  'CC0-1.0',
  'CC-BY-4.0',
  'LicenseRef-NCBI-public-data',
  'LicenseRef-seqx-generated',
] as const;

type Json = null | boolean | number | string | Json[] | { [key: string]: Json };
export type Obj = { [key: string]: Json };

/** The entry's members the views read; every other (comments, references, cross references…) is removed. */
export const EntryMembers = ['entryType', 'primaryAccession', 'uniProtkbId', 'organism', 'sequence', 'features'];
/** A feature's members that are evidence, not what a view reads. */
export const FeatureDropped = ['evidences', 'featureCrossReferences', 'ligand'];
/** Natural variants are the variants fixture's (a thousand rows); every other feature type stays. */
export const EntryDroppedTypes = ['Natural variant'];

export function trimEntry(entry: Obj): Obj {
  const out: Obj = {};
  for (const k of EntryMembers) if (k in entry) out[k] = entry[k];
  out.features = (entry.features as Obj[])
    .filter((f) => !EntryDroppedTypes.includes(f.type as string))
    .map((f) => Object.fromEntries(Object.entries(f).filter(([k]) => !FeatureDropped.includes(k))));
  return out;
}

/**
 * The variants kept are UniProt's own (`sourceType` `uniprot`, or `mixed`: UniProt's, large-scale studies report them
 * too; CC-BY-4.0), with their UniProt and ClinVar xrefs by name; a record only large-scale studies report is removed,
 * its ClinVar xref or not (no licence is recorded for it).
 */
export const KeptSources = ['UniProt', 'ClinVar'];
export const VariantMembers = [
  'begin',
  'end',
  'xrefs',
  'consequenceType',
  'wildType',
  'mutatedType',
  'sourceType',
  'clinicalSignificances',
];
export const ResponseMembers = [
  'accession',
  'entryName',
  'proteinName',
  'geneName',
  'organismName',
  'taxid',
  'sequence',
  'features',
];

export function trimVariation(response: Obj): Obj {
  const out: Obj = {};
  for (const k of ResponseMembers) if (k in response) out[k] = response[k];
  const features: Obj[] = [];
  // the canonical sequence's positions: the response also holds the changes of other transcripts' proteins
  const length = (response.sequence as string).length;
  for (const f of response.features as Obj[]) {
    const names = ((f.xrefs ?? []) as Obj[]).map((x) => x.name as string).filter((n) => KeptSources.includes(n));
    // the views read single-residue changes with a source of known licence
    const at = Number(f.begin);
    if (
      f.sourceType === 'large_scale_study' ||
      !names.length ||
      f.begin !== f.end ||
      !f.wildType ||
      !f.mutatedType ||
      at < 1 ||
      at > length
    )
      continue;
    const kept: Obj = {};
    for (const k of VariantMembers) if (k in f) kept[k] = f[k];
    // one xref a source, by its name; a significance by its type
    kept.xrefs = [...new Set(names)].map((name) => ({ name }));
    if (Array.isArray(f.clinicalSignificances))
      kept.clinicalSignificances = (f.clinicalSignificances as Obj[]).map((c) => ({ type: c.type as string }));
    features.push(kept);
  }
  out.features = features;
  return out;
}

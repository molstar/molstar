/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { Diff } from './baseline.ts';
import type { LinkKind } from './candidates.ts';

/** A run's result per component, in `<tier>.json` (plan §4, §5.1). */
export type Bucket = 'pass' | 'fail' | 'unsanitizable' | 'known' | 'excluded' | 'error';
/** L0 bond by bond, L1 canonical SMILES, L2 InChI, L3 InChI with keto-enol and 1,5 tautomerism (plan §5) */
export type Level = 'L0' | 'L1' | 'L2' | 'L3';
export interface Row {
  /** the component; in a dataset, `<entry>_<component>` (plan §6) */
  id: string;
  /** a dataset's row: its component */
  comp?: string;
  /** a dataset's pair judged at another tier than the run's, from the dataset file (plan §6: tier B) */
  judgedAt?: string;
  bucket: Bucket;
  level?: Level;
  reason?: string;
  /** leaving atoms without coordinates, compared without them (tier B) */
  leavingMissing?: number;
  expected?: string;
  perceived?: string;
  l0?: string[];
  /** the L0 bonds by their atoms' positions in the SDFs (0-based), for the report's highlights */
  l0Atoms?: [number, number][];
  diagnostics?: string[];
  /** tier C: the entry, the copy compared, its link kind and number of caps (plan §3.2, §4.2) */
  entry?: string;
  copy?: string;
  kind?: LinkKind;
  /** tier C: a copy covalently linked to another residue was found in the entry (excluded or not) */
  linked?: boolean;
  caps?: number;
}

/** `<tier>.json` */
export interface Run {
  /** `ccd`, the whole CCD, or a dataset's name (plan §6) */
  dataset: string;
  /** `all`, `filter` or `sample`: a baseline is made from `all` only */
  selection: string;
  tier: string;
  connectivity: string;
  hydrogens: string;
  components: number;
  counts: Record<string, number>;
  levels: Record<string, number>;
  byKind?: Record<string, Record<string, number>>;
  /** tier C: the non-polymer components linked in their entry, every bucket (the reports' covalent-ligand flag) */
  covalentLigands?: string[];
  /** against the baseline, when there is one (plan §6) */
  diff?: Diff;
  rows: Row[];
}

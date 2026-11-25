/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { initRDKit, isRDKitParseError, type Mol } from 'rdkit-esm';

export type Identifiers = { inchi: string; smiles: string } | { error: string };

export type RDKit = Awaited<ReturnType<typeof initRDKit>>;
let rdkit: Promise<RDKit> | undefined;

/** RDKit, loaded once (compare and report share it). */
export function loadRDKit(): Promise<RDKit> {
  if (!rdkit) {
    // rdkit-esm 0.8.0's Node path reads the free variables `require` and `__dirname`, absent from an ES module
    const g = globalThis as { require?: NodeJS.Require; __dirname?: string };
    g.require ??= createRequire(import.meta.url);
    g.__dirname ??= path.dirname(fileURLToPath(import.meta.resolve('rdkit-esm')));
    rdkit = initRDKit();
  }
  return rdkit;
}

/**
 * A mol block's standard InChI and canonical SMILES, or RDKit's message when it cannot sanitize it (never autofixed,
 * plan §4.1).
 */
export async function createIdentifiers() {
  const rdkit = await loadRDKit();
  return (molBlock: string): Identifiers => {
    try {
      const mol = rdkit.getMol(molBlock, { throwOnError: true });
      const ids = { inchi: mol.getInchi(), smiles: mol.getSmiles() };
      mol.delete();
      return ids;
    } catch (e) {
      if (isRDKitParseError(e)) return { error: e.rdkitMessage };
      throw e;
    }
  };
}

/**
 * InChI's keto-enol (`KET`) and 1,5-shift (`15T`) tautomerism options, beyond standard InChI's mobile H: a perceived
 * graph that differs from the reference only so passes at L3 (plan §5). Asked only after the standard InChIs differ.
 */
const Tautomers: NonNullable<Parameters<Mol['getInchi']>[0]> = ['KET', '15T'];

/** A mol block's InChI with the tautomerism options (it sanitized once already, as `createIdentifiers` checked). */
export async function createTautomerInchi() {
  const rdkit = await loadRDKit();
  return (molBlock: string) => {
    const mol = rdkit.getMol(molBlock, { throwOnError: true });
    const inchi = mol.getInchi(Tautomers);
    mol.delete();
    return inchi;
  };
}

/**
 * The layers compared (plan §5): formula, connections and H. Stereo (`/b`, `/t`, `/m`, `/s`) is left out: the CCD's
 * does not always agree with its own coordinates, and an order error shows in the formula or the `/h` layer anyway.
 * Non-standard layers (`/f` fixed H, `/i` isotopes, `/r` reconnected) end the comparison. `/q` and `/p` are dropped
 * unless `charges` (tier 0): without H protonation and charge are unknowable, and with H they are not what perception
 * is measured on.
 */
export function comparableInchi(inchi: string, charges: boolean) {
  const [, ...layers] = inchi.split('/');
  const kept: string[] = [];
  for (const layer of layers) {
    const kind = layer.charAt(0);
    if (kind === 'f' || kind === 'i' || kind === 'r') break;
    if ('btms'.includes(kind) || (!charges && (kind === 'q' || kind === 'p'))) continue;
    kept.push(layer);
  }
  return kept.join('/');
}

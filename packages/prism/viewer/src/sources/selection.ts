/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { structureElementStatsLabel } from '@molstar/model/model/label';
import { Loci } from '@molstar/model/model/loci';
import { QueryContext, StructureElement, StructureSelection } from '@molstar/model/model/structure';
import { isSingleLigandLoci } from '@molstar/model/model/structure/structure/util/superposition-ligand';
import { polymer } from '@molstar/plugin/state/queries/structure/type';
import { summarize } from '../targets/index.js';

// The plugin's selection as the selection bar says it (contexts.md §7: the plugin is read here, in `sources/`): what is
// selected, in Mol*'s own words and as granular as the selection is (a residue's label, a chain's, "12 Residues"),
// and nothing when nothing is.

/**
 * The summary of selections (a structure's each): Mol*'s label of their statistics, condensed, without markup;
 * undefined for none. One thing is named (a residue, an atom, a chain, a ligand), more are counted by kind ("12 Residues").
 * Several structures' selections are counted too (Mol*'s statistics keep the names of a single one only).
 */
export function selectionSummary(selections: readonly StructureElement.Loci[]): string | undefined {
  const held = selections.filter((l) => !Loci.isEmpty(l));
  if (held.length === 0) return undefined;
  let stats = StructureElement.Stats.ofLoci(held[0]);
  for (let i = 1; i < held.length; i++) {
    const next = StructureElement.Stats.create();
    StructureElement.Stats.add(next, stats, StructureElement.Stats.ofLoci(held[i]));
    stats = next;
  }
  // One thing (an atom, a residue, a chain, a structure) is named; more are counted by kind.
  const things =
    stats.elementCount + stats.conformationCount + stats.residueCount + stats.chainCount + stats.structureCount;
  const named = held.length === 1 && things === 1;
  const label = structureElementStatsLabel(stats, { countsOnly: !named, condensed: true, htmlStyling: false }).trim();
  // Mol* capitalizes its kinds in a count ("12 Residues"); in a sentence of the bar they read in lower case.
  return named
    ? label
    : label.replace(
        /(\d) (Residues?|Chains?|Atoms?|Elements?|Structures?|Conformations?|Units?)/g,
        (_, n, kind) => `${n} ${kind.toLowerCase()}`,
      );
}

// --- what a selection can be superposed by -----------------------------------------------------------------------

/** The polymer part of a selection: its atoms in polymer entities (waters, ligands and ions left out). */
export function polymerPart(loci: StructureElement.Loci): StructureElement.Loci {
  const sub = StructureElement.Loci.toStructure(loci);
  const part = StructureSelection.toLociWithSourceUnits(polymer.query(new QueryContext(sub)));
  return StructureElement.Loci.remap(part, loci.structure);
}

/** What one structure's selection holds, for the Superpose section to say what a method would make of it. */
export interface PickedStructure {
  /** The structure's state ref. */
  readonly ref: string;
  readonly label: string;
  /** "42 residues (350 atoms) in chains A and B". */
  readonly line: string;
  readonly atoms: number;
  readonly residues: number;
  /** The author chains it is in, in structure order. */
  readonly chains: readonly string[];
  /** Atoms in polymer entities, and the chains (units) they are in. */
  readonly polymerAtoms: number;
  readonly polymerChains: number;
  /** The selection is one ligand residue. */
  readonly ligand: boolean;
}

export function pickedOf(ref: string, label: string, loci: StructureElement.Loci): PickedStructure {
  const part = polymerPart(loci);
  const selected = summarize(loci);
  return {
    ref,
    label,
    line: selected.line,
    atoms: selected.atoms,
    residues: selected.residues,
    chains: selected.chains,
    polymerAtoms: StructureElement.Loci.size(part),
    polymerChains: part.elements.length,
    ligand: !!isSingleLigandLoci(loci),
  };
}

export const samePicked = (a: PickedStructure, b: PickedStructure) =>
  a.ref === b.ref &&
  a.label === b.label &&
  a.line === b.line &&
  a.atoms === b.atoms &&
  a.residues === b.residues &&
  a.chains.length === b.chains.length &&
  a.chains.every((c, i) => c === b.chains[i]) &&
  a.polymerAtoms === b.polymerAtoms &&
  a.polymerChains === b.polymerChains &&
  a.ligand === b.ligand;

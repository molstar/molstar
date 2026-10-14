/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { structureElementStatsLabel } from '@molstar/model/model/label';
import { Loci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';

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

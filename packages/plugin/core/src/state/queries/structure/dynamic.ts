/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { StructureElement, StructureProperties, type Structure } from '@molstar/model/model/structure';
import {
  AminoAcidNamesL,
  DnaBaseNames,
  RnaBaseNames,
  WaterNames,
  type ElementSymbol,
} from '@molstar/model/model/structure/model/types';
import { ElementNames } from '@molstar/model/model/structure/model/properties/atomic/types';
import { SetUtils } from '@molstar/core/util/set';
import { StructureSelectionQuery } from './query.js';

export function ResidueQuery([names, label]: [string[], string], category: string, priority = 0) {
  const description =
    names.length === 1 && !StandardResidues.has(names[0]) ? `[${names[0]}] ${label}` : `${label} (${names.join(', ')})`;
  return StructureSelectionQuery(
    description,
    MS.struct.modifier.union([
      MS.struct.generator.atomGroups({
        'residue-test': MS.core.set.has([MS.set(...names), MS.ammp('auth_comp_id')]),
      }),
    ]),
    { category, priority, description },
  );
}

export function ElementSymbolQuery([names, label]: [string[], string], category: string, priority: number) {
  const description = `${label} (${names.join(', ')})`;
  return StructureSelectionQuery(
    description,
    MS.struct.modifier.union([
      MS.struct.generator.atomGroups({
        'atom-test': MS.core.set.has([MS.set(...names), MS.acp('elementSymbol')]),
      }),
    ]),
    { category, priority, description },
  );
}

export function EntityDescriptionQuery([names, label]: [string[], string], category: string, priority: number) {
  const description = `${label}`;
  return StructureSelectionQuery(
    `${label}`,
    MS.struct.modifier.union([
      MS.struct.generator.atomGroups({
        'entity-test': MS.core.list.equal([MS.list(...names), MS.ammp('entityDescription')]),
      }),
    ]),
    { category, priority, description },
  );
}

const StandardResidues = SetUtils.unionMany(AminoAcidNamesL, RnaBaseNames, DnaBaseNames, WaterNames);

export function getElementQueries(structures: Structure[]) {
  const uniqueElements = new Set<ElementSymbol>();
  for (const structure of structures) {
    structure.uniqueElementSymbols.forEach((e) => uniqueElements.add(e));
  }

  const queries: StructureSelectionQuery[] = [];
  uniqueElements.forEach((e) => {
    const label = ElementNames[e] || e;
    queries.push(ElementSymbolQuery([[e], label], 'Element Symbol', 0));
  });
  return queries;
}

export function getNonStandardResidueQueries(structures: Structure[]) {
  const residueLabels = new Map<string, string>();
  const uniqueResidues = new Set<string>();
  for (const structure of structures) {
    structure.uniqueResidueNames.forEach((r) => uniqueResidues.add(r));
    for (const m of structure.models) {
      structure.uniqueResidueNames.forEach((r) => {
        const comp = m.properties.chemicalComponentMap.get(r);
        if (comp) residueLabels.set(r, comp.name);
      });
    }
  }

  const queries: StructureSelectionQuery[] = [];
  SetUtils.difference(uniqueResidues, StandardResidues).forEach((r) => {
    const label = residueLabels.get(r) || r;
    queries.push(ResidueQuery([[r], label], 'Ligand/Non-standard Residue', 200));
  });
  return queries;
}

export function getPolymerAndBranchedEntityQueries(structures: Structure[]) {
  const uniqueEntities = new Map<string, string[]>();
  const l = StructureElement.Location.create();
  for (const structure of structures) {
    l.structure = structure;
    for (const ug of structure.unitSymmetryGroups) {
      l.unit = ug.units[0];
      l.element = ug.elements[0];
      const entityType = StructureProperties.entity.type(l);
      if (entityType === 'polymer' || entityType === 'branched') {
        const description = StructureProperties.entity.pdbx_description(l);
        uniqueEntities.set(description.join(', '), description);
      }
    }
  }

  const queries: StructureSelectionQuery[] = [];
  uniqueEntities.forEach((v, k) => {
    queries.push(EntityDescriptionQuery([v, k], 'Polymer/Carbohydrate Entities', 300));
  });
  return queries;
}

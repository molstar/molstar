/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '@molstar/model/model/structure';
import type {
  InteractionElementSchema,
  InteractionInfo,
  StructureInteractionElement,
  StructureInteractions,
} from '@molstar/interactions-extension/model';

export function getCustomInteractionData(
  interactions: InteractionElementSchema[],
  structures: { [ref: string]: Structure },
): StructureInteractions {
  const elements: StructureInteractionElement[] = [];

  for (const schema of interactions) {
    let info: InteractionInfo;
    if (schema.kind === 'covalent') {
      info = { kind: schema.kind, degree: schema.degree };
    } else {
      info = { kind: schema.kind };
    }
    elements.push({
      sourceSchema: schema,
      info,
      aStructureRef: schema.aStructureRef,
      a: resolveLoci(structures[schema.aStructureRef!], schema.a),
      bStructureRef: schema.bStructureRef,
      b: resolveLoci(structures[schema.bStructureRef!], schema.b),
    });
  }

  return { kind: 'structure-interactions', elements };
}

function resolveLoci(structure: Structure, schema: StructureElement.Schema) {
  // (a structure has many interactions, which are looked up in the indices of the model instead of querying every atom for each of them)
  return StructureElement.Schema.toLoci(structure, schema, { useIndex: true });
}

/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { type ElementIndex, Structure, StructureElement } from '@molstar/model/model/structure';
import { structureIntersect } from '@molstar/model/model/structure/query/utils/structure-set';
import {
  StaticStructureComponentTypes,
  createStructureComponent,
} from '@molstar/plugin/state/helpers/structure-component';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { MolScriptBuilder } from '@molstar/query-language/language/builder';
import { Expression } from '@molstar/query-language/language/expression';
import { mapArrayToObject, pickObjectKeys } from '@molstar/core/util/object';
import { substructureFromRows } from '@molstar/mvs/helpers/row-structure';
import type { MVSAnnotationRow } from '@molstar/mvs/helpers/schemas';
import { Choice } from '@molstar/core/util/param-choice';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { capitalize } from '@molstar/core/util/string';
import {
  MVSAnnotationStructureComponentParams,
  createMVSAnnotationStructureComponent,
} from './annotation-structure-component.js';

/** Allowed values for a static selector */
export const StaticSelectorChoice = new Choice(
  mapArrayToObject(StaticStructureComponentTypes, (t) => capitalize(t)),
  'all',
);
export type StaticSelectorChoice = Choice.Values<typeof StaticSelectorChoice>;

/** Parameter definition for specifying a part of structure (kinda extension of `StructureComponentParams` from mol-plugin-state/helpers/structure-component) */
export const SelectorParams = PD.MappedStatic(
  'static',
  {
    static: StaticSelectorChoice.PDSelect(),
    expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.all),
    bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
    script: PD.Script({ language: 'mol-script', expression: '(sel.atom.all)' }),
    annotation: PD.Group(
      pickObjectKeys(MVSAnnotationStructureComponentParams, ['annotationId', 'fieldName', 'fieldValues']),
    ),
    /** Selection criteria (rows of `label_asym_id`, `auth_seq_id`, ... fields), an element is selected if it satisfies any of them.
     * Resolved through the indices of the model, which is much faster than a `expression` equivalent for many rows. */
    schema: PD.Value<MVSAnnotationRow[]>([]),
  },
  {
    description:
      'Define a part of the structure where this layer applies (use Static:all to apply to the whole structure)',
  },
);

/** Parameter values for specifying a part of structure */
export type Selector = PD.Values<{ selector: typeof SelectorParams }>['selector'];

/** `Selector` for selecting the whole structure */
export const SelectorAll = { name: 'static', params: 'all' } satisfies Selector;

/** Decide whether a selector is `SelectorAll` */
export function isSelectorAll(props: Selector): props is typeof SelectorAll {
  return props.name === 'static' && props.params === 'all';
}

/** Data structure for fast lookup of a structure element location in a substructure */
export type ElementSet = { [unitId: number]: SortedArray<ElementIndex> };

export const ElementSet = {
  /** Create an `ElementSet` from a structure */
  fromStructure(structure: Structure | undefined): ElementSet {
    if (!structure) return {};
    const out: ElementSet = {};
    for (const unit of structure.units) {
      out[unit.id] = unit.elements;
    }
    return out;
  },
  /** Create an `ElementSet` from the substructure of `structure` defined by `selector` */
  fromSelector(structure: Structure | undefined, selector: Selector): ElementSet {
    if (!structure) return {};
    const selection = substructureFromSelector(structure, selector); // (rows are resolved through the indices of the model)
    return this.fromStructure(selection);
  },
  /** Decide if the element set `set` contains structure element location `location` */
  has(set: ElementSet, location: StructureElement.Location): boolean {
    const array = set[location.unit.id];
    return array ? SortedArray.has(array, location.element) : false;
  },
};

/** Return a substructure of `structure` defined by `selector` */
export function substructureFromSelector(structure: Structure, selector: Selector): Structure {
  if (selector.name === 'schema') return substructureFromRows(structure, selector.params);
  const pso =
    selector.name === 'annotation'
      ? createMVSAnnotationStructureComponent(structure, { ...selector.params, label: '', nullIfEmpty: false }, {})
      : createStructureComponent(structure, { type: selector, label: '', nullIfEmpty: false }, { source: structure });
  return PluginStateObject.Molecule.Structure.is(pso) ? pso.data : Structure.Empty;
}

/** Select against the root structure, then restrict the result to the supplied substructure. */
export function substructureFromRootSelector(structure: Structure, selector: Selector): Structure {
  return structureIntersect(structure, substructureFromSelector(structure.root, selector));
}

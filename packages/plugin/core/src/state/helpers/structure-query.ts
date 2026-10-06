/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Expression } from '@molstar/model/script/language/expression';
import { type QueryFn, Structure, StructureSelection as Sel, QueryContext } from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import { compile } from '@molstar/model/script/runtime/query/compiler';
import type { PluginStateObject as SO } from '../objects.js';

export { StructureQueryHelper };
namespace StructureQueryHelper {
  export interface CacheEntry {
    script?: Script;
    expression: Expression;
    compiled: QueryFn<Sel>;
    originalStructure: Structure;
    currentStructure: Structure;
  }

  export function isUnchanged(entry: CacheEntry, query: Script | Expression, structure: Structure) {
    if (entry.currentStructure !== structure) return false;
    if (Script.is(query)) {
      return !!entry.script && Script.areEqual(entry.script, query);
    }
    return entry.expression === query;
  }

  export function create(structure: Structure, query: Script | Expression): CacheEntry {
    const script = Script.is(query) ? query : void 0;
    const expression = Script.is(query) ? Script.toExpression(query) : query;
    const compiled = compile<Sel>(expression);

    return { script, expression, compiled, originalStructure: structure, currentStructure: structure };
  }

  export function run(entry: CacheEntry, structure: Structure) {
    return entry.compiled(new QueryContext(structure));
  }

  export function createAndRun(structure: Structure, query: Script | Expression) {
    const entry = create(structure, query);
    return { entry, selection: run(entry, structure) };
  }

  export function updateStructure(entry: CacheEntry, structure: Structure) {
    entry.currentStructure = structure;
    return entry.compiled(new QueryContext(structure));
  }

  export function updateStructureObject(obj: SO.Molecule.Structure, selection: Sel, label?: string) {
    const s = Sel.unionStructure(selection);
    obj.label = `${label || 'Selection'}`;
    obj.description = Structure.elementDescription(s);
    obj.data = s;
  }
}

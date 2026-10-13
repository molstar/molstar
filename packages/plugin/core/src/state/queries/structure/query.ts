/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { QueryContext, StructureSelection } from '@molstar/model/model/structure';
import type { Structure, StructureQuery } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { Expression } from '@molstar/query-language/language/expression';
import { compile } from '@molstar/model/script/runtime/query/compiler';
import type { RuntimeContext } from '@molstar/core/task';

export enum StructureSelectionCategory {
  Type = 'Type',
  Structure = 'Structure Property',
  Atom = 'Atom Property',
  Bond = 'Bond Property',
  Residue = 'Residue Property',
  AminoAcid = 'Amino Acid',
  NucleicBase = 'Nucleic Base',
  Manipulate = 'Manipulate Selection',
  Validation = 'Validation',
  Misc = 'Miscellaneous',
  Internal = 'Internal',
}

export { StructureSelectionQuery };

interface StructureSelectionQuery {
  readonly label: string;
  readonly expression: Expression;
  readonly description: string;
  readonly category: string;
  readonly isHidden: boolean;
  readonly priority: number;
  readonly referencesCurrent: boolean;
  readonly query: StructureQuery;
  readonly ensureCustomProperties?: (ctx: CustomProperty.Context, structure: Structure) => Promise<void>;
  getSelection(plugin: PluginContext, runtime: RuntimeContext, structure: Structure): Promise<StructureSelection>;
}

interface StructureSelectionQueryProps {
  description?: string;
  category?: string;
  isHidden?: boolean;
  priority?: number;
  referencesCurrent?: boolean;
  ensureCustomProperties?: (ctx: CustomProperty.Context, structure: Structure) => Promise<void>;
}

function StructureSelectionQuery(
  label: string,
  expression: Expression,
  props: StructureSelectionQueryProps = {},
): StructureSelectionQuery {
  let _query: StructureQuery;
  return {
    label,
    expression,
    description: props.description || '',
    category: props.category ?? StructureSelectionCategory.Misc,
    isHidden: !!props.isHidden,
    priority: props.priority || 0,
    referencesCurrent: !!props.referencesCurrent,
    get query() {
      if (!_query) _query = compile<StructureSelection>(expression);
      return _query;
    },
    ensureCustomProperties: props.ensureCustomProperties,
    async getSelection(plugin, runtime, structure) {
      const current = plugin.managers.structure.selection.getStructure(structure);
      const currentSelection = current
        ? StructureSelection.Sequence(structure, [current])
        : StructureSelection.Empty(structure);
      if (props.ensureCustomProperties) {
        await props.ensureCustomProperties({ runtime, assetManager: plugin.managers.asset }, structure);
      }
      if (!_query) _query = compile<StructureSelection>(expression);
      return _query(new QueryContext(structure, { currentSelection }));
    },
  };
}

/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { omitObjectKeys, pickObjectKeys } from '../../helpers/utils';
import { ConversionRules, addDefaults, condenseTree, convertTree } from '../generic/tree-utils';
import { MolstarKind, MolstarNode, MolstarTree } from './molstar-tree';
import { FullMVSTree, MVSTree, MVSTreeSchema } from '../mvs/mvs-tree';
import { MVSDefaults } from '../mvs/mvs-defaults';
import { MolstarParseFormatT, ParseFormatT } from '../mvs/param-types';


/** Convert `format` parameter of `parse` node in `MolstarTree`
 * into `format` and `is_binary` parameters in `MolstarTree` */
export const ParseFormatMvsToMolstar = {
    mmcif: { format: 'cif', is_binary: false },
    bcif: { format: 'cif', is_binary: true },
    pdb: { format: 'pdb', is_binary: false },
} satisfies { [p in ParseFormatT]: { format: MolstarParseFormatT, is_binary: boolean } };

/** Conversion rules for conversion from `MVSTree` (with all parameter values) to `MolstarTree` */
const mvsToMolstarConversionRules: ConversionRules<FullMVSTree, MolstarTree> = {
    'download': node => [],
    'parse': (node, parent) => {
        const { format, is_binary } = ParseFormatMvsToMolstar[node.params.format];
        const convertedNode: MolstarNode<'parse'> = { kind: 'parse', params: { ...node.params, format } };
        if (parent?.kind === 'download') {
            return [
                { kind: 'download', params: { ...parent.params, is_binary } },
                convertedNode,
            ] satisfies MolstarNode[];
        } else {
            console.warn('"parse" node is not being converted, this is suspicious');
            return [convertedNode] satisfies MolstarNode[];
        }
    },
    'structure': (node, parent) => {
        if (parent?.kind !== 'parse') throw new Error('Parent of "structure" must be "parse".');
        const { format } = ParseFormatMvsToMolstar[parent.params.format];
        return [
            { kind: 'trajectory', params: { format, ...pickObjectKeys(node.params, ['block_header', 'block_index']) } },
            { kind: 'model', params: pickObjectKeys(node.params, ['model_index']) },
            { kind: 'structure', params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index']) },
        ] satisfies MolstarNode[];
    },
};

/** Node kinds in `MolstarTree` that it makes sense to condense */
const molstarNodesToCondense = new Set<MolstarKind>(['download', 'parse', 'trajectory', 'model'] satisfies MolstarKind[]);

/** Convert MolViewSpec tree into MolStar tree */
export function convertMvsToMolstar(mvsTree: MVSTree): MolstarTree {
    const full = addDefaults<typeof MVSTreeSchema>(mvsTree, MVSDefaults) as FullMVSTree;
    const converted = convertTree<FullMVSTree, MolstarTree>(full, mvsToMolstarConversionRules);
    if (converted.kind !== 'root') throw new Error("Root's type is not 'root' after conversion from MVS tree to Molstar tree.");
    const condensed = condenseTree<MolstarTree>(converted, molstarNodesToCondense);
    return condensed;
}

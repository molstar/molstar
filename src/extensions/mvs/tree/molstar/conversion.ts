/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { omitObjectKeys, pickObjectKeys } from '../../../../mol-util/object';
import { ConversionRules, addDefaults, condenseTree, convertTree, dfs, resolveUris } from '../generic/tree-utils';
import { FullMVSTree, MVSTree, MVSTreeSchema } from '../mvs/mvs-tree';
import { MolstarParseFormatT, ParseFormatT } from '../mvs/param-types';
import { MolstarKind, MolstarNode, MolstarTree } from './molstar-tree';


/** Convert `format` parameter of `parse` node in `MolstarTree`
 * into `format` and `is_binary` parameters in `MolstarTree` */
export const ParseFormatMvsToMolstar = {
    // trajectory
    mmcif: { format: 'cif', is_binary: false },
    bcif: { format: 'cif', is_binary: true },
    pdb: { format: 'pdb', is_binary: false },
    pdbqt: { format: 'pdbqt', is_binary: false },
    gro: { format: 'gro', is_binary: false },
    xyz: { format: 'xyz', is_binary: false },
    mol: { format: 'mol', is_binary: false },
    sdf: { format: 'sdf', is_binary: false },
    mol2: { format: 'mol2', is_binary: false },
    lammpstrj: { format: 'lammpstrj', is_binary: false },
    // coordinates
    xtc: { format: 'xtc', is_binary: true },
    // maps
    map: { format: 'map', is_binary: true },
} satisfies { [p in ParseFormatT]: { format: MolstarParseFormatT, is_binary: boolean } };

/** Conversion rules for conversion from `MVSTree` (with all parameter values) to `MolstarTree` */
const mvsToMolstarConversionRules: ConversionRules<FullMVSTree, MolstarTree> = {
    'download': node => ({ subtree: [] }),
    'parse': (node, parent) => {
        const { format, is_binary } = ParseFormatMvsToMolstar[node.params.format];
        if (parent?.kind === 'download') {
            return {
                subtree: [
                    { kind: 'download', params: { ...parent.params, is_binary }, custom: parent.custom, ref: parent.ref },
                    { kind: 'parse', params: { ...node.params, format }, custom: node.custom, ref: node.ref }
                ] satisfies MolstarNode[]
            };
        } else {
            console.warn('"parse" node is not being converted, this is suspicious');
            return {
                subtree: [
                    { kind: 'parse', params: { ...node.params, format }, custom: node.custom, ref: node.ref }
                ] satisfies MolstarNode[]
            };
        }
    },
    'coordinates': (node, parent) => {
        if (parent?.kind !== 'parse') throw new Error(`Parent of "coordinates" must be "parse", not "${parent?.kind}".`);
        const { format } = ParseFormatMvsToMolstar[parent.params.format];
        return {
            subtree: [
                { kind: 'coordinates', params: { format }, custom: node.custom, ref: node.ref }
            ] satisfies MolstarNode[]
        };
    },
    'structure': (node, parent) => {
        if (parent?.kind !== 'parse') throw new Error(`Parent of "structure" must be "parse", not "${parent?.kind}".`);
        const { format } = ParseFormatMvsToMolstar[parent.params.format];

        if (node.params.coordinates_ref) {
            return {
                subtree: [
                    { kind: 'trajectory', params: { format, ...pickObjectKeys(node.params, ['block_header', 'block_index']) } },
                    { kind: 'model', params: { model_index: 0 } },
                    { kind: 'trajectory_with_coordinates', params: { coordinates_ref: node.params.coordinates_ref } },
                    { kind: 'model', params: pickObjectKeys(node.params, ['model_index']) },
                    { kind: 'structure', params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index', 'coordinates_ref']), custom: node.custom, ref: node.ref },
                ] satisfies MolstarNode[]
            };
        } else {
            return {
                subtree: [
                    { kind: 'trajectory', params: { format, ...pickObjectKeys(node.params, ['block_header', 'block_index']) } },
                    { kind: 'model', params: pickObjectKeys(node.params, ['model_index']) },
                    { kind: 'structure', params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index', 'coordinates_ref']), custom: node.custom, ref: node.ref },
                ] satisfies MolstarNode[]
            };
        }
    },
};

/** Node kinds in `MolstarTree` that it makes sense to condense */
const molstarNodesToCondense = new Set<MolstarKind>(['download', 'parse', 'trajectory', 'model'] satisfies MolstarKind[]);

/** Convert MolViewSpec tree into MolStar tree */
export function convertMvsToMolstar(mvsTree: MVSTree, sourceUrl: string | undefined): MolstarTree {
    const full = addDefaults<typeof MVSTreeSchema>(mvsTree, MVSTreeSchema) as FullMVSTree;
    if (sourceUrl) resolveUris(full, sourceUrl, ['uri', 'url']);
    const converted = convertTree<FullMVSTree, MolstarTree>(full, mvsToMolstarConversionRules);
    if (converted.kind !== 'root') throw new Error("Root's type is not 'root' after conversion from MVS tree to Molstar tree.");
    const condensed = condenseTree<MolstarTree>(converted, molstarNodesToCondense);
    return condensed;
}


type FileExtension = `.${Lowercase<string>}`;
function fileExtensionMatches(filename: string, extensions: (FileExtension | '*')[]): boolean {
    filename = filename.toLowerCase();
    return extensions.some(ext => ext === '*' || filename.endsWith(ext));
}

const StructureFormatExtensions: Record<ParseFormatT, (FileExtension | '*')[]> = {
    // trajectory
    mmcif: ['.cif', '.mmif'],
    bcif: ['.bcif'],
    pdb: ['.pdb', '.ent'],
    pdbqt: ['.pdbqt'],
    gro: ['.gro'],
    xyz: ['.xyz'],
    mol: ['.mol'],
    sdf: ['.sdf'],
    mol2: ['.mol2'],
    lammpstrj: ['.lammpstrj'],
    // coordinates
    xtc: ['.xtc'],
    // volumes
    map: ['.map', '.ccp4', '.mrc', '.mrcs'],
};

/** Run some sanity check on a MVSTree. Return a list of potential problems (`undefined` if there are none) */
export function mvsSanityCheckIssues(tree: MVSTree): string[] | undefined {
    const result: string[] = [];
    dfs(tree, (node, parent) => {
        if (node.kind === 'parse' && parent?.kind === 'download') {
            const source = parent.params.url;
            const extensions = StructureFormatExtensions[node.params.format];
            if (!fileExtensionMatches(source, extensions)) {
                result.push(`Parsing data from ${source} as ${node.params.format} format might be a mistake. The file extension doesn't match recommended file extensions (${extensions.join(', ')})`);
            }
        }
    });
    return result.length > 0 ? result : undefined;
}

/** Run some sanity check on a MVSTree and print potential issues to the console. */
export function mvsSanityCheck(tree: MVSTree): void {
    const issues = mvsSanityCheckIssues(tree);
    if (issues) {
        console.warn('There are potential issues in the MVS tree:');
        for (const issue of issues) {
            console.warn(' ', issue);
        }
    }
}

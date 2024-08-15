/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { treeValidationIssues } from './tree/generic/tree-schema';
import { treeToString } from './tree/generic/tree-utils';
import { Root, createMVSBuilder } from './tree/mvs/mvs-builder';
import { MVSTree, MVSTreeSchema } from './tree/mvs/mvs-tree';


// TODO Clearly divide what belongs to SnapshotMetadata and what to GlobalMetadata

/** Metadata describing global properties relating to the format specification. */
interface FormatMetadata {
    /** Timestamp when this view was exported */
    timestamp: string,
    /** Version of the spec used to write this tree */
    version: string,
    /** Name of this MVSData */
    title?: string,
    /** Detailed description of this view */
    description?: string,
    /** Format of the description */
    description_format?: 'markdown' | 'plaintext',
}
/** Metadata describing details of an individual state/snapshot */
export interface SnapshotMetadata {
    /** Detailed description of this view */
    description?: string,
    /** Format of the description */
    description_format?: 'markdown' | 'plaintext',
    /** Unique identifier of this state, useful when working with collections of states */
    key?: string,
    /** Name of this view */
    title?: string,
    /** How long to linger on one snapshot. Leave empty to not transition automatically */
    lingerDurationMs?: number,
    /** Timespan for the animation to the next snapshot. Leave empty to skip animations */
    transitionDurationMs?: number,
}
interface Snapshot {
    root: MVSTree,
    /** Associated metadata */
    metadata: SnapshotMetadata,
}
export interface MVSData_State extends Snapshot {
    kind?: 'single',
    /** Associated metadata */
    metadata: SnapshotMetadata & FormatMetadata,
}
export interface MVSData_States {
    kind: 'multiple',
    /** MolViewSpec tree */
    snapshots: Snapshot[],
    /** Associated metadata */
    metadata: FormatMetadata,
}

/** Top level of the MolViewSpec (MVS) data format. */
export type MVSData = MVSData_State | MVSData_States

// /** Top level of the MolViewSpec (MVS) data format. */
// export interface MVSData {
//     /** MolViewSpec tree */
//     root: MVSTree,
//     /** Associated metadata */
//     metadata: MVSMetadata,
// }
// interface MVSMetadata {
//     /** Version of the spec used to write this tree */
//     version: string,
//     /** Name of this view */
//     title?: string,
//     /** Detailed description of this view */
//     description?: string,
//     /** Format of the description */
//     description_format?: 'markdown' | 'plaintext',
//     /** Timestamp when this view was exported */
//     timestamp: string,
// }

export const MVSData = {
    /** Currently supported major version of MolViewSpec format (e.g. 1 for version '1.0.8') */
    SupportedVersion: 1,

    /** Parse MVSJ (MolViewSpec-JSON) format to `MVSData`. Does not include any validation. */
    fromMVSJ(mvsjString: string): MVSData {
        const result: MVSData = JSON.parse(mvsjString);
        const major = majorVersion(result?.metadata?.version);
        if (major === undefined) {
            console.error('Loaded MVS does not contain valid version info.');
        } else if (major > (majorVersion(MVSData.SupportedVersion) ?? 0)) {
            console.warn(`Loaded MVS is of higher version (${result.metadata.version}) than currently supported version (${MVSData.SupportedVersion}). Some features may not work as expected.`);
        }
        return result;
    },

    /** Encode `MVSData` to MVSJ (MolViewSpec-JSON) string. Use `space` parameter to control formatting (as with `JSON.stringify`). */
    toMVSJ(mvsData: MVSData, space?: string | number): string {
        return JSON.stringify(mvsData, undefined, space);
    },

    /** Validate `MVSData`. Return `true` if OK; `false` if not OK.
     * If `options.noExtra` is true, presence of any extra node parameters is treated as an issue. */
    isValid(mvsData: MVSData, options: { noExtra?: boolean } = {}): boolean {
        return MVSData.validationIssues(mvsData, options) === undefined;
    },

    /** Validate `MVSData`. Return `undefined` if OK; list of issues if not OK.
     * If `options.noExtra` is true, presence of any extra node parameters is treated as an issue. */
    validationIssues(mvsData: MVSData, options: { noExtra?: boolean } = {}): string[] | undefined {
        const version = mvsData?.metadata?.version;
        if (typeof version !== 'string') return [`MVSData.metadata.version must be a string, not ${typeof version}: ${version}`];
        if (mvsData.kind === 'single' || mvsData.kind === undefined) {
            // if (mvsData.root === undefined) return [`"root" missing in MVS`];
            // return treeValidationIssues(MVSTreeSchema, mvsData.root, options);
            return snapshotValidationIssues(mvsData, options);
        } else if (mvsData.kind === 'multiple') {
            if (mvsData.snapshots === undefined) return [`"snapshots" missing in MVS`];
            const issues: string[] = [];
            for (const snapshot of mvsData.snapshots) { // would use .flatMap if it didn't work in a completely unpredictable way
                const snapshotIssues = snapshotValidationIssues(snapshot, options);
                if (snapshotIssues) issues.push(...snapshotIssues);
            }
            if (issues.length > 0) return issues;
            else return undefined;
        } else {
            return [`MVSData.kind must be 'single' or 'multiple', not ${mvsData.kind}`];
        }
    },

    /** Return a human-friendly textual representation of `mvsData`. */
    toPrettyString(mvsData: MVSData): string {
        const type = mvsData.kind === 'multiple' ? 'multiple states' : 'single state';
        const title = mvsData.metadata.title !== undefined ? ` "${mvsData.metadata.title}"` : '';
        const trees = mvsData.kind === 'multiple' ?
            mvsData.snapshots.map((s, i) => `[Snapshot #${i}]\n${treeToString(s.root)}`).join('\n')
            : treeToString(mvsData.root);
        return `MolViewSpec ${type}${title} (version ${mvsData.metadata.version}, created ${mvsData.metadata.timestamp}):\n${trees}`;
    },

    /** Create a new MolViewSpec builder containing only a root node. Example of MVS builder usage:
     *
     * ```
     * const builder = MVSData.createBuilder();
     * builder.canvas({ background_color: 'white' });
     * const struct = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og2_updated.cif' }).parse({ format: 'mmcif' }).modelStructure();
     * struct.component().representation().color({ color: '#3050F8' });
     * console.log(MVSData.toPrettyString(builder.getState()));
     * ```
     */
    createBuilder(): Root {
        return createMVSBuilder();
    },
};


/** Get the major version from a semantic version string, e.g. '1.0.8' -> 1 */
function majorVersion(semanticVersion: string | number): number | undefined {
    if (typeof semanticVersion === 'string') return parseInt(semanticVersion.split('.')[0]);
    if (typeof semanticVersion === 'number') return Math.floor(semanticVersion);
    console.error(`Version should be a string, not ${typeof semanticVersion}: ${semanticVersion}`);
    return undefined;
}

function snapshotValidationIssues(snapshot: Snapshot, options: { noExtra?: boolean } = {}): string[] | undefined {
    if (snapshot.root === undefined) return [`"root" missing in snapshot`];
    return treeValidationIssues(MVSTreeSchema, snapshot.root, options);
}

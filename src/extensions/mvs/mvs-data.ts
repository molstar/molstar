/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { treeValidationIssues } from './tree/generic/tree-schema';
import { treeToString } from './tree/generic/tree-utils';
import { MVSAnimationSchema, MVSAnimationTree } from './tree/animation/animation-tree';
import { Root, createMVSBuilder } from './tree/mvs/mvs-builder';
import { MVSTree, MVSTreeSchema } from './tree/mvs/mvs-tree';

/** Top-level metadata for a MVS file (single-state or multi-state). */
export interface GlobalMetadata {
    /** Name of this MVSData */
    title?: string,
    /** Detailed description of this view */
    description?: string,
    /** Format of `description`. Default is 'markdown'. */
    description_format?: 'markdown' | 'plaintext',
    /** Timestamp when this view was exported. */
    timestamp: string,
    /** Version of MolViewSpec used to write this file. */
    version: string,
}
export const GlobalMetadata = {
    create(metadata?: Pick<GlobalMetadata, 'title' | 'description' | 'description_format'>): GlobalMetadata {
        return {
            ...metadata,
            version: `${MVSData.SupportedVersion}`,
            timestamp: utcNowISO(),
        };
    },
};

/** Metadata for an individual snapshot. */
export interface SnapshotMetadata {
    /** Name of this snapshot. */
    title?: string,
    /** Detailed description of this snapshot. */
    description?: string,
    /** Format of `description`. Default is 'markdown'. */
    description_format?: 'markdown' | 'plaintext',
    /** Unique identifier of this state, useful when working with collections of states. */
    key?: string,
    /** Timespan for snapshot. */
    linger_duration_ms: number,
    /** Timespan for the animation to the next snapshot. Leave empty to skip animations. */
    transition_duration_ms?: number,
}

export interface Snapshot {
    /** Root of the node tree */
    root: MVSTree,
    /** Associated metadata */
    metadata: SnapshotMetadata,
    /** Optional animation */
    animation?: MVSAnimationTree,
}

/** MVSData with a single state */
export interface MVSData_State {
    kind?: 'single',
    /** Root of the node tree */
    root: MVSTree,
    /** Associated metadata */
    metadata: GlobalMetadata,
}

/** MVSData with multiple states (snapshots) */
export interface MVSData_States {
    kind: 'multiple',
    /** Ordered collection of individual states */
    snapshots: Snapshot[],
    /** Associated metadata */
    metadata: GlobalMetadata,
}

/** Top level of the MolViewSpec (MVS) data format. */
export type MVSData = MVSData_State | MVSData_States


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

    /** Create a multi-state MVS data from a list of snapshots. */
    createMultistate(snapshots: Snapshot[], metadata?: Pick<GlobalMetadata, 'title' | 'description' | 'description_format'>): MVSData_States {
        return {
            kind: 'multiple',
            snapshots: [...snapshots],
            metadata: GlobalMetadata.create(metadata),
        };
    },

    /** Convert single-state MVSData into multi-state MVSData with one state. */
    stateToStates(state: MVSData_State): MVSData_States {
        return {
            kind: 'multiple',
            metadata: state.metadata,
            snapshots: [{
                metadata: {
                    title: state.metadata.title,
                    description: state.metadata.description,
                    description_format: state.metadata.description_format,
                    key: undefined,
                    linger_duration_ms: 1000,
                    transition_duration_ms: 250,
                },
                root: state.root,
            }],
        };
    },
};


/** Get the major version from a semantic version string, e.g. '1.0.8' -> 1 */
function majorVersion(semanticVersion: string | number): number | undefined {
    if (typeof semanticVersion === 'string') return parseInt(semanticVersion.split('.')[0]);
    if (typeof semanticVersion === 'number') return Math.floor(semanticVersion);
    console.error(`Version should be a string, not ${typeof semanticVersion}: ${semanticVersion}`);
    return undefined;
}

function snapshotValidationIssues(snapshot: MVSData_State | Snapshot, options: { noExtra?: boolean } = {}): string[] | undefined {
    if (snapshot.root === undefined) return [`"root" missing in snapshot`];
    const state = treeValidationIssues(MVSTreeSchema, snapshot.root, options);
    const animation = 'animation' in snapshot && snapshot.animation !== undefined
        ? treeValidationIssues(MVSAnimationSchema, snapshot.animation, options)
        : undefined;
    if (state && animation) return [...state, ...animation];
    if (state) return state;
    if (animation) return animation;
    return undefined;
}

/** Return the current universal time, in ISO format, e.g. '2023-11-24T10:45:49.873Z' */
function utcNowISO(): string {
    return new Date().toISOString();
}

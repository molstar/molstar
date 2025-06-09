/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { dfs } from 'molviewspec/lib/tree/generic/tree-utils';
import { MVSTree } from 'molviewspec/lib/tree/mvs/mvs-tree';
import { MolstarTree } from './molstar-tree';

/** Convert a MVS tree to a Molstar tree (i.e. expand it with Molstar-specific nodes). */
export function convertMvsToMolstar(tree: MVSTree, sourceUrl?: string): MolstarTree {
    const result: MolstarTree = {
        kind: 'root',
        children: [],
    };

    if (tree.kind === 'root') {
        result.children = tree.children?.map(child => convertMvsNode(child, sourceUrl)) || [];
    } else {
        result.children = [convertMvsNode(tree, sourceUrl)];
    }

    return result;
}

function convertMvsNode(node: any, sourceUrl?: string): any {
    // Handle download nodes with URL resolution
    if (node.kind === 'download') {
        return {
            ...node,
            params: {
                ...node.params,
                url: resolveUrl(node.params.url, sourceUrl),
            },
        };
    }

    // For all other nodes, recursively convert children
    const result = { ...node };
    if (node.children) {
        result.children = node.children.map((child: any) => convertMvsNode(child, sourceUrl));
    }

    return result;
}

function resolveUrl(url: string, baseUrl?: string): string {
    if (!baseUrl || url.match(/^https?:\/\//)) {
        return url;
    }

    try {
        return new URL(url, baseUrl).href;
    } catch {
        return url;
    }
}

/** Run sanity checks on a MVS tree and print potential issues to the console. */
export function mvsSanityCheck(tree: MVSTree): void {
    const issues: string[] = [];

    dfs(tree, (node) => {
        // Check for common issues
        if (node.kind === 'download' && !node.params?.url) {
            issues.push(`Download node missing URL: ${JSON.stringify(node)}`);
        }

        if (node.kind === 'component' && node.params?.selector === undefined) {
            issues.push(`Component node with undefined selector: ${JSON.stringify(node)}`);
        }
    });

    if (issues.length > 0) {
        console.warn('MVS Sanity Check Issues:');
        issues.forEach(issue => console.warn(`  - ${issue}`));
    }
}

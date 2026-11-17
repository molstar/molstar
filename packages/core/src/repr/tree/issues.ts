/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The coded issue (_spec/README.md "Issues"): one shape everywhere, a stable namespaced code callers branch on, a
// severity with one meaning, the place as a step path and, where the place has one, an address. Written from meso*'s
// tree core, with the severity added and the codes moved under `repr.*`.

/** A position in a tree: keys and child indices from the root (`['children', 3, 'params', 'count']`). */
export type TreePath = readonly (string | number)[];

/**
 * One meaning everywhere: `error` the input is invalid, `warning` something is lost (kept, skipped, not reported),
 * `info` never changes what is drawn.
 */
export type Severity = 'error' | 'warning' | 'info';

export interface Issue {
  /** Stable and namespaced (`repr.codec.json.syntax`, `repr.schema.param.type`); callers branch on it, never on `message`. */
  readonly code: string;
  /** The code's severity, from its catalog. */
  readonly severity: Severity;
  /** Human text, free to improve; never quotes the input (an offset, a path or an address says where). */
  readonly message: string;
  /** The actionable next step: the names that would fit, the expected shape. */
  readonly hint?: string;
  /** Where in the text (a parse issue): 0-based offset. */
  readonly offset?: number;
  /** Where in the tree, step by step from the root (`['children', 3, 'params', 'count']`). */
  readonly path?: TreePath;
  /** Where in the tree, by address (`/data/p53`), when the place has one. */
  readonly address?: string;
}

export interface IssueCodeInfo {
  readonly severity: Severity;
  readonly description: string;
  /** The usual next step, when one fits every occurrence. */
  readonly hint?: string;
  /** Retired codes stay reserved: a code is API. */
  readonly retired?: boolean;
}

/** A module's catalog of codes, all in one namespace. */
export interface IssueCatalog<C extends string = string> {
  readonly namespace: string;
  readonly codes: Readonly<Record<C, IssueCodeInfo>>;
}

export function defineIssues<C extends string>(
  namespace: string,
  codes: Readonly<Record<C, IssueCodeInfo>>,
): IssueCatalog<C> {
  for (const code of Object.keys(codes)) {
    if (!code.startsWith(`${namespace}.`)) throw new Error(`issue code ${code} is outside its namespace ${namespace}`);
  }
  return { namespace, codes };
}

/** What an issue says besides its code, severity and message. */
export interface IssuePlace {
  readonly hint?: string;
  readonly offset?: number;
  readonly path?: TreePath;
  readonly address?: string;
}

/** An issue of a cataloged code, its severity the catalog's; absent places are left out. */
export function issueOf<C extends string>(
  catalog: IssueCatalog<C>,
  code: C,
  message: string,
  place: IssuePlace = {},
): Issue {
  const out: { -readonly [K in keyof Issue]: Issue[K] } = { code, severity: catalog.codes[code].severity, message };
  if (place.hint !== undefined) out.hint = place.hint;
  if (place.offset !== undefined) out.offset = place.offset;
  if (place.path !== undefined) out.path = place.path;
  if (place.address !== undefined) out.address = place.address;
  return out;
}

/** Whether any issue is an error (the input is invalid). */
export function hasErrors(issues: readonly Issue[]): boolean {
  for (const i of issues) if (i.severity === 'error') return true;
  return false;
}

/** The tree's codes: nodes, values, columns, addresses and the readers' depth bound. */
export const TreeIssues = defineIssues('repr.tree', {
  'repr.tree.depth_exceeded': {
    severity: 'error',
    description: 'A text nests deeper than the readers build (MaxDepth levels)',
  },
  'repr.tree.node.shape': {
    severity: 'error',
    description: 'A node is not an object with a `kind`, or has a key other than kind, id, params, custom, children, args',
  },
  'repr.tree.node.kind': {
    severity: 'error',
    description: 'A kind is not a word (a lowercase letter, then a-z, 0-9, _ or -), optionally after a namespace and `:`',
  },
  'repr.tree.node.id': {
    severity: 'error',
    description: 'An id is not one address segment: a letter, digit or _, then letters, digits, _, ., + or -',
  },
  'repr.tree.value.shape': {
    severity: 'error',
    description: 'A value is none of: string, number, bigint, boolean, null, list, record, node, typed column',
  },
  'repr.tree.value.reserved_key': {
    severity: 'error',
    description: 'A param or record key is `kind`, `__proto__`, or starts with `$` (tagged atoms)',
  },
  'repr.tree.value.tag': { severity: 'error', description: 'A tagged atom (`$number`, `$bigint`) has a malformed payload' },
  'repr.tree.address.syntax': {
    severity: 'error',
    description: 'Not an address: /, then ids or query forms separated by /, optionally ending in /*',
  },
  'repr.tree.address.form': {
    severity: 'error',
    description: 'A query segment is not (kind :key value …): it has an id or children, or a record value',
  },
  'repr.tree.address.unknown_kind': {
    severity: 'error',
    description: "A query segment names a kind the language (the caller's kinds) lacks",
  },
  'repr.tree.address.unknown_param': {
    severity: 'error',
    description: "A query segment names a param its kind lacks (the caller's params): it could match nothing",
  },
  'repr.tree.address.resolver': {
    severity: 'error',
    description:
      'A segment cannot apply to what its prefix names: an id, a `/*` or a form without a resolver after a non-node result, or a resolver that failed',
  },
  'repr.tree.address.cancelled': { severity: 'error', description: "A selection was cancelled by the caller's signal" },
  'repr.tree.address.deadline_exceeded': {
    severity: 'error',
    description: "A selection ran past the caller's deadline",
  },
  'repr.tree.address.option': {
    severity: 'error',
    description: 'A selection option is out of range: a limit or maxDepth that is not a non-negative integer or Infinity',
  },
  'repr.tree.address.limit_exceeded': {
    severity: 'error',
    description: 'An address names more nodes than the caller allows (its limit)',
  },
  'repr.tree.address.depth_exceeded': {
    severity: 'error',
    description: 'An address has more than MaxDepth segments, or its walk goes deeper than the caller allows',
  },
});

/** The issue paths as JSON pointers (`/children/3/params/count`), for messages and logs. */
export function pathToPointer(path: TreePath): string {
  return path.map((p) => `/${String(p).replaceAll('~', '~0').replaceAll('/', '~1')}`).join('');
}

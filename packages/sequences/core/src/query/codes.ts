/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's issue codes and their severities (seqspec:seqql#errors): one severity per code, wherever it is reported.
// An error is a query written wrong or one that cannot run; a warning, a part of what was written that selects
// nothing, or a result that may be incomplete; info, a note on how a valid query ran (why its result is empty).
import type { QueryIssue } from './types.js';

const ERRORS = `read.unclosed read.extra_close read.trailing read.string read.reserved read.keyword read.empty read.depth
  read.size read.head reply.empty reply.nested reply.none reply.unclosed reply.unopened symbol.unknown
  symbol.duplicate symbol.namespace symbol.unavailable alias.collision alias.reserved arg.count arg.missing arg.order
  arg.repeated arg.unknown type slot sel.syntax sel.mixed macro.depth current field.unknown slice.range ref.unknown
  ref.kind numbering.unknown numbering.ambiguous classes.unknown input.unknown limit.positions limit.runs limit.groups limit.time aborted`;
const WARNINGS = `lazy letters.columns label.unknown ids.unknown onto.kind onto.ambiguous onto.unmapped field.near
  slice.clipped`;
const INFOS = 'empty';

const table = (names: string, severity: QueryIssue['severity']) =>
  names
    .trim()
    .split(/\s+/)
    .map((c) => ({ code: `seqql.${c}`, severity }));

/** Every SeqQL code with its severity, in the reference's order. */
export const SeqqlIssueCodes: readonly { readonly code: string; readonly severity: QueryIssue['severity'] }[] = [
  ...table(ERRORS, 'error'),
  ...table(WARNINGS, 'warning'),
  ...table(INFOS, 'info'),
];

const SEVERITY: ReadonlyMap<string, QueryIssue['severity']> = new Map(SeqqlIssueCodes.map((c) => [c.code, c.severity]));

/** A code's registered severity; an unlisted code (an extension's own) is an error. */
export function seqqlSeverity(code: string): QueryIssue['severity'] {
  return SEVERITY.get(code) ?? 'error';
}

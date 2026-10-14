/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seq-core/query: SeqQL, the language of Seq* queries (seqspec:seqql#overview). A leaf: the types,
// the base table with full and short names, expressions and their JSON, the registry and the builder (QL1); the
// reader, the binder and the canonical printer (QL2a); the checker, `sel`, and the generated skill and reference
// (QL2b); queries in replies (seqspec:seqql#queries-in-text).
export { bindQuery, type QueryBinding } from "./bind.ts";
export { Q, type QueryError, type SymbolPath, type SymbolTree, type UnionToIntersection } from "./builder.ts";
export { checkQuery, checkText, type QueryCheck } from "./check.ts";
export { SeqqlIssueCodes, seqqlSeverity } from "./codes.ts";
export { extractQueries, type QueryBlock, type QueryExtraction, SEQQL_CLOSE, SEQQL_OPEN } from "./extract.ts";
export { type QueryIssuePlace, queryIssue, renderIssue, suggest } from "./issue.ts";
export { createLanguage, type Language, type LanguageExtension } from "./language.ts";
export { printQuery } from "./print.ts";
export { MAX_DEPTH, NUMBER, type QueryNode, type QueryRead, readQuery, type TextPosition } from "./read.ts";
export { GUIDANCE, referenceText, type SkillOptions, skillText } from "./skill.ts";
export { defineSymbol } from "./symbol.ts";
export { BASE_SYMBOLS, type BaseSymbolName } from "./table.ts";
export {
  type Apply,
  type ArgDef,
  type Expression,
  expressionTooDeep,
  isApply,
  isExpression,
  type Keywords,
  MAX_EXPRESSION_DEPTH,
  type MacroReport,
  type QueryIssue,
  SEQQL_VERSION,
  type SeqQLQuery,
  type Slot,
  type SymbolDef,
  type SymbolExample,
  type SymbolHint,
  type SymbolKind,
  type SymbolSpec,
  type Type,
} from "./types.ts";

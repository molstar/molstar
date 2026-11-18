/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Applying a patch (_spec/README.md "Patches", Applying), in order: the checks (result, then base); the ops, all or
// nothing, by structural sharing; normalization of what they made; refs following a move or a rename; validation of
// what changed, by how each object came to be; link integrity, the old root against the new; the language's
// validators; the result check. The caller's root is never touched.

import {
  type Language,
  NoChanges,
  normalize,
  normalizeValue,
  validateChanges,
  type ValidateOptions,
} from '../schema/index.js';
import {
  type Address,
  hasErrors,
  hashNode,
  type Issue,
  issueOf,
  type Node,
  own,
  printAddress,
  type TreePath,
  type Value,
} from '../tree/index.js';
import {
  type ApplyOptions,
  MaxOps,
  type Op,
  type Patch,
  type PatchCode,
  PatchIssues,
  type PatchResult,
} from './patch.js';
import { RefPass, refIndexOf, rewriteRef, type Site } from './refs.js';
import { type Draft, type Normalizer, Work } from './work.js';

const refuse = (...issues: Issue[]): PatchResult => ({ ok: false, issues });
const issue = (code: PatchCode, message: string, path: TreePath) =>
  issueOf(PatchIssues, code, message, { path, hint: PatchIssues.codes[code].hint });

function normalizerOf(language: Language, options: ValidateOptions): Normalizer {
  const languages = (id: string) => (id === language.id ? language : options.languages?.(id));
  return {
    node: (n) => normalize(language, n, { languages }),
    value: (kind, key, v) => {
      const spec = own(language.kinds, kind);
      const t = spec?.params ? own(spec.params, key) : undefined;
      return t ? normalizeValue(t, v, { languages }) : v;
    },
  };
}

/** Applies a patch to a root: a new revision with the patch's inverse and what changed, or the issues refusing it. */
export function applyPatch(root: Node, patch: Patch, options: ApplyOptions = {}): PatchResult {
  const revision = hashNode(root);
  const language = options.language;
  const refs = language ? {} : { refs: 'not-rewritten' as const };
  if (patch.result !== undefined && patch.result === revision)
    return { ok: true, changed: false, root, revision, changes: NoChanges, issues: [], ...refs };
  if (patch.base !== undefined && patch.base !== revision)
    return refuse(issue('repr.patch.stale', 'The patch was made against another revision', ['params', 'base']));
  if (options.requireBase && patch.base === undefined)
    return refuse(issue('repr.patch.base', 'The patch names no revision it was made against', ['params']));
  if (!Array.isArray(patch.ops)) return refuse(issue('repr.patch.op', 'A patch has a list of ops', ['children']));
  const maxOps = options.maxOps ?? MaxOps;
  if (patch.ops.length > maxOps) return refuse(issue('repr.patch.limit', `More than ${maxOps} ops`, ['children']));

  const validate = options.validate ?? {};
  const work = new Work(root, language ? normalizerOf(language, validate) : undefined);
  let moves = 0;
  for (let i = 0; i < patch.ops.length; i++) {
    const failed = work.apply(patch.ops[i], i);
    if (failed) return refuse(failed);
    // refs follow a move or a rename, op by op: a ref an earlier op wrote follows a later move
    for (; language && moves < work.moves.length; moves++) {
      const { from, to } = work.moves[moves];
      const index = refIndexOf(work.finalize(), language, options.validate);
      for (const site of index.under(from)) {
        const next = rewriteRef(site.value, from, to);
        if (next !== site.value) work.rewrite(site.path, next);
      }
    }
  }

  let after = work.finalize();
  let next = hashNode(after);
  const changed = next !== revision;
  if (!changed) {
    after = root;
    next = revision;
  }
  const changes = changed ? work.changes() : NoChanges;
  const issues: Issue[] = [];
  let dangled = false;
  const rebound: Op[] = [];
  if (language && changed) {
    const found = validateChanges(language, after, changes, validate);
    if (hasErrors(found)) return refuse(...found);
    issues.push(...found);
    const links = linkIntegrity(root, after, work, language, options);
    if (hasErrors(links.issues)) return refuse(...links.issues);
    issues.push(...links.issues);
    dangled = links.dangled;
    rebound.push(...links.fixups);
  }
  if (changed) {
    for (const v of options.validators ?? []) {
      const found = v(root, after, changes);
      if (hasErrors(found)) return refuse(...found);
      issues.push(...found);
    }
  }
  if (patch.result !== undefined && patch.result !== next)
    return refuse(issue('repr.patch.result', 'The patch gave another revision than it says', ['params', 'result']));
  if (!changed) return { ok: true, changed: false, root, revision, changes, issues, ...refs };

  const inverse: Patch = { ops: [...work.inverse].reverse().concat(rebound), base: next, result: revision };
  const needs: { allowDangling?: true; allowRebind?: true } = {};
  if (dangled) needs.allowRebind = true;
  if (rebound.length) needs.allowDangling = true;
  return {
    ok: true,
    changed: true,
    root: after,
    revision: next,
    inverse,
    ...(Object.keys(needs).length ? { inverseNeeds: needs } : {}),
    changes,
    issues,
    ...refs,
  };
}

/**
 * Link integrity (Applying, step 6), judged on the old root against the new. Dangling: a ref of the new root that
 * names no node, and named one in the old root or was written by this patch. Rebound: a ref the patch did not change
 * that named no node in the old root and names one in the new. Without an op that changes an address, only the refs
 * the patch wrote are read.
 */
function linkIntegrity(
  before: Node,
  after: Node,
  work: Work,
  language: Language,
  options: ApplyOptions,
): { issues: Issue[]; dangled: boolean; fixups: Op[] } {
  const oldPass = new RefPass(before, language, options.validate);
  const newPass = new RefPass(after, language, options.validate);
  const issues: Issue[] = [];
  const fixups: Op[] = [];
  let dangled = false;

  const visit = (
    n: Node,
    d: Draft | undefined,
    path: TreePath,
    at: Address | undefined,
    origin: Node | undefined,
    oldAddress: Address | undefined,
    made: boolean,
  ) => {
    const written = d?.written;
    if (work.structural || written) {
      const sites = newPass.nodeSites(n, path, at);
      let old: Map<string, Site> | undefined;
      const oldOf = (rel: TreePath) => {
        if (!origin) return undefined;
        old ??= new Map(oldPass.nodeSites(origin, path, oldAddress).map((s) => [s.path.join('\u0000'), s]));
        return old.get(rel.join('\u0000'));
      };
      for (const s of sites) {
        const key = s.path[path.length + 1] as string;
        const wrote = made || !!written?.has(key);
        if (!work.structural && !wrote) continue;
        const names = newPass.names(s);
        const was = wrote ? undefined : oldOf(s.path);
        if (!names) {
          if (wrote || !was || oldPass.names(was)) {
            dangled = true;
            const place = { path: s.path, ...(at ? { address: printAddress(at) } : {}) };
            const i = issueOf(PatchIssues, 'repr.patch.dangling', 'The ref names no node', {
              ...place,
              hint: PatchIssues.codes['repr.patch.dangling'].hint,
            });
            issues.push(options.allowDangling ? { ...i, severity: 'warning' } : i);
          }
        } else if (was && was.value === s.value && !oldPass.names(was)) {
          if (!options.allowRebind) {
            issues.push(
              issueOf(PatchIssues, 'repr.patch.rebound', 'A ref that named no node names one now', {
                path: s.path,
                ...(at ? { address: printAddress(at) } : {}),
                hint: PatchIssues.codes['repr.patch.rebound'].hint,
              }),
            );
          } else if (oldAddress && origin?.params) {
            // its value restored after the moves back, which would take it along with the refs they rewrite
            const value = own(origin.params, key) as Value;
            if (!fixups.some((f) => f.op === 'set' && f.key === key && printAddress(f.at) === printAddress(oldAddress)))
              fixups.push({ op: 'set', at: oldAddress, key, value });
          }
        }
      }
    }
    const children = n.children ?? [];
    const kids = d?.kids;
    if (!work.structural && !kids) return;
    children.forEach((c, i) => {
      const k = kids?.[i];
      const cAt = at && c.id !== undefined ? [...at, c.id] : undefined;
      // an unmaterialized child is the old root's own when its parent's children are
      const same = d ? d.base === d.origin : origin === n;
      const cOrigin = k ? k.origin : same ? c : undefined;
      const cOld = k ? k.oldAddress : oldAddress && c.id !== undefined && cOrigin ? [...oldAddress, c.id] : undefined;
      visit(c, k, [...path, 'children', i], cAt, cOrigin, cOld, k ? k.made : made);
    });
  };
  const r = work.root.draft;
  visit(after, r, [], [], r.origin, r.oldAddress, r.made);
  return { issues, dangled, fixups };
}

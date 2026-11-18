/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Applying a patch (_spec/README.md "Patches", Applying), in order: the checks (result, then base); the ops, all or
// nothing, by structural sharing; normalization of what they made; refs following a move or a rename; validation of
// what changed, by how each object came to be; link integrity, the old root against the new; the language's
// validators; the result check. The caller's root is never touched.

import { NoChanges } from '../schema/changes.js';
import {
  type Changes,
  type Language,
  normalize,
  normalizeValue,
  type ValidateOptions,
  validateChanges,
} from '../schema/index.js';
import {
  type Address,
  hasErrors,
  hashNode,
  type Issue,
  isNode,
  issueOf,
  type Node,
  own,
  printAddress,
  resolveAddress,
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
import { deferDerive, LiveRefs, type LiveSite, RefPass, rewriteRef, type Step } from './refs.js';
import { type Draft, type Normalizer, sameAddress, Work } from './work.js';

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
  // one ref index for the apply, derived across the ops
  const live = language ? new LiveRefs(work, language, options.validate) : undefined;
  work.listener = live;
  let moves = 0;
  for (let i = 0; i < patch.ops.length; i++) {
    if (live && patch.ops[i]?.op === 'move') live.retype();
    const failed = work.apply(patch.ops[i], i);
    if (failed) return refuse(failed);
    // refs follow a move or a rename, op by op: a ref an earlier op wrote follows a later move
    for (; live && moves < work.moves.length; moves++) live.move(work.moves[moves].from, work.moves[moves].to);
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
  const needs: { allowDangling?: true; allowRebind?: true } = {};
  let fixups: Op[] = [];
  if (language && live && changed) {
    const found = validateChanges(language, after, changes, validate);
    if (hasErrors(found)) return refuse(...found);
    issues.push(...found);
    const links = linkIntegrity(after, work, live, changes, options);
    if (hasErrors(links.issues)) return refuse(...links.issues);
    issues.push(...links.issues);
    // the allowances turned around: refs left dangling rebind on the way back, and the other way
    if (links.dangled) needs.allowRebind = true;
    if (links.rebound) needs.allowDangling = needs.allowRebind = true;
    const restoring = restoringSets(root, work, live);
    fixups = restoring.fixups;
    if (restoring.dangling) needs.allowDangling = true;
    if (restoring.rebinds) needs.allowRebind = true;
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

  // the new root's ref index is derived from the old root's, when the old one exists, once the next apply reads it:
  // a chain of applies never builds from scratch, and a run of sets that read no index derives nothing
  if (language) deferDerive(root, after, language, options.validate);

  const inverse: Patch = { ops: [...work.inverse].reverse().concat(fixups), base: next, result: revision };
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

/** The value at `rel` (`['params', key, …]`) in a node, or undefined. */
function valueAt(n: Node | undefined, rel: TreePath): Value | undefined {
  let v: unknown = n;
  for (const k of rel) {
    if (v === null || typeof v !== 'object') return undefined;
    v = Array.isArray(v) ? v[k as number] : Object.hasOwn(v, k) ? (v as Record<string, unknown>)[k] : undefined;
  }
  return v as Value | undefined;
}

/** A stored patch's step as the old root had it: its node there (a rewrite made the patch's node anew). */
function oldStep(h: Draft, step: Step | undefined): Step | undefined {
  if (!step) return undefined;
  const slot = valueAt(h.origin, step.rel);
  return isNode(slot) ? { ...step, slot } : step;
}

/** A draft's step path in the new root. */
function stepPath(d: Draft): TreePath {
  const out: (string | number)[] = [];
  for (let at = d; at.parent; at = at.parent) out.unshift('children', at.parent.kids!.indexOf(at));
  return out;
}

/**
 * Link integrity (Applying, step 6), judged on the old root against the new, through the apply's ref index: the refs
 * the ops wrote, made or rewrote, and those of the old root under an address the patch changed (removed, made, moved
 * from or to). Dangling: a ref that names no node, and named one in the old root or was written by this patch (a ref
 * that was in the old root with the same value at the same place and dangled there is not written). Rebound: a ref
 * the patch did not change whose target changed identity: it named no node and names one now, or the node it named
 * was removed or moved away and its address names another node.
 */
function linkIntegrity(
  after: Node,
  work: Work,
  live: LiveRefs,
  changes: Changes,
  options: ApplyOptions,
): { issues: Issue[]; dangled: boolean; rebound: boolean } {
  const oldPass = live.oldPass;
  const newPass = new RefPass(after, live.language, options.validate);
  const issues: Issue[] = [];
  let dangled = false;
  let rebound = false;

  const report = (code: 'repr.patch.dangling' | 'repr.patch.rebound', message: string, s: LiveSite) => {
    const at = work.addressOf(s.holder);
    const i = issueOf(PatchIssues, code, message, {
      path: [...stepPath(s.holder), ...s.rel],
      ...(at ? { address: printAddress(at) } : {}),
      hint: PatchIssues.codes[code].hint,
    });
    if (code === 'repr.patch.dangling') {
      dangled = true;
      issues.push(options.allowDangling ? { ...i, severity: 'warning' } : i);
    } else {
      rebound = true;
      if (!options.allowRebind) issues.push(i);
    }
  };
  const names = (pass: RefPass, run: Address, step?: Step) => pass.names(run, step);
  /** A ref the patch did not change: dangling if it named a node, rebound if its target changed identity. */
  const untouched = (s: LiveSite, before: boolean) => {
    const now = names(newPass, s.run, s.step);
    if (!now && before) report('repr.patch.dangling', 'The ref names no node', s);
    else if (now && !before) report('repr.patch.rebound', 'A ref that named no node names one now', s);
    else if (now && before && !s.step) {
      const origin = work.originAt(s.run);
      if (!origin || !sameAddress(origin, s.run))
        report('repr.patch.rebound', 'The node the ref named is gone, and its address names another', s);
    }
  };

  for (const s of live.written()) {
    const h = s.holder;
    const key = String(s.rel[1]);
    if (h.made) {
      const at = work.addressOf(h);
      const old = at && valueAt(resolveAddress(work.old, at), s.rel);
      // the same ref at the same place, dangling before: not introduced
      if (old === s.value && !names(oldPass, s.run, s.step)) untouched(s, false);
      else if (!names(newPass, s.run, s.step)) report('repr.patch.dangling', 'The ref names no node', s);
    } else if (h.written?.has(key)) {
      if (!names(newPass, s.run, s.step)) report('repr.patch.dangling', 'The ref names no node', s);
    } else {
      // a key a move rewrote: its refs the move left alone are untouched
      const old = valueAt(h.origin, s.rel);
      if (old === s.value) untouched(s, names(oldPass, s.run, oldStep(h, s.step)));
      else {
        // a ref a move took along: judged by what it named before and names now, as one left alone
        const was = typeof old === 'string' ? live.runOf(old) : undefined;
        const before = !!was && names(oldPass, was, oldStep(h, s.step));
        const now = names(newPass, s.run, s.step);
        if (!now && (before || !was)) report('repr.patch.dangling', 'The ref names no node', s);
        else if (now && !before) report('repr.patch.rebound', 'A ref that named no node names one now', s);
        else if (now && was && !s.step) {
          const origin = work.originAt(s.run);
          if (!origin || !sameAddress(origin, was))
            report('repr.patch.rebound', 'The node the ref named is gone, and its address names another', s);
        }
      }
    }
  }
  if (work.structural) {
    const seen = new Map<Draft, Set<string>>();
    const prefixes = [...changes.removed, ...changes.made, ...changes.moved.flatMap((m) => [m.from, m.to])];
    for (const p of prefixes) {
      for (const s of live.under(p, true, false)) {
        let set = seen.get(s.holder);
        if (!set) seen.set(s.holder, (set = new Set()));
        const k = s.rel.join('\u0000');
        if (set.has(k)) continue;
        set.add(k);
        untouched(s, names(oldPass, s.run, s.step));
      }
    }
  }
  return { issues, dangled, rebound };
}

/**
 * The nodes of a subtree taken out (removed, or made over by a replace) that the old root had: each as it was when
 * taken out, with its old root's node and address.
 */
function restored(work: Work, d: Draft, visit: (node: Node, origin: Node, at: Address) => void) {
  const nodes = (n: Node, at: Address) => {
    visit(n, n, at);
    for (const c of n.children ?? []) if (c.id !== undefined) nodes(c, [...at, c.id]);
  };
  const drafts = (x: Draft) => {
    if (!x.kids && x.origin && x.base === x.origin && x.oldAddress) {
      // untouched below: the old root's own nodes, the draft's params aside
      visit(work.finalize(x), x.origin, x.oldAddress);
      for (const c of x.base.children ?? []) if (c.id !== undefined) nodes(c, [...x.oldAddress, c.id]);
      return;
    }
    if (x.origin && x.oldAddress && !x.made) visit(work.finalize(x), x.origin, x.oldAddress);
    for (const k of x.kids ?? []) drafts(k);
  };
  drafts(d);
}

/**
 * The inverse's restoring sets (The inverse and the allowances): every ref the inverse's moves would take along
 * without the forward patch having rewritten it is set back at the end of the inverse. Each ref is followed from where
 * the inverse has it (the new root, or the op that restores its value or node) through the inverse's later moves; one
 * that does not end as it was in the old root gets a `set` of its key's old value at its node's old address. Also
 * whether the inverse writes a ref that names no node (content restored with refs that dangled): it needs
 * `allowDangling`.
 */
function restoringSets(
  before: Node,
  work: Work,
  live: LiveRefs,
): { fixups: Op[]; dangling: boolean; rebinds: boolean } {
  const n = work.inverse.length;
  const moves: { pos: number; from: Address; to: Address }[] = [];
  work.inverseMeta.forEach((m, i) => {
    if (m?.kind === 'move') moves.push({ pos: n - 1 - i, from: m.from, to: m.to });
  });
  moves.sort((a, b) => a.pos - b.pos);
  const oldPass = live.oldPass;
  // the addresses in stored patches' `set` values, untyped: the inverse may type them over another state
  const loosePass = new RefPass(work.old, live.language, live.options);
  loosePass.loose = true;
  type Entry = { pos: number; sites: { value: string; run: Address; step?: Step; expected: Value | undefined }[] };
  /** By the old address of a node and a key: where the inverse last gives that key its value. */
  const entries = new Map<string, { at: Address; key: string; entry: Entry }>();
  const put = (at: Address, key: string, entry: Entry) => {
    const k = `${printAddress(at)}\u0000${key}`;
    const had = entries.get(k);
    if (!had || had.entry.pos < entry.pos) entries.set(k, { at, key, entry });
    else if (had.entry.pos === entry.pos) had.entry.sites.push(...entry.sites);
  };
  let dangling = false;

  // the values and nodes the inverse restores, from the op that restores them
  const firstKey = new Set<string>();
  work.inverseMeta.forEach((m, i) => {
    const op = work.inverse[i];
    const pos = n - 1 - i;
    if (m?.kind === 'key') {
      const d = m.draft;
      if (d.made || !d.origin || !d.oldAddress) return;
      const k = `${printAddress(d.oldAddress)}\u0000${m.key}`;
      if (firstKey.has(k)) return; // the first op on a key restores the old root's value: the last in the inverse
      firstKey.add(k);
      if (op.op !== 'set' || (op.where ?? 'params') !== 'params') return;
      const sites: Entry['sites'] = [];
      const origin = d.origin;
      const site = (rel: TreePath, value: string, step?: Step) => {
        const run = live.runOf(value);
        if (run) sites.push({ value, run, step, expected: valueAt(origin, rel) });
      };
      oldPass.walkKey(origin.kind, m.key, op.value, site);
      loosePass.walkKey(origin.kind, m.key, op.value, (rel, value, step) => step?.loose && site(rel, value, step));
      put(d.oldAddress, m.key, { pos, sites });
    } else if (m?.kind === 'node') {
      if (op.op !== 'insert' && op.op !== 'replace') return;
      // the content the inverse restores, by lineage: each node of it the old root had, at its old address
      restored(work, m.draft, (node, origin, at) => {
        const site = (rel: TreePath, value: string, step?: Step) => {
          const run = live.runOf(value);
          if (run) put(at, String(rel[1]), { pos, sites: [{ value, run, step, expected: valueAt(origin, rel) }] });
        };
        oldPass.walkNode(node, site);
        loosePass.walkNode(node, (rel, value, step) => step?.loose && site(rel, value, step));
      });
    }
  });
  // the refs of the new root the inverse's moves may take along, from its first op
  if (moves.length) {
    const seen = new Map<Draft, Set<string>>();
    for (const m of work.moves) {
      for (const s of live.under(m.to)) {
        const h = s.holder;
        const key = String(s.rel[1]);
        if (h.made || !h.origin || !h.oldAddress || h.byOp?.has(key)) continue;
        let set = seen.get(h);
        if (!set) seen.set(h, (set = new Set()));
        const k = s.rel.join('\u0000');
        if (set.has(k)) continue;
        set.add(k);
        put(h.oldAddress, key, {
          pos: -1,
          sites: [{ value: s.value, run: s.run, step: s.step, expected: valueAt(h.origin, s.rel) }],
        });
      }
    }
    // a stored patch's `set` value is typed by the state the inverse's move finds, which the new root's typing does
    // not tell: every address in it is taken as one the inverse's moves may take along
    for (const { holder: h, key } of live.storedKeys()) {
      if (h.made || !h.origin || !h.oldAddress || h.byOp?.has(key)) continue;
      const value = h.params ? own(h.params, key) : undefined;
      if (value === undefined) continue;
      const origin = h.origin;
      const oldAddress = h.oldAddress;
      loosePass.walkKey(h.kind, key, value, (rel, v, step) => {
        const run = step?.loose && live.runOf(v);
        if (run) put(oldAddress, key, { pos: -1, sites: [{ value: v, run, step, expected: valueAt(origin, rel) }] });
      });
    }
  }

  const fixups: Op[] = [];
  for (const { at, key, entry } of entries.values()) {
    let off = false;
    for (const s of entry.sites) {
      let v = s.value;
      for (const m of moves) if (m.pos > entry.pos) v = rewriteRef(v, m.from, m.to);
      if (v !== s.expected) off = true;
    }
    // the inverse writes the key when an op of it restores the key or its node, or a set at its end does
    if (!off && entry.pos < 0) continue;
    const node = resolveAddress(before, at);
    const value = node?.params ? own(node.params, key) : undefined;
    if (value === undefined || !isNode(node)) continue;
    if (off) fixups.push({ op: 'set', at, key, value });
    // a ref it writes that dangled in the old root, typed as the old root types it (a stored patch's `set` value by
    // the old root's dry run, which the inverse's link integrity reads once the key is back)
    oldPass.walkKey(node.kind, key, value, (_rel, v, step) => {
      const run = live.runOf(v);
      if (run && !oldPass.names(run, step)) dangling = true;
    });
  }
  // content the inverse makes again that a ref of the old root named: the ref names a node the inverse made (its
  // identity changed, from the inverse's side) unless the inverse writes it
  let rebinds = false;
  for (const m of work.inverseMeta) {
    if (rebinds || m?.kind !== 'node') continue;
    // by lineage: a node moved into the content before it was taken out counts at its own old address
    restored(work, m.draft, (_node, _origin, at) => {
      if (!rebinds && live.namedInOld(at)) rebinds = true;
    });
  }
  return { fixups, dangling, rebinds };
}

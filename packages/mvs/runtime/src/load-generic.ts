/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { PluginContext } from '@molstar/plugin/context';
import { PluginState } from '@molstar/plugin/state';
import {
  State,
  StateObject,
  StateObjectSelector,
  StateTransform,
  StateTransformer,
  StateTree,
} from '@molstar/core/state';
import { UUID } from '@molstar/core/util';
import { stringHash } from '@molstar/mvs/helpers/utils';
import type { Kind, Subtree, SubtreeOfKind, Tree } from '@molstar/mvs-builder/tree/generic/tree-schema';

/** Function responsible for loading a tree node `node` into Mol*.
 * Should apply changes within `updateParent.update` but not commit them.
 * Should modify `context` accordingly, if it is needed for loading other nodes later.
 * `updateParent.selector` is the result of loading the node's parent into Mol* state hierarchy (or the hierarchy root in case of root node). */
export type LoadingAction<TNode extends Tree, TContext> = (
  updateParent: UpdateTarget,
  node: TNode,
  context: TContext,
) => UpdateTarget | undefined;

/** Loading actions for loading a tree into Mol*, per node kind. */
export type LoadingActions<TTree extends Tree, TContext> = {
  [kind in Kind<Subtree<TTree>>]?: LoadingAction<SubtreeOfKind<TTree, kind>, TContext>;
};

/** Type for defining custom behavior when loading trees, usually based on node custom properties. */
export interface LoadingExtension<TTree extends Tree, TContext, TExtensionContext> {
  id: string;
  description: string;
  /** Runs before the tree is loaded */
  createExtensionContext: (tree: TTree, context: TContext) => TExtensionContext;
  /** Runs after the tree is loaded */
  disposeExtensionContext?: (extensionContext: TExtensionContext, tree: TTree, context: TContext) => void;
  /**
   * The action reads only its node (kind, params, custom) and keeps nothing in the extension or loading context: what it
   * does to a subtree can be replayed when the subtree did not change (`load.md §4`). An extension without it turns
   * subtree reuse off.
   */
  local?: boolean;
  /** Runs on every node of the tree */
  action: (
    updateTarget: UpdateTarget,
    node: Subtree<TTree>,
    context: TContext,
    extensionContext: TExtensionContext,
  ) => void;
}

export interface LoadTreeVirtualOptions<TTree extends Tree, TContext> {
  /** The tree replaces the plugin's data tree (the snapshot holds the root and the loaded nodes only). */
  replaceExisting?: boolean;
  extensions?: LoadingExtension<TTree, TContext, any>[];
  /** Ref hashes kept between calls (`load.md §2`): one map for all the snapshots and frames of a load. */
  refHashCache?: Map<string, string>;
  /**
   * Serialized trees whose transforms may be reused (`load.md §3`), in order of preference (the previous frame, then
   * the snapshot): a transform equal to one of them (ref, transformer, parent, tags, dependsOn, state, and params as
   * JSON) is that object, version included, instead of a new one with a newly hashed version.
   */
  reuse?: StateTree.Serialized[];
  /** Replays the unit subtrees that did not change since an earlier tree (`load.md §4`). */
  subtreeMemo?: SubtreeMemo<TTree>;
}

/**
 * What loading each "unit" subtree did, to replay it when the next tree has the same subtree object (`load.md §4`). The
 * caller sets `units` before each load: the subtrees whose identity means "unchanged" (`convertMvsToMolstarMemo`), none
 * nested in another; `isReusable` excludes units whose loading reads outside them.
 */
export class SubtreeMemo<TTree extends Tree = Tree> {
  /** Counts over all memos, for tests and the bench: units replayed, units loaded, loads redone after a mismatch. */
  static readonly stats = { replayed: 0, loaded: 0, fallbacks: 0 };
  readonly records = new WeakMap<object, SubtreeRecord>();
  units: Set<object> = new Set();
  /** Loads replayed in full after a mismatch (`SubtreeReuseMismatch`): 0 when the rules of `load.md §4` hold. */
  fallbacks = 0;
  /** Recorded in the current load: their serialized transforms are known after `serialize`. */
  readonly pending: SubtreeRecord[] = [];
  constructor(readonly isReusable: (node: Subtree<TTree>) => boolean) {}
}

export interface SubtreeRecord {
  /** The target the subtree was loaded under. */
  parentRef: StateTransform.Ref;
  /** The transforms the subtree added, in order, and what they serialized to. */
  applies: { parent: StateTransform.Ref; transformer: StateTransformer; ref: StateTransform.Ref }[];
  serialized?: StateTransform.Serialized[];
  /** Ref counters taken: hash, number, hash, number, ... */
  hashes: (string | number)[];
  /** MVS ref → transform ref */
  refMap: [string, string][];
  created: Set<StateTransform.Ref>;
  /** False when the loading of the subtree touched transforms it did not add. */
  reusable: boolean;
}

/** A replayed subtree does not fit the tree it was replayed into; the load is done again without reuse. */
export class SubtreeReuseMismatch extends Error {}

export function loadTreeVirtual<TTree extends Tree, TContext>(
  plugin: PluginContext,
  tree: TTree,
  loadingActions: LoadingActions<TTree, TContext>,
  context: TContext,
  options?: LoadTreeVirtualOptions<TTree, TContext>,
) {
  const updateRoot: UpdateTarget = UpdateTarget.create(
    plugin,
    options?.replaceExisting ?? false,
    options?.refHashCache,
  );
  const memo = options?.subtreeMemo;
  if (memo && !(options?.extensions ?? []).every((e) => e.local)) {
    options = { ...options, subtreeMemo: undefined };
  }
  let serialized: StateTree.Serialized;
  try {
    loadTreeInUpdate(updateRoot, tree, loadingActions, context, options);
    serialized = updateRoot.builder.serialize(options?.reuse);
    if (memo) finishRecords(memo, serialized);
  } catch (e) {
    if (!(e instanceof SubtreeReuseMismatch) || !memo) throw e;
    memo.fallbacks++;
    SubtreeMemo.stats.fallbacks++;
    memo.pending.length = 0;
    return loadTreeVirtual(plugin, tree, loadingActions, context, { ...options, subtreeMemo: undefined });
  }
  const stateSnapshot: State.Snapshot = { tree: serialized };
  const pluginStateSnapshot: PluginState.Snapshot = { id: UUID.create22(), data: stateSnapshot };
  return pluginStateSnapshot;
}

function loadTreeInUpdate<TTree extends Tree, TContext>(
  updateRoot: UpdateTarget,
  tree: TTree,
  loadingActions: LoadingActions<TTree, TContext>,
  context: TContext,
  options?: LoadTreeVirtualOptions<TTree, TContext>,
) {
  const mapping = new Map<Subtree<TTree>, UpdateTarget | undefined>();
  if (options?.replaceExisting) {
    UpdateTarget.deleteChildren(updateRoot);
  }
  const extensionContexts = (options?.extensions ?? []).map((ext) => ({
    ext,
    extCtx: ext.createExtensionContext(tree, context),
  }));
  const mvsRefMap = new Map<string, string>();
  const { builder, targetManager } = updateRoot;
  const visitNode = (node: Subtree<TTree>, parent: Subtree<TTree> | undefined) => {
    const kind: Kind<typeof node> = node.kind;
    let msNode: UpdateTarget | undefined;
    const updateParent = parent ? mapping.get(parent) : updateRoot;
    const action = loadingActions[kind] as LoadingAction<typeof node, TContext> | undefined;
    if (action) {
      if (updateParent) {
        msNode = action(updateParent, node, context);
        if (msNode && node.ref) {
          UpdateTarget.tag(msNode, mvsRefTags(node.ref));
          mvsRefMap.set(node.ref, msNode.selector.ref);
          builder.recording?.refMap.push([node.ref, msNode.selector.ref]);
        }
        mapping.set(node, msNode);
      } else {
        console.warn(`No target found for this "${node.kind}" node`);
        return;
      }
    }
    if (updateParent) {
      for (const { ext, extCtx } of extensionContexts) {
        ext.action(msNode ?? updateParent, node, context, extCtx);
      }
    }
  };

  const memo = options?.subtreeMemo;
  const visit = (node: Subtree<TTree>, parent: Subtree<TTree> | undefined) => {
    const updateParent = memo?.units.has(node) && parent ? mapping.get(parent) : undefined;
    if (!updateParent) {
      visitNode(node, parent);
      for (const child of node.children ?? []) visit(child as Subtree<TTree>, node);
      return;
    }
    // a unit (`load.md §4`): replay what loading it did last time, or load it and record that
    const record = memo!.records.get(node);
    if (record?.reusable && record.serialized && record.parentRef === updateParent.selector.ref) {
      for (let i = 0; i < record.hashes.length; i += 2) {
        targetManager.takeRef(record.hashes[i] as string, record.hashes[i + 1] as number);
      }
      for (let i = 0; i < record.applies.length; i++) {
        const { parent: p, transformer, ref } = record.applies[i];
        builder.applySerialized(p, transformer, ref, record.serialized[i]);
      }
      for (const [mvsRef, ref] of record.refMap) mvsRefMap.set(mvsRef, ref);
      SubtreeMemo.stats.replayed++;
      return;
    }
    const recording: SubtreeRecord = {
      parentRef: updateParent.selector.ref,
      applies: [],
      hashes: [],
      refMap: [],
      created: new Set(),
      reusable: memo!.isReusable(node),
    };
    builder.recording = recording;
    targetManager.recording = recording;
    try {
      visitNode(node, parent);
      for (const child of node.children ?? []) visit(child as Subtree<TTree>, node);
    } finally {
      builder.recording = undefined;
      targetManager.recording = undefined;
    }
    memo!.records.set(node, recording);
    memo!.pending.push(recording);
    SubtreeMemo.stats.loaded++;
  };
  visit(tree as Subtree<TTree>, undefined);

  for (const target of updateRoot.targetManager.allTargets) {
    UpdateTarget.dependsOn(target, mvsRefMap);
  }

  extensionContexts.forEach((e) => e.ext.disposeExtensionContext?.(e.extCtx, tree, context));
}

/** The serialized transforms of the subtrees recorded in this load. */
function finishRecords(memo: SubtreeMemo<any>, tree: StateTree.Serialized) {
  if (!memo.pending.length) return;
  const byRef = new Map(tree.transforms.map((t) => [t.ref, t]));
  for (const record of memo.pending) {
    const serialized = record.applies.map((a) => byRef.get(a.ref)!);
    if (serialized.every((s) => s)) record.serialized = serialized;
    else record.reusable = false;
  }
  memo.pending.length = 0;
}

/** A wrapper for building a Mol* state tree, with deterministic transform refs (`load.md §2`).
 * ```
 * updateTarget = UpdateTarget.create(plugin, true); // like update = plugin.build();
 * UpdateTarget.apply(updateTarget, transformer, params); // like update.to(selector).apply(transformer, params);
 * updateTarget.builder.serialize(); // like StateTree.toJSON(update.getTree({ useHashVersion: true }))
 * ```
 */
export interface UpdateTarget {
  readonly builder: VirtualStateBuilder;
  readonly selector: StateObjectSelector;
  readonly targetManager: TargetManager;
  readonly mvsDependencyRefs: Set<string>;

  readonly transformer?: StateTransformer;
  readonly transformParams?: any;
}
export const UpdateTarget = {
  /** Create a new update, with `selector` pointing to the root. */
  create(plugin: PluginContext, replaceExisting: boolean, refHashCache?: Map<string, string>): UpdateTarget {
    const builder = new VirtualStateBuilder(plugin.state.data, replaceExisting);
    return {
      builder,
      selector: new StateObjectSelector(builder.rootRef, plugin.state.data),
      targetManager: new TargetManager(plugin, replaceExisting, refHashCache),
      mvsDependencyRefs: new Set(),
    };
  },
  /** Add a child node to `target.selector`, return a new `UpdateTarget` pointing to the new child. */
  apply<A extends StateObject, B extends StateObject, P extends {}>(
    target: UpdateTarget,
    transformer: StateTransformer<A, B, P>,
    params?: Partial<P>,
    options?: Partial<StateTransform.Options>,
  ): UpdateTarget {
    let refSuffix: string = transformer.id;
    if (transformer.id === StructureRepresentation3D.id) {
      const reprType = (params as any)?.type?.name ?? '';
      refSuffix += `:${reprType}`;
    }
    const ref = target.targetManager.getChildRef(target.selector, refSuffix);
    target.builder.apply(target.selector.ref, transformer, params, { ...options, ref });
    const result: UpdateTarget = {
      ...target,
      selector: new StateObjectSelector(ref, target.selector.state),
      mvsDependencyRefs: new Set(),
      transformer,
      transformParams: params,
    };
    target.targetManager.allTargets.push(result);
    return result;
  },
  setMvsDependencies(target: UpdateTarget, refs: string[] | Set<string>): UpdateTarget {
    refs.forEach((ref) => target.mvsDependencyRefs.add(ref));
    return target;
  },
  dependsOn(target: UpdateTarget, mapping: Map<string, string>): UpdateTarget {
    if (!target.mvsDependencyRefs.size) return target;
    const dependsOn = Array.from(target.mvsDependencyRefs)
      .map((d) => mapping.get(d)!)
      .filter((d) => d);
    if (!dependsOn.length) return target;
    target.builder.dependsOn(target.selector.ref, dependsOn);
    return target;
  },
  /** Add tags to `target.selector` */
  tag(target: UpdateTarget, tags: string[]): UpdateTarget {
    if (tags.length > 0) {
      target.builder.tag(target.selector.ref, tags);
    }
    return target;
  },
  /** Update the state (e.g. `isHidden`) of `target.selector` */
  updateState(target: UpdateTarget, state: Partial<StateTransform.State>): UpdateTarget {
    target.builder.updateState(target.selector.ref, state);
    return target;
  },
  /** Delete all children of `target.selector`. */
  deleteChildren(target: UpdateTarget): UpdateTarget {
    target.builder.deleteChildren(target.selector.ref);
    return target;
  },
};

interface VirtualNode {
  transform: StateTransform;
  children: StateTransform.Ref[];
  /** Replayed (`load.md §4`): what it serializes to. */
  serialized?: StateTransform.Serialized;
}

/**
 * The subset of `StateBuilder` the MVS loader uses, building a serialized tree directly (`load.md §3`): no immutable
 * maps, no state, nothing to commit. Same semantics, so the same output as `StateTree.toJSON(builder.getTree({
 * useHashVersion: true }))`: decorators are inserted between a node and its children, other transformers are applied
 * under the node's decorator root, tags and dependsOn are sorted unions, the transforms added or changed get the hash
 * of their params as version, and the order is the tree's pre-order with children in insertion order.
 */
export class VirtualStateBuilder {
  private readonly nodes = new Map<StateTransform.Ref, VirtualNode>();
  /** Added or changed: their version is the hash of the params. */
  private readonly affected = new Set<StateTransform.Ref>();
  readonly rootRef: StateTransform.Ref;
  /** The unit subtree being loaded, if any (`load.md §4`). */
  recording: SubtreeRecord | undefined = undefined;

  constructor(state: State, replaceExisting: boolean) {
    const tree = state.tree;
    this.rootRef = tree.root.ref;
    this.nodes.set(tree.root.ref, { transform: tree.root, children: [] });
    if (!replaceExisting) {
      StateTree.doPreOrder(tree, tree.root, this.nodes, (t, _, nodes) => {
        if (t.ref === tree.root.ref) return;
        nodes.set(t.ref, { transform: t, children: [] });
        nodes.get(t.parent)!.children.push(t.ref);
      });
    }
  }

  private node(ref: StateTransform.Ref) {
    const node = this.nodes.get(ref);
    if (!node) throw new Error(`Could not find node '${ref}'.`);
    return node;
  }

  /** A change of an existing transform: not replayable when the unit being recorded did not add it. */
  private changed(ref: StateTransform.Ref) {
    if (this.nodes.get(ref)?.serialized) throw new SubtreeReuseMismatch(`Replayed transform '${ref}' changed.`);
    if (this.recording && !this.recording.created.has(ref)) this.recording.reusable = false;
  }

  private decoratorRoot(ref: StateTransform.Ref): StateTransform.Ref {
    const { children } = this.node(ref);
    if (children.length !== 1) return ref;
    const child = this.nodes.get(children[0]);
    if (child?.transform.transformer.definition.isDecorator) return this.decoratorRoot(child.transform.ref);
    return ref;
  }

  private add(t: StateTransform) {
    const existing = this.nodes.get(t.ref);
    if (existing) {
      if (existing.transform.parent !== t.parent) throw new Error(`Transform '${t.ref}' is already present.`);
      existing.transform = t;
    } else {
      this.node(t.parent).children.push(t.ref);
      this.nodes.set(t.ref, { transform: t, children: [] });
    }
    this.affected.add(t.ref);
  }

  /** `StateBuilder.To.apply` on `parent`; returns the ref of the new transform. */
  apply(
    parent: StateTransform.Ref,
    transformer: StateTransformer,
    params: any,
    options: Partial<StateTransform.Options>,
  ): StateTransform.Ref {
    if (transformer.definition.isDecorator && this.nodes.get(parent)?.serialized) {
      throw new SubtreeReuseMismatch(`Replayed transform '${parent}' changed.`);
    }
    const ref = this.place(parent, transformer, params, options);
    if (this.recording) {
      this.recording.applies.push({ parent, transformer, ref });
      this.recording.created.add(ref);
    }
    return ref;
  }

  /** `apply` of a transform recorded earlier, whose serialized form is known (`load.md §4`). */
  applySerialized(
    parent: StateTransform.Ref,
    transformer: StateTransformer,
    ref: StateTransform.Ref,
    serialized: StateTransform.Serialized,
  ) {
    this.place(parent, transformer, undefined, { ref });
    this.node(ref).serialized = serialized;
    this.affected.delete(ref);
  }

  private place(
    parent: StateTransform.Ref,
    transformer: StateTransformer,
    params: any,
    options: Partial<StateTransform.Options>,
  ): StateTransform.Ref {
    if (transformer.definition.isDecorator) {
      // StateBuilder.To.insert: the new transform takes over the children of `parent`
      if (this.recording && !this.recording.created.has(parent)) this.recording.reusable = false;
      const node = this.node(parent);
      const children = node.children;
      const t = transformer.apply(parent, params, options);
      node.children = [];
      this.add(t);
      const inserted = this.node(t.ref);
      for (const c of children) {
        const child = this.node(c);
        child.transform = StateTransform.withParent(child.transform, t.ref);
        inserted.children.push(c);
      }
      this.affected.add(parent);
      return t.ref;
    }
    const t = transformer.apply(this.decoratorRoot(parent), params, options);
    this.add(t);
    return t.ref;
  }

  tag(ref: StateTransform.Ref, tags: string[]) {
    const node = this.node(ref);
    const updated = StateTransform.withTags(node.transform, stringArrayUnion(node.transform.tags, tags));
    if (arrayEqual(node.transform.tags, updated.tags)) return;
    this.changed(ref);
    node.transform = updated;
    this.affected.add(ref);
  }

  dependsOn(ref: StateTransform.Ref, dependsOn: string[]) {
    const node = this.node(ref);
    const updated = StateTransform.withDependsOn(node.transform, stringArrayUnion(node.transform.dependsOn, dependsOn));
    if (arrayEqual(node.transform.dependsOn, updated.dependsOn)) return;
    this.changed(ref);
    node.transform = updated;
    this.affected.add(ref);
  }

  updateState(ref: StateTransform.Ref, state: Partial<StateTransform.State>) {
    const node = this.node(ref);
    if (!StateTransform.isStateChange(node.transform.state, state)) return;
    this.changed(ref);
    node.transform = StateTransform.withState(node.transform, state);
    this.affected.add(ref);
  }

  deleteChildren(ref: StateTransform.Ref) {
    const node = this.node(ref);
    const remove = (r: StateTransform.Ref) => {
      for (const c of this.node(r).children) remove(c);
      this.nodes.delete(r);
      this.affected.delete(r);
    };
    for (const c of node.children) remove(c);
    node.children = [];
  }

  /** The serialized tree in pre-order; transforms equal to one in `reuse` are that one (`load.md §3`). */
  serialize(reuse?: StateTree.Serialized[]): StateTree.Serialized {
    const transforms: StateTransform.Serialized[] = [];
    const reuseMaps: (Map<string, StateTransform.Serialized> | undefined)[] = reuse ? reuse.map(() => undefined) : [];
    const visit = (ref: StateTransform.Ref) => {
      const { transform, children, serialized } = this.node(ref);
      if (serialized) {
        if (serialized.parent !== transform.parent) {
          throw new SubtreeReuseMismatch(`Replayed transform '${ref}' is under another parent.`);
        }
        transforms.push(serialized);
      } else if (!this.affected.has(ref)) {
        transforms.push(StateTransform.toJSON(transform));
      } else {
        let found: StateTransform.Serialized | undefined;
        let closest: StateTransform.Serialized | undefined;
        let shared: any;
        if (reuse) {
          const index = transforms.length;
          for (let i = 0; i < reuse.length && !found; i++) {
            let candidate: StateTransform.Serialized | undefined = reuse[i].transforms[index];
            if (candidate?.ref !== ref) {
              reuseMaps[i] ??= new Map(reuse[i].transforms.map((t) => [t.ref, t]));
              candidate = reuseMaps[i]!.get(ref);
            }
            if (!candidate) continue;
            if (!closest) {
              closest = candidate;
              if (canShareParams(transform, candidate)) {
                // one walk: equal → the candidate; else the params with its equal parts (memory of long animations)
                shared = shareEqualParts(transform.params, candidate.params);
                if (shared === candidate.params && isSameTransformMeta(transform, candidate)) found = candidate;
                continue;
              }
            }
            if (isSameTransform(transform, candidate)) found = candidate;
          }
        }
        if (found) {
          transforms.push(found);
        } else {
          let t = transform;
          if (shared !== undefined && shared !== t.params) t = { ...t, params: shared };
          if (closest) {
            if (t.tags && arrayEqual(t.tags, closest.tags)) t = { ...t, tags: closest.tags };
            if (t.dependsOn && arrayEqual(t.dependsOn, closest.dependsOn)) t = { ...t, dependsOn: closest.dependsOn };
          }
          StateTransform.setParamsHashVersion(t);
          transforms.push(StateTransform.toJSON(t));
        }
      }
      for (const c of children) visit(c);
    };
    visit(this.rootRef);
    return { transforms };
  }
}

/** `t` would serialize to `s` with the same hash version (`load.md §3`). */
function isSameTransform(t: StateTransform, s: StateTransform.Serialized) {
  if (!isSameTransformMeta(t, s)) return false;
  // the version is the hash of `JSON.stringify(t.params ?? {})`; `s.params` is `t.params`, or undefined for none
  if ((t.params === undefined || t.params === null) !== (s.params === undefined || s.params === null)) return false;
  return jsonEqual(t.params ?? {}, s.params ?? {});
}

/** All but the params of `isSameTransform`. */
function isSameTransformMeta(t: StateTransform, s: StateTransform.Serialized) {
  if (t.transformer.id !== s.transformer || t.parent !== s.parent) return false;
  // (params serialized by the transformer: the hash is of the raw params, which `s` does not keep)
  if (t.transformer.definition.customSerialization) return false;
  if (!arrayEqual(t.tags, s.tags) || !arrayEqual(t.dependsOn, s.dependsOn)) return false;
  return isSameSerializedState(t.state, s.state);
}

/** `s.params` is what `t.params` would serialize to (no custom serialization), and both are objects. */
function canShareParams(t: StateTransform, s: StateTransform.Serialized) {
  return (
    t.transformer.id === s.transformer &&
    !t.transformer.definition.customSerialization &&
    typeof t.params === 'object' &&
    t.params !== null &&
    typeof s.params === 'object' &&
    s.params !== null
  );
}

function isSameSerializedState(state: StateTransform.State, serialized: StateTransform.Serialized['state']) {
  let n = 0;
  for (const k of Object.keys(state)) {
    if (!(state as any)[k]) continue;
    n++;
    if (!serialized || !(serialized as any)[k]) return false;
  }
  return n === (serialized ? Object.keys(serialized).length : 0);
}

/** `JSON.stringify(a) === JSON.stringify(b)`, without building the strings. */
export function jsonEqual(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  if (typeof (a as any)?.toJSON === 'function' || typeof (b as any)?.toJSON === 'function') {
    return JSON.stringify(a) === JSON.stringify(b);
  }
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null) {
    return jsonPrimitive(a) === jsonPrimitive(b);
  }
  const arrA = Array.isArray(a);
  if (arrA !== Array.isArray(b)) return false;
  if (arrA) {
    const x = a as unknown[];
    const y = b as unknown[];
    if (x.length !== y.length) return false;
    for (let i = 0; i < x.length; i++) {
      // (in arrays, values JSON leaves out are null)
      const u = isJsonOmitted(x[i]) ? null : x[i];
      const v = isJsonOmitted(y[i]) ? null : y[i];
      if (!jsonEqual(u, v)) return false;
    }
    return true;
  }
  const keysA = Object.keys(a);
  const keysB = Object.keys(b);
  let i = 0;
  let j = 0;
  while (true) {
    while (i < keysA.length && isJsonOmitted((a as any)[keysA[i]])) i++;
    while (j < keysB.length && isJsonOmitted((b as any)[keysB[j]])) j++;
    if (i === keysA.length || j === keysB.length) return i === keysA.length && j === keysB.length;
    if (keysA[i] !== keysB[j] || !jsonEqual((a as any)[keysA[i]], (b as any)[keysB[j]])) return false;
    i++;
    j++;
  }
}

/**
 * `fresh` with each part that is `jsonEqual` to the same part of `old` replaced by that part of `old`; plain objects and
 * arrays on the way are copied, never written. Same JSON as `fresh`.
 */
export function shareEqualParts(fresh: any, old: any): any {
  if (fresh === old) return old;
  if (!isPlainContainer(fresh) || !isPlainContainer(old) || Array.isArray(fresh) !== Array.isArray(old)) {
    return jsonEqual(fresh, old) ? old : fresh;
  }
  // one walk: `equal` stays true while every part is `old`'s (then `fresh` is JSON-equal to `old`)
  let copy: any;
  let equal: boolean;
  if (Array.isArray(fresh)) {
    equal = fresh.length === old.length;
    for (let i = 0; i < fresh.length && i < old.length; i++) {
      // (in arrays, values JSON leaves out are null: compared as such)
      const f = isJsonOmitted(fresh[i]) ? null : fresh[i];
      const o = isJsonOmitted(old[i]) ? null : old[i];
      const v = shareEqualParts(f, o);
      if (v !== o) equal = false;
      if (v !== f) (copy ??= fresh.slice())[i] = v;
    }
  } else {
    const keysF = Object.keys(fresh);
    const keysO = Object.keys(old);
    let i = 0;
    let j = 0;
    equal = true;
    while (i < keysF.length || j < keysO.length) {
      while (i < keysF.length && isJsonOmitted(fresh[keysF[i]])) i++;
      while (j < keysO.length && isJsonOmitted(old[keysO[j]])) j++;
      if (i === keysF.length && j === keysO.length) break;
      if (i === keysF.length || j === keysO.length || keysF[i] !== keysO[j]) {
        equal = false;
        break;
      }
      i++;
      j++;
    }
    for (const k of keysF) {
      const f = fresh[k];
      if (isJsonOmitted(f) || !(k in old)) continue;
      const v = shareEqualParts(f, old[k]);
      if (v !== old[k]) equal = false;
      if (v !== f) (copy ??= { ...fresh })[k] = v;
    }
  }
  return equal ? old : (copy ?? fresh);
}

function isPlainContainer(x: unknown): boolean {
  if (typeof x !== 'object' || x === null) return false;
  if (Array.isArray(x)) return true;
  const proto = Object.getPrototypeOf(x);
  return proto === Object.prototype || proto === null;
}

function isJsonOmitted(x: unknown) {
  return x === undefined || typeof x === 'function' || typeof x === 'symbol';
}

function jsonPrimitive(x: unknown) {
  if (typeof x === 'number') return Number.isFinite(x) ? x : null;
  if (isJsonOmitted(x)) return undefined;
  return x;
}

function arrayEqual(a: string[] | undefined, b: string[] | undefined) {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** As `StateBuilder`'s: the union in order of first appearance. */
function stringArrayUnion(...arrays: (string[] | undefined)[]): string[] {
  const set = new Set<string>();
  const ret: string[] = [];
  for (const xs of arrays) {
    if (!xs) continue;
    for (const x of xs) {
      if (set.has(x)) continue;
      set.add(x);
      ret.push(x);
    }
  }
  return ret;
}

/** Manages transform refs in a deterministic way. Uses refs like !mvs:3ce3664304d32c5d:0 */
class TargetManager {
  /** For each hash (e.g. 3ce3664304d32c5d), store the number of already used refs with that hash. */
  private _counter: Record<string, number> = {};
  /** `hashBase → hash`; shared between the calls of a load (refs are the same in every frame) */
  private readonly hashCache: Map<string, string>;
  constructor(plugin: PluginContext, replaceExisting: boolean, hashCache?: Map<string, string>) {
    this.hashCache = hashCache ?? new Map();
    if (!replaceExisting) {
      plugin.state.data.cells.forEach((cell) => {
        const ref = cell.transform.ref;
        if (ref.startsWith('!mvs:')) {
          const [_, hash, idNumber] = ref.split(':');
          const nextIdNumber = parseInt(idNumber) + 1;
          if (nextIdNumber > (this._counter[hash] ?? 0)) {
            this._counter[hash] = nextIdNumber;
          }
        }
      });
    }
  }
  /** The unit subtree being loaded, if any (`load.md §4`). */
  recording: SubtreeRecord | undefined = undefined;

  /** Return ref for a new node with given `hash`; update the counter accordingly. */
  private nextRef(hash: string): string {
    this._counter[hash] ??= 0;
    const idNumber = this._counter[hash]++;
    this.recording?.hashes.push(hash, idNumber);
    return `!mvs:${hash}:${idNumber}`;
  }

  /** Take the ref `!mvs:<hash>:<n>` again in a replay; it must be the next one (`load.md §4`). */
  takeRef(hash: string, n: number) {
    if ((this._counter[hash] ?? 0) !== n) throw new SubtreeReuseMismatch(`Ref counter of ${hash} is not ${n}.`);
    this._counter[hash] = n + 1;
  }
  /** Return ref for a new node based on parent and desired suffix. */
  getChildRef(parent: StateObjectSelector, suffix: string): string {
    const hashBase = parent.ref.replace(/^!mvs:/, '') + ':' + suffix;
    let hash = this.hashCache.get(hashBase);
    if (hash === undefined) {
      hash = stringHash(hashBase);
      this.hashCache.set(hashBase, hash);
    }
    const result = this.nextRef(hash);
    return result;
  }

  readonly allTargets: UpdateTarget[] = [];
}

/** Create node tags based of MVS node.ref */
export function mvsRefTags(mvsNodeRef: string | undefined): string[] {
  if (mvsNodeRef === undefined) return [];
  else return [`mvs-ref:${mvsNodeRef}`];
}

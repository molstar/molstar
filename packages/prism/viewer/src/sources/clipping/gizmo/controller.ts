/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GizmoController } from '@molstar/graphics/util/gizmo/controller';
import type { PluginContext } from '@molstar/plugin/context';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { BehaviorSubject } from 'rxjs';
import { type Cut, type CutObject, representationsSphere } from '../cuts.js';
import { ClipBinding, type ClipItem, type GizmoCommit } from './binding.js';

export type { GizmoCommit } from './binding.js';

// The clipping gizmos (clipping.md §4): the gizmo library's controller (`@molstar/graphics/util/gizmo`) on the cuts'
// binding (`ClipBinding`), for the active cut (its editor open: `handleOnly`), shown while the section is (`show`,
// `hide`), its clicks claimed (`InteractivityManager.clickOwners`: no focus, no select). The active cut outlives `hide`: shown
// again, its gizmo comes back. A cut is known by its number, which shifts when one before it goes: the active cut is
// followed by its object (its tag), or by its number when its object changed and the count did not (an edit).

/**
 * The active cuts `before` (by number) in `cuts`: each by its object (its tag), else by its number when its object is
 * gone and the count is the same as `followedCount` (it was edited; a removal shifts the numbers, and the next cut is
 * not the removed one)
 */
export function followCuts(
  before: ReadonlyMap<number, CutObject>,
  cuts: readonly Pick<Cut, 'index' | 'object'>[],
  followedCount: number,
): Map<number, CutObject> {
  const next = new Map<number, CutObject>();
  const sameCount = cuts.length === followedCount;
  for (const cut of cuts) {
    const same = [...before.values()].some((o) => RepresentationClip.sameCut(o, cut.object));
    const old = before.get(cut.index);
    const edited = sameCount && !!old && !cuts.some((c) => RepresentationClip.sameCut(c.object, old));
    if (same || edited) next.set(cut.index, cut.object);
  }
  return next;
}

export class ClipGizmos {
  readonly #shown = observableState(false);
  /** The active cuts by number (one: the section's), and their objects as last seen (to follow them) */
  readonly #handled = observableState<readonly number[]>([]);
  #handledObjects = new Map<number, CutObject>();
  /** How many cuts there were when the active ones were last followed: an edit keeps the count, a removal does not */
  #followedCount = 0;
  readonly #fadeNormal = observableState(false);
  readonly #items = new BehaviorSubject<readonly ClipItem[]>([]);
  #controller: GizmoController<CutObject> | undefined;
  #binding: ClipBinding | undefined;
  #subs: { unsubscribe(): void }[] = [];
  #releaseClicks: (() => void) | undefined;

  constructor(
    private readonly plugin: PluginContext,
    private readonly cuts: ReadonlyState<readonly Cut[]>,
  ) {}

  get shown(): ReadonlyState<boolean> {
    return this.#shown.asReadonly;
  }

  /** The cuts whose gizmos show (while `shown`), by number */
  get handled(): ReadonlyState<readonly number[]> {
    return this.#handled.asReadonly;
  }

  /** Whether the slide arrow fades out as it turns toward the viewer; by default always drawn (clipping.md §4) */
  get fadeNormal(): ReadonlyState<boolean> {
    return this.#fadeNormal.asReadonly;
  }

  setFadeNormal(fade: boolean): void {
    this.#fadeNormal.set(fade);
    this.#controller?.setFadeSlides(fade);
  }

  /** The active cut: its gizmo alone, or none */
  handleOnly(index: number | undefined): void {
    const cut = index === undefined ? undefined : this.cuts.value.find((c) => c.index === index);
    this.#handledObjects = new Map(cut ? [[cut.index, cut.object] as const] : []);
    this.#followedCount = this.cuts.value.length;
    this.#publish(this.cuts.value);
  }

  /** Shows the gizmos; a drag's release commits through `commit`. Without a canvas (headless): false */
  show(commit: GizmoCommit): boolean {
    if (this.#binding) {
      this.#binding.setCommit(commit);
      return true;
    }
    const canvas = this.plugin.canvas3d;
    if (!canvas) return false;
    const binding = new ClipBinding(this.plugin, this.#items);
    binding.setCommit(commit);
    const controller = new GizmoController(canvas, binding, { fadeSlides: this.#fadeNormal.value });
    this.#binding = binding;
    this.#controller = controller;
    this.#releaseClicks = this.plugin.managers.interactivity.clickOwners.claim((r) => controller.owns(r));
    this.#subs.push(this.cuts.subscribe((cuts) => this.#publish(cuts)));
    this.#shown.set(true);
    return true;
  }

  hide(): void {
    if (!this.#controller) return;
    for (const s of this.#subs) s.unsubscribe();
    this.#subs = [];
    this.#releaseClicks?.();
    this.#releaseClicks = undefined;
    this.#controller.dispose();
    this.#binding?.dispose();
    this.#controller = undefined;
    this.#binding = undefined;
    this.#shown.set(false);
  }

  /** The drag in progress, undone */
  cancel(): void {
    this.#controller?.cancel();
  }

  /** The active cuts followed into `cuts`, and published as the gizmos' items */
  #publish(cuts: readonly Cut[]) {
    // during a drag the cut is the drag's (its number and object are followed again once it ends)
    if (this.#controller?.dragging.value === undefined) this.#follow(cuts);
    const items: ClipItem[] = [];
    for (const cut of cuts) {
      if (!this.#handledObjects.has(cut.index)) continue;
      const bounds = representationsSphere(
        cut.targets.flatMap((t) => {
          const cell = this.plugin.state.data.cells.get(t.ref);
          return cell ? [cell] : [];
        }),
      );
      if (!bounds) continue;
      items.push({ key: cut.object.tag || `#${cut.index}`, value: cut.object, label: `Cut ${cut.index}`, bounds });
    }
    this.#items.next(items);
  }

  /** The active cuts followed into `cuts` (`followCuts`) */
  #follow(cuts: readonly Cut[]) {
    const next = followCuts(this.#handledObjects, cuts, this.#followedCount);
    this.#followedCount = cuts.length;
    this.#handledObjects = next;
    const ids = [...next.keys()];
    const prev = this.#handled.value;
    if (ids.length !== prev.length || ids.some((x, i) => x !== prev[i])) this.#handled.set(ids);
  }
}

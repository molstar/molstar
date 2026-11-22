/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { asyncScheduler, Subject } from 'rxjs';
import { debounceTime, throttleTime } from 'rxjs/operators';
import type { ParticleTileLoader } from '@molstar/model/model/particles/particle-tile-loader';
import {
  type State,
  type StateObject,
  StateSelection,
  type StateObjectCell,
  type StateTransformer,
} from '@molstar/core/state';

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginBehavior } from '../behavior.js';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import {
  isParticleTilesViewBehind,
  setParticleTilesCameraChanged,
  ParticleTilesView,
  setParticleTilesStep,
} from '@molstar/plugin/state/transforms/particles/tiles';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ParticleTileDrawings } from '@molstar/model/model/particles/particle-tile-drawings';

import type { ParticleTileSet } from '@molstar/model/model/particles/particle-tiles';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';

// Particle tile views by the camera that follow it (`ParticleTilesView`'s camera mode, `follow` on, no camera given):
// as the canvas' camera moves (at most every `intervalMs`) and as tiles load, a view behind what the camera picks
// (`isParticleTilesViewBehind`) updated a step at a time (its `budget` of particles added a step), a frame drawn
// between steps: the tiles streamed in, the frames not blocked. A view fixed to a camera given (a snapshot) is not
// followed: tests and comparisons need one.

type TilesViewParams = StateTransformer.Params<ParticleTilesView>;

/** whether a view's params follow the camera */
function isFollowing(params: TilesViewParams | undefined) {
  const tiles = params?.tiles;
  return tiles?.name === 'camera' && !!tiles.params.follow && !tiles.params.camera;
}

/** a view's particles (a slot list's live ones) */
const particlesOf = (list: ParticleList | undefined) => (list ? ParticleSlots.getLiveCount(list) : 0);

/**
 * the particles of a step that takes about `stepMs` at `rate` particles a ms (the last steps'), whole: none where no
 * rate is known yet (the budget's)
 */
export function getParticleTilesStepParticles(rate: number | undefined, stepMs: number) {
  if (rate === undefined || !(rate > 0) || !(stepMs > 0)) return undefined;
  return Math.max(1, Math.round(rate * stepMs));
}

/**
 * A camera's changes that move what it sees: its view or its projection but for the depth range (near and far, the
 * projection's 10th and 14th entries, perspective or orthographic), which follows the scene's bounds (a view's own
 * tiles as they load) and which a pick does not read. Why: a change read as a move picks the quicker way (fewer rays)
 * and in full once settled; where those differ on a tile, each step it brings moves the bounds again: a step every
 * few frames with the camera still, without end (a sample's jitter did the same: `Camera.update`'s silent)
 */
export function createParticleTilesCameraMoved() {
  const view = new Float64Array(16).fill(Number.NaN);
  const projection = new Float64Array(16).fill(Number.NaN);
  return (camera: { readonly view: ArrayLike<number>; readonly projection: ArrayLike<number> }) => {
    let moved = false;
    for (let i = 0; i < 16; ++i) {
      if (view[i] !== camera.view[i]) moved = true;
      view[i] = camera.view[i];
      if (i === 10 || i === 14) continue;
      if (projection[i] !== camera.projection[i]) moved = true;
      projection[i] = camera.projection[i];
    }
    return moved;
  };
}

/** the least time between passes while what draws a set changes (a cut dragged), in ms */
const DrawingsIntervalMs = 250;

/** the next frame (a step's update drawn before the next), a timeout without one (headless) */
const nextFrame = (f: () => void) =>
  typeof requestAnimationFrame !== 'undefined' ? requestAnimationFrame(() => f()) : setTimeout(f, 0);

export const ParticleTilesFollowCamera = PluginBehavior.create<{ intervalMs: number }>({
  name: 'particle-tiles-follow-camera',
  category: 'interaction',
  display: {
    name: 'Particle Tiles Follow the Camera',
    description: 'Tile views by the camera pick their tiles again as it moves, adding them in steps as they load.',
  },
  params: () => ({
    intervalMs: PD.Numeric(
      100,
      { min: 0, max: 2000, step: 10 },
      { description: 'How often the tiles are picked again while the camera moves, ms.' },
    ),
  }),
  ctor: class extends PluginBehavior.Handler<{ intervalMs: number }> {
    /** the tile sets' loaders followed: a tile loaded wakes the views of its set; a set's let go as it goes */
    private loaders = new Map<ParticleTileLoader, () => void>();
    /** the tile sets' drawings followed (their representations' clip objects): a change wakes the views */
    private drawings = new Map<ParticleTileSet, () => void>();
    private readonly drawingsChanged = new Subject<void>();
    private updating = false;
    /** a pass asked for while an update runs: scheduled as it ends */
    private again = false;

    private scheduled = false;
    private disposed = false;
    /**
     * the next pass picks for the camera where it is (the camera changed, what draws a set, a set or view made); else
     * (a tile loaded, a step's next) the views step toward their sets' latest picks: a pick (a pass over the tiles) at
     * most a camera's pass, every `intervalMs` while it moves and once as it stops, not a frame
     */
    private fresh = true;
    /** per tile set its steps' particles a ms (a running mean): what its next step's particles are by */
    private rates = new WeakMap<ParticleTileSet, number>();

    register(): void {
      this.disposed = false;
      // a tile set removed or made again: its loader no longer followed (else it and its set kept for the session)
      const events = this.ctx.state.data.events.object;
      this.subscribeObservable(events.removed, ({ obj }) => this.forget(obj));
      this.subscribeObservable(events.updated, ({ state, ref, action, obj, oldObj, oldData }) => {
        if (action === 'recreate') this.forget(oldObj);
        else if (SO.Particle.TileSet.is(obj) && oldData && oldData !== obj.data) this.forgetData(oldData);
        if (action === 'recreate' || (SO.Particle.TileSet.is(obj) && oldData !== obj.data)) this.wake(state, ref, obj);
      });
      this.subscribeObservable(events.created, ({ state, ref, obj }) => this.wake(state, ref, obj));
      this.subscribeObservable(this.ctx.behaviors.canvas3d.initialized, () => {
        const camera = this.ctx.canvas3d?.camera;
        if (!camera) return;
        // when it last moved (its picks the quicker ones while it moves), a pass at most every `intervalMs`, and one
        // in full once it settled (`ParticleTileDrawings.SettleMs` after its last move); its depth range's changes not
        // a move (`createParticleTilesCameraMoved`)
        // one subscription reading the camera (the check keeps the last matrices: a subscriber each would see a move once)
        const isMoved = createParticleTilesCameraMoved();
        const moved = new Subject<void>();
        this.subscribeObservable(camera.changed, () => {
          if (!isMoved(camera)) return;
          setParticleTilesCameraChanged(camera);
          moved.next();
        });
        const changed = moved.pipe(
          throttleTime(this.params.intervalMs, asyncScheduler, { leading: true, trailing: true }),
        );
        this.subscribeObservable(changed, () => this.schedule(true));
        this.subscribeObservable(moved.pipe(debounceTime(ParticleTileDrawings.SettleMs + 10)), () =>
          this.schedule(true),
        );
      });
      // what draws the sets changing (a cut dragged): a pass at most every `intervalMs`, `DrawingsIntervalMs` at least
      // (a pick a pass on the main thread: a drag's frames kept), its picks the quicker ones, and once it settled
      // (`ParticleTileDrawings.SettleMs`: picked again in full)
      const interval = Math.max(this.params.intervalMs, DrawingsIntervalMs);
      const throttle = throttleTime(interval, asyncScheduler, { leading: true, trailing: true });
      this.subscribeObservable(this.drawingsChanged.pipe(throttle), () => this.schedule(true));
      this.subscribeObservable(this.drawingsChanged.pipe(debounceTime(ParticleTileDrawings.SettleMs + 10)), () =>
        this.schedule(true),
      );
    }

    /**
     * a tile set or a view made (again, e.g. a set of other params with a new loader, a view moved to another set): its
     * loader followed and a pass scheduled, the camera perhaps still (else the view stays at its first step until it
     * moves); a view's update by a step does not (the pass after it does)
     */
    private wake(state: State, ref: string, obj: StateObject | undefined) {
      if (obj && SO.Particle.TileSet.is(obj)) {
        this.watch(obj.data);
        this.schedule(true);
      } else if (state.cells.get(ref)?.transform.transformer === ParticleTilesView) {
        this.schedule(true);
      }
    }

    /** a pass over the following views at the next frame, once however often asked; `fresh`: a pick for the camera */
    private schedule(fresh = false) {
      if (fresh) this.fresh = true;
      if (this.scheduled || this.disposed) return;
      this.scheduled = true;
      nextFrame(() => {
        this.scheduled = false;
        // a step that failed (a view's update throwing) logged, not an unhandled rejection; no step again at once (a
        // failure kept a pass a frame), the next tile loaded or camera change scheduling it
        if (!this.disposed) this.follow().catch((e) => this.ctx.log.error(`Particle tiles: ${e}`));
      });
    }

    /** the following views behind the camera updated a step; again at the next frame while any is behind */
    private async follow() {
      if (this.updating) {
        // after the update in progress (its end schedules the next pass): not a pass a frame meanwhile
        this.again = true;
        return;
      }
      const fresh = this.fresh;
      this.fresh = false;
      const state = this.ctx.state.data;
      const cells = state.select(StateSelection.Generators.ofTransformer(ParticleTilesView));
      const update = state.build();
      const updated: StateObjectCell[] = [];
      const steps: { set: ParticleTileSet; cell: StateObjectCell; before: number; stepMs: number }[] = [];
      for (const cell of cells) {
        const params = cell.transform.params as TilesViewParams | undefined;
        if (!params || !isFollowing(params) || !cell.obj) continue;
        const parent = state.cells.get(cell.transform.parent)?.obj;
        if (!parent || !SO.Particle.TileSet.is(parent)) continue;
        this.watch(parent.data);
        if (!isParticleTilesViewBehind(parent.data, params, cell.obj.data, this.ctx, cell.cache, !fresh)) continue;
        // the same params, a new object: the view applied again, a step toward what the camera picks now
        update.to(cell).update({ ...params, tiles: { ...params.tiles } });
        updated.push(cell);
        const stepMs = params.tiles.name === 'camera' ? (params.tiles.params.stepMs ?? 0) : 0;
        steps.push({ set: parent.data.set, cell, before: particlesOf(cell.obj.data as ParticleList), stepMs });
      }
      if (updated.length === 0) return;
      // a pass that picked for the camera (`fresh`): its step at the next frame, toward that pick (the set's latest),
      // not in this one. Why: a pick and a step in one frame (a few ms each) miss a display's frame where either alone
      // does not
      if (fresh) {
        this.schedule();
        return;
      }
      this.updating = true;
      const start = performance.now();
      try {
        await this.ctx.runTask(state.updateTree(update, { doNotUpdateCurrent: true, doNotLogTiming: true }));
      } finally {
        this.updating = false;
        // a pass asked for meanwhile (a camera change, a tile loaded): also where the update failed
        if (this.again) {
          this.again = false;
          this.schedule();
        }
      }
      // the step's time (the views' and their representations' updates) by the particles each set's views added: a
      // set's share of the step's particles a ms, its next step's particles toward its `stepMs` (steps by time; the
      // sets of a step together about its time)
      const ms = performance.now() - start;
      const addedBySet = new Map<ParticleTileSet, number>();
      for (const step of steps) {
        const added = Math.max(0, particlesOf(step.cell.obj?.data as ParticleList) - step.before);
        addedBySet.set(step.set, (addedBySet.get(step.set) ?? 0) + added);
      }
      if (ms > 0) {
        for (const step of steps) {
          const added = addedBySet.get(step.set)!;
          if (added === 0) continue;
          const measured = added / ms;
          const before = this.rates.get(step.set);
          const rate = before === undefined ? measured : (before + measured) / 2;
          this.rates.set(step.set, rate);
          setParticleTilesStep(step.set, getParticleTilesStepParticles(rate, step.stepMs));
        }
      }
      // the next step after a frame drawn (the tiles ready not added yet, or loaded since)
      this.schedule();
    }

    /** a tile set's object `obj` gone: its loader unsubscribed */
    private forget(obj: StateObject | undefined) {
      if (obj && SO.Particle.TileSet.is(obj)) this.forgetData(obj.data);
    }

    private forgetData(data: SO.Particle.TileSet['data'] | undefined) {
      if (!data) return;
      this.loaders.get(data.loader)?.();
      this.loaders.delete(data.loader);
      this.drawings.get(data.set)?.();
      this.drawings.delete(data.set);
    }

    /**
     * a tile set's changes wake the views: its loader's (a tile loaded), what draws it (a cut moved: a representation's
     * clip objects, a gizmo's drag), the latter at most every `intervalMs` and once more as it settles
     */
    private watch(data: SO.Particle.TileSet['data']) {
      const { loader, set } = data;
      if (!this.loaders.has(loader))
        this.loaders.set(
          loader,
          loader.subscribe(() => this.schedule()),
        );
      if (!this.drawings.has(set)) {
        this.drawings.set(
          set,
          ParticleTileDrawings.subscribe(set, () => this.drawingsChanged.next()),
        );
      }
    }

    update(p: { intervalMs: number }) {
      const changed = p.intervalMs !== this.params.intervalMs;
      this.params = p;
      // the camera's subscription by the new interval: its subscriptions made again
      if (changed) {
        this.dispose();
        this.register();
      }
      return changed;
    }

    dispose() {
      super.dispose();
      this.disposed = true;
      for (const unsubscribe of this.loaders.values()) unsubscribe();
      this.loaders.clear();
      for (const unsubscribe of this.drawings.values()) unsubscribe();
      this.drawings.clear();
    }
  },
});

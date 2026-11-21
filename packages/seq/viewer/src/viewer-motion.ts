/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { resolveColor, type Viewport as RenderViewport } from '@molstar/seqx-core/render';
import type {
  Frame as StoryFrame,
  PlanOptions,
  StepViewport,
  StoryState,
  TransitionPlan,
  ViewState,
} from '@molstar/seqx-core/story';
import { axisTotal, nearestOffsetOf, type SpecModel, type ViewModel } from '@molstar/seqx-viewer/models';
import type { AnimationRequest, SeqAnimation, SeqAnimationParticipant } from '@molstar/seqx-viewer/ports';
import type { SeqStoryPart } from './story.js';
import type { SeqViewer } from './viewer.js';
import type { FrameClock } from './viewer-clock.js';

// The viewer's motion (_spec/README.md "Composition", Motion and Stories): the frame slot (a plan's frame shown over the
// state, its viewports by section), a jump's motion from the viewports drawn to the model's, a story step entered as a
// transition, each run on the application's clock (the animation port) or cut; what the frames drew, where the
// next motion starts. The viewports a frame reads as the slot shows them are `ViewerViewports` (`viewer-viewports.ts`).

export class ViewerMotion {
  /** The plan and frame shown over the state, its viewports by section. */
  shown: { plan: TransitionPlan; frame: StoryFrame; views?: ViewState } | undefined;
  /** A request on the animation port running a plan (a jump's, with its count). */
  private _request: { request: AnimationRequest; plan: TransitionPlan; jumps?: number } | undefined;
  /** Each leading section's viewport the last frame drew, by address: where a motion starts. */
  readonly drawnViews = new Map<string, StepViewport>();
  /** The window the last frame drew (NaN before one). */
  readonly drawnWindow = { start: Number.NaN, end: Number.NaN };
  /** The spec the last frame drew, the step shown, and whether `playStep` is entering one (no second transition). */
  private _drawnRoot: Node | undefined;
  private _shownStep: string | undefined;
  private _stepping = false;
  private readonly _drawnFrames = new State(0);
  /** The application's clock (the animation port): Seq*'s own motions are plans requested of it. */
  private readonly _animation: SeqAnimation;
  /** The participant face of the animation port: what the application's clock drives. */
  readonly participant: SeqAnimationParticipant;
  private readonly spec: SpecModel;
  private readonly view: ViewModel;
  constructor(
    private readonly viewer: SeqViewer,
    private readonly clock: FrameClock,
    /** The stories and transitions the spec lists (`story`): none, every jump a cut and no step entered. */
    private readonly part: SeqStoryPart | undefined,
    animation: (participant: SeqAnimationParticipant) => SeqAnimation,
  ) {
    this.spec = viewer.spec;
    this.view = viewer.view;
    // the animation port (_spec/README.md "Ports", Animation): the application's clock, the page's frames by default
    this.participant = {
      progress: (plan, t) => {
        if (this.part) this._showFrame(plan, this.part.sample(plan, t));
      },
      show: (plan, frame) => this._showFrame(plan, frame),
      clear: () => this._clearFrame(),
      draw: () => {
        this.viewer.draw();
      },
      drawn: this._drawnFrames.asReadonly,
      durationOf: (plan) => plan.duration,
    };
    this._animation = animation(this.participant);
  }

  /** Cancels the running request and disposes the port's clock unless it outlives the viewer. */
  dispose(): void {
    this._request?.request.cancel();
    this._request = undefined;
    if (this._animation.disposeWithViewer !== false) this._animation.dispose?.();
  }

  /** Whether a motion (a jump's, a story step's) is drawing. */
  get moving(): boolean {
    return this._request !== undefined;
  }

  /**
   * Whether a motion, a jump's or a story's, is a cut: `view.motion` (`cut`; `platform`: the platform's reduced motion;
   * `animate`: never), one item for both (_spec/README.md, Configuration).
   */
  private get _cut(): boolean {
    const motion = this.viewer.config.view.motion;
    return motion === 'cut' || (motion === 'platform' && this.viewer.reducedMotion);
  }

  /** Each leading section's viewport in the model now (what a jump goes to), by address. */
  private _modelViews(): ViewState {
    const sections: Record<string, StepViewport> = {};
    for (const s of this.view.resolved.value?.sections ?? []) {
      if (this.view.isFollowing(s.address)) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) continue;
      // the bands' top and zoom (display rows, px a row): what a step tweens; the folds are the session's, never a step's
      const own = this.view.bandsOf(s.address);
      let bands: Record<string, { bandTop: number; rowScale?: number }> | undefined;
      for (const track in own ?? {}) {
        const b = own![track];
        (bands ??= {})[track] = { bandTop: b.bandTop, ...(b.rowScale !== undefined ? { rowScale: b.rowScale } : {}) };
      }
      sections[s.address] = {
        start: v.start,
        end: v.end,
        ...(v.top ? { top: v.top } : {}),
        ...(v.scale !== undefined ? { scale: v.scale } : {}),
        ...(bands ? { bands } : {}),
      };
    }
    return { sections };
  }

  /** What a plan of this viewer needs besides its states: the cyclic sections' totals (a window's shortest arc). */
  private _planOptions(): Pick<PlanOptions, 'cyclic' | 'config' | 'resolveColor'> {
    const cyclic: Record<string, number> = {};
    for (const s of this.view.resolved.value?.sections ?? [])
      if (this.view.wraps(s.address)) {
        const total = axisTotal(s.axis);
        if (total) cyclic[s.address] = total;
      }
    // a token blends through the color it draws in (the viewer's token table, light or dark as shown)
    const tokens = this.viewer.tokens;
    const color = (c: string) => {
      const r = resolveColor(c, tokens);
      // RGBA: an opaque color as hex; one with alpha is no color a plan blends
      return r === undefined || (r & 0xff) !== 0xff ? undefined : `#${(r >>> 8).toString(16).padStart(6, '0')}`;
    };
    return { cyclic, config: this.viewer.config.story, resolveColor: color };
  }

  /**
   * A jump moved the view (_spec/README.md "Composition", Motion): a plan from every leading section's viewport drawn
   * (window, letter size) to the model's, requested of the animation port over `view.jumpMs`, unless motion is cut.
   */
  jumped(jumps: number): void {
    const { jumpMs } = this.viewer.config.view;
    this.end();
    const motion = this.part;
    if (!motion || this.viewer.retired || this._cut || !(jumpMs > 0) || !this.drawnViews.size) return;
    const from: ViewState = { sections: Object.fromEntries(this.drawnViews) };
    const shown = this.spec.shown.value;
    const p = motion.plan(
      { root: shown, view: from },
      { root: shown, view: this._modelViews() },
      { ...this._planOptions(), duration: jumpMs },
    );
    if (p.channels.length) this._run(p, jumps);
  }

  /** Runs a plan on the application's clock, shown from its start now; the frame left at its end. */
  private _run(plan: TransitionPlan, jumps?: number): void {
    this._request?.request.cancel();
    if (!this.part) return;
    this._showFrame(plan, this.part.sample(plan, 0));
    const request = this._animation.request(plan, {
      duration: plan.duration,
      done: () => {
        if (this._request?.plan !== plan) return;
        this._request = undefined;
        this._clearFrame();
      },
    });
    this._request = { request, plan, ...(jumps !== undefined ? { jumps } : {}) };
  }

  /** Ends a running motion where the model is (no frame left over it). */
  end(): void {
    if (!this._request) return;
    this._request.request.cancel();
    this._request = undefined;
    this._clearFrame();
  }

  /** Shows a frame over the state (the frame slot): its spec when the plan changes the spec, its viewports. */
  private _showFrame(plan: TransitionPlan, frame: StoryFrame): void {
    const motion = this.part;
    if (this.viewer.retired || !motion) return;
    // a plan between two equal roots (a jump's) shows the state itself
    this.spec.setFrame(plan.from.root === plan.to.root ? undefined : motion.frameRoot(plan, frame));
    const views = motion.frameView(plan, frame);
    this.shown = { plan, frame, ...(views ? { views } : {}) };
    this.clock.moved();
    this.viewer.schedule();
  }

  /** Leaves the frame: the state is shown and drawn again. */
  private _clearFrame(): void {
    if (!this.shown && this.spec.frame === undefined) return;
    this.shown = undefined;
    this.spec.setFrame(undefined);
    this.viewer.schedule();
  }

  /**
   * A section's viewport moved: a move by hand ends a motion, a story's too, and clears a frame a host shows (a jump's
   * own move is followed by its count moving; a step's view is placed before its motion starts).
   */
  handMoved(): void {
    const r = this._request;
    if (r && (r.jumps === undefined || r.jumps === this.view.jumps.value)) this.end();
    else if (!r && this.shown) this._clearFrame();
  }

  /** The spec's head moved: an edit, an undo or a load under a frame that changes the spec ends it (made from the old root). */
  headMoved(): void {
    const f = this.shown;
    if (!f || f.plan.from.root === f.plan.to.root) return;
    this.end();
    this._clearFrame();
  }

  /** Stops a motion or a host's frame where it is drawn (`SeqViewer.stopMotion`). */
  stop(): void {
    const r = this._request;
    if (!r && !this.shown) return;
    const drawn = [...this.drawnViews];
    r?.request.cancel();
    this._request = undefined;
    this._clearFrame();
    for (const [address, v] of drawn) this.view.setViewport(address, v);
    this.viewer.draw();
  }

  /** Enters a story step, or leaves the step shown, as a transition (`SeqViewer.playStep`). */
  playStep(id: string | undefined, options: { readonly run?: boolean }): { plan?: TransitionPlan; issue?: Issue } {
    this.end();
    const from = this._fromState(this.spec.shown.value);
    this._stepping = true;
    try {
      if (id === undefined) this.spec.leaveStep();
      else {
        const issue = this.spec.enterStep(id);
        if (issue) return { issue };
      }
    } finally {
      this._stepping = false;
    }
    this._shownStep = this.spec.step.value?.id;
    return { plan: this._toStep(from, options.run ?? true) };
  }

  /** What a transition starts from: a root and the viewports drawn last (the model's before a frame). */
  private _fromState(root: Node): StoryState {
    const views = this.drawnViews.size ? { sections: Object.fromEntries(this.drawnViews) } : this._modelViews();
    return { root, view: { ...views, ...this._viewAddress() } };
  }

  /** A step entered elsewhere (the `step.enter` command, a host): the same transition, from what was drawn. */
  stepChanged(): void {
    const id = this.spec.step.value?.id;
    if (this._stepping || this.viewer.retired || id === this._shownStep) return;
    this._shownStep = id;
    this.end();
    this._toStep(this._fromState(this._drawnRoot ?? this.spec.shown.value), true);
  }

  /** The step shown (its view placed) and the plan to it from `from`, run or shown at its start (none: no motion). */
  private _toStep(from: StoryState, run: boolean): TransitionPlan | undefined {
    const st = this.spec.step.value?.state;
    if (st?.view) {
      if (st.view.view !== undefined) this.view.show(st.view.view);
      const resolved = this.view.resolved.value;
      for (const [address, v] of Object.entries(st.view.sections)) {
        const axis = resolved?.sections.find((s) => s.address === address)?.axis;
        const s = axis && nearestOffsetOf(axis, v.start);
        const e = axis && nearestOffsetOf(axis, v.end);
        if (s === undefined || e === undefined) continue;
        this.view.setViewport(address, {
          start: s,
          end: e,
          ...(v.top !== undefined ? { top: v.top } : {}),
          ...(v.scale !== undefined ? { scale: v.scale } : {}),
        });
        // a step's bands: their top and zoom placed in the view model, the session's folds kept
        for (const track in v.bands ?? {}) {
          const b = v.bands![track];
          this.view.setBand(track, {
            ...(b.bandTop !== undefined ? { bandTop: b.bandTop } : {}),
            ...(b.rowScale !== undefined ? { rowScale: b.rowScale } : {}),
          });
        }
      }
    }
    const motion = this.part;
    if (!motion) {
      this.viewer.schedule();
      return undefined;
    }
    const to = { root: this.spec.shown.value, view: { ...this._modelViews(), ...this._viewAddress() } };
    const timing = st?.timing;
    const p = motion.plan(from, to, {
      ...this._planOptions(),
      reducedMotion: this._cut,
      ...(timing ? { duration: timing.duration, easing: timing.easing } : {}),
    });
    if (this._cut || !p.channels.length) {
      this.viewer.schedule();
      return p;
    }
    if (run) this._run(p);
    else this._showFrame(p, motion.sample(p, 0));
    return p;
  }

  /** The view shown's address, when one is shown by address. */
  private _viewAddress(): { view?: string } {
    const v = this.view.view.value;
    return v !== undefined ? { view: v } : {};
  }

  /**
   * A leading section's viewport as drawn now: the frame's (an animation's), else the model's. `eye` false: a jump's
   * motion is the eye's way only (the SVG and `viewport` read the model, the destination), a story's frame is shown.
   */
  shownViewport(address: string, eye = true): StepViewport | undefined {
    const jump =
      this._request !== undefined && this._request.jumps !== undefined && this._request.plan === this.shown?.plan;
    const f = eye || !jump ? this.shown?.views?.sections[address] : undefined;
    if (f) return f;
    const v = this.view.get(address);
    return v && !v.pending ? v : undefined;
  }

  /**
   * What a frame drew: its window and each leading section's viewport, where the next motion starts (kept in place:
   * `Map.clear` makes a new table a frame); a motion is not settled while it runs.
   */
  laid(window: RenderViewport | undefined, viewports: Readonly<Record<string, RenderViewport>> | undefined): void {
    if (window) {
      this.drawnWindow.start = window.start;
      this.drawnWindow.end = window.end;
    }
    if (viewports) {
      for (const k of this.drawnViews.keys()) if (!(k in viewports)) this.drawnViews.delete(k);
      for (const k in viewports) this.drawnViews.set(k, viewports[k]);
    } else if (this.drawnViews.size) this.drawnViews.clear();
    if (this._request) this.clock.moved();
  }

  /** A frame was drawn: the state's root it drew, or one more frame of the slot (what the port's `drawn` counts). */
  drawn(): void {
    if (!this.shown) this._drawnRoot = this.spec.shown.value;
    if (this.shown) this._drawnFrames.set(this._drawnFrames.value + 1);
  }
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { clampSectionTops } from './section-scroll.js';
import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { Issue } from '@molstar/core/repr/tree';
import { type Frame, type FrameInput, type Placement, Renderer } from '@molstar/seqx-core/render';
import type { SeqTokens } from '@molstar/seqx-core/render';
import type { HoverCardBehavior, RowsBehavior } from '@molstar/seqx-viewer/behaviors';
import type { CanvasBackend } from '@molstar/seqx-viewer/canvas';
import type { SeqViewer } from './viewer.js';
import type { FrameClock } from './viewer-clock.js';
import type { ViewerEmphasis } from './viewer-emphasis.js';
import type { ViewerMotion } from './viewer-motion.js';
import type { ViewerRails } from './viewer-rails.js';
import type { ViewerViewports } from './viewer-viewports.js';

// The viewer's frame (_spec/README.md "Composition"): the renderer, made again when the model or the view changes; a
// frame lays the view out when anything it reads changed (else only the emphasis passes are redrawn, a hover), asks
// the store for what the layout wants, follows the content's height, re-hits a resting pointer, and the Canvas2D
// backend draws what changed. One reused input record: a warm frame allocates nothing.

type Mutable<T> = { -readonly [K in keyof T]: T[K] };

export class ViewerFrame {
  private _backend: CanvasBackend | undefined;
  private _renderer: Renderer | undefined;
  private _frame: Frame | undefined;
  /** The backend draws the base passes again on the next frame (a new backend, a token table). */
  private _baseDirty = true;
  private readonly _rendererFor: unknown[] = [];
  private _wanter: { want(w: readonly any[]): void; dispose(): void } | undefined;
  private _contentHeight = 0;
  private readonly _size = { width: 0, height: 0 };
  private readonly _input = { size: this._size } as Mutable<FrameInput>;
  /** Bumped after every frame drawn: a navigator redraws over it. */
  private readonly _drawn = new State(0);
  readonly drawn: ReadonlyState<number> = this._drawn.asReadonly;

  constructor(
    private readonly viewer: SeqViewer,
    private readonly fixedDpr: number | undefined,
    private readonly clock: FrameClock,
    private readonly motion: ViewerMotion,
    /** The window and the viewports a frame reads, as the frame slot shows them. */
    private readonly viewports: ViewerViewports,
    private readonly rails: ViewerRails,
    private readonly emphasis: ViewerEmphasis,
    /** The layout's issues of the last frame, set when it lays out (`SeqViewer.layoutIssues`). */
    private readonly issues: State<readonly Issue[]>,
  ) {}

  /** The renderer of the last frame (made again when the model or the view changes). */
  get renderer(): Renderer | undefined {
    return this._renderer;
  }

  /** The last frame drawn. */
  get frame(): Frame | undefined {
    return this._frame;
  }

  /** Attaches the backend that draws the lists (none: the frame lays out and hit-tests only); it draws everything next. */
  useBackend(backend: CanvasBackend | undefined): void {
    this._backend = backend;
    this._baseDirty = true;
  }

  /** The token table changed: the renderer keeps its layout key's token version, the backend redraws. */
  setTokens(tokens: SeqTokens): void {
    if (this._renderer) this._renderer.tokens = tokens;
    this._baseDirty = true;
  }

  /** A section's placement in the last frame (its strips, what its layout's arithmetic reads), by address. */
  placementOf(section: string): Placement | undefined {
    const sections = this._frame?.layout.sections;
    if (!sections) return undefined;
    for (let i = 0; i < sections.length; i++) if (sections[i].section.address === section) return sections[i].placement;
    return undefined;
  }

  /** Lets go of what the store was asked for (the viewer retired). */
  dispose(): void {
    this._wanter?.dispose();
  }

  /** Draws a frame now (`SeqViewer.draw`). */
  draw(): Frame | undefined {
    const viewer = this.viewer;
    if (viewer.retired) return undefined;
    // the model the view is resolved over (a source's complexes in it: core's _spec "Overlay", Opened sources)
    const model = viewer.view.model.value;
    const resolved = viewer.view.resolved.value;
    if (!resolved || viewer.width <= 0) return undefined;
    const config = viewer.config;
    const store = viewer.store;
    const context = viewer.context;
    const rr = this._rendererFor;
    if (!this._renderer || rr[0] !== model || rr[1] !== resolved) {
      this._renderer = new Renderer({
        model,
        view: resolved,
        data: store.store,
        tokens: viewer.tokens,
        text: context.text,
        config: context.renderConfig,
        // the context's as they stand: an extension's coming or going resolves the view again, so a new renderer
        types: context.types,
        ...(context.schemes ? { theme: { schemes: context.schemes } } : {}),
        layouts: viewer.layouts,
        sizings: context.sizings,
        ...(context.trackSizing ? { trackSizing: context.trackSizing } : {}),
      });
      rr[0] = model;
      rr[1] = resolved;
      this._wanter?.dispose();
      this._wanter = store.store.wanter();
    }
    const renderer = this._renderer;
    const dpr =
      this.fixedDpr ?? (config.canvas.dpr === 'device' ? (globalThis.devicePixelRatio ?? 1) : config.canvas.dpr);
    viewer.dpr = dpr;

    const { height: h } = config.canvas;
    const fixed = h !== 'content';
    const size = this._size;
    size.width = viewer.width;
    size.height = fixed ? h : Math.min(config.canvas.maxHeight, Math.max(this._contentHeight, 32));
    this.emphasis.update();
    const input = this._input;
    // every leading section its own viewport, read by its layout (_spec/README.md "Composition"); the first's once
    const vp = this.viewports.window(resolved);
    const viewports = this.viewports.all(resolved);
    if (viewports) input.viewports = viewports;
    else delete input.viewports;
    const rails = this.rails;
    rails.update(resolved);
    // the rows' zoom the host configured (the rows rails' bounds and grips)
    input.rowsZoom = config.rows.zoom;
    if (rails.rails) {
      input.rails = rails.rails;
      input.railsActive = rails.active;
    } else {
      delete input.rails;
      delete input.railsActive;
    }
    // what each leading section showed: where the next motion starts (a motion is not settled while it runs)
    this.motion.laid(vp, viewports);
    input.dpr = dpr;
    input.settled = this.clock.quiet && store.store.settled.value;
    input.emphasis = this.emphasis.emphasis;
    input.cvdSafe = viewer.legends.cvdSafe.value;
    // the header entries' parts (the `headers` behavior's, by name: a viewer without it bundles none of it)
    const headers = viewer.behaviors.get('headers')?.partsOf;
    if (headers) input.headers = headers;
    else if (input.headers) delete input.headers;
    // tier 2: the store's dynamic nodes (its version is part of the frame's key)
    input.overlay = store.overlay(viewer.spec.model.value, resolved, context);

    let frame = renderer.frame(input);
    let rebuilt = frame.rebuilt;
    if (!fixed) {
      const last = frame.layout.sections.at(-1);
      const content = last
        ? Math.ceil(last.y + last.height + context.renderConfig.inset + (frame.layout.pinnedPx ?? 0))
        : 32;
      if (content !== this._contentHeight) {
        this._contentHeight = content;
        size.height = Math.min(config.canvas.maxHeight, Math.max(content, 32));
        frame = renderer.frame(input);
        rebuilt ||= frame.rebuilt;
      }
    }
    viewer.height = size.height;
    if (rebuilt) {
      viewer.counts.layouts++;
      // what the models read of the placement (a text's lines in view) follows the frame
      viewer.view.placed();
      // a scrolled band's stored top kept inside its lanes (by name: a viewer without the rows behavior bundles none)
      (viewer.behaviors.get('rows') as RowsBehavior | undefined)?.framed(frame.layout);
      // a section's stored scroll kept within its content (its tracks hidden or shrunk while scrolled)
      clampSectionTops(viewer.view, frame.layout);
      if (!this.clock.restDue) this.clock.changed();
      this._wanter?.want(frame.wants);
      this.issues.set(frame.layout.issues);
      // legends on settle only (a pan's frames derive none), compared structurally
      if (input.settled) viewer.legends.update(frame.layout);
      store.holdShown(frame.layout.shown);
    }
    this._frame = frame;
    // the layout placed the view (its width): a rail that went idle or live with it is drawn so now (the window pass
    // alone: which rails can move is outside the key)
    const wasActive = rails.active;
    rails.update(resolved);
    if (rails.rails && rails.active !== wasActive) {
      input.railsActive = rails.active;
      frame = renderer.frame(input);
      this._frame = frame;
    }

    // a resting pointer is hit again once a frame the view moves under it (EC-143)
    // (what that frame laid out is still to be drawn: `rebuilt` is the layout's, not the emphasis redraw's)
    if (rebuilt && config.pointer.rehit) {
      this.clock.rehitting = true;
      const changed = viewer.input.rehit();
      this.clock.rehitting = false;
      if (changed) {
        this.emphasis.update();
        frame = renderer.frame(input);
        this._frame = frame;
      }
    }

    const backend = this._backend;
    if (backend) {
      backend.fonts.setStacks(renderer.tokens.fonts);
      this._baseDirty = backend.resize(size.width, size.height, dpr) || this._baseDirty;
      backend.draw(frame.lists, rebuilt || frame.base || this._baseDirty);
      this._baseDirty = false;
    }
    viewer.counts.frames++;
    this._drawn.set(this._drawn.value + 1);
    // the hover card follows the frame: a cursor's card after the frame that revealed it, none outliving a new layout
    // (by name: a viewer without the card bundles none of it)
    (viewer.behaviors.get('hoverCard') as HoverCardBehavior | undefined)?.framed(rebuilt);
    // the headers' geometry and parts follow the frame (by name, as the card)
    viewer.behaviors.get('headers')?.framed();
    this.motion.drawn();
    // the settled picture is drawn: it rests now, once
    this.clock.drawn();
    return frame;
  }
}

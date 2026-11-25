/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Emphasis, EmphasisHighlight, Hit as RenderHit, RenderSvgOptions } from '@molstar/seqx-core/render';
import type { Location } from '@molstar/seqx-core/query';
import { isNothing, Projector } from '@molstar/seqx-viewer/models';
import { highlightColorOf } from './legends.js';
import type { SeqViewer } from './viewer.js';

// The emphasis the passes draw (_spec/README.md "Emphasis"): the selection's Location and its projection onto the
// linked sections, the cursor, the highlight channels (each with its projection) and the hover with its projection,
// in reused records set before each frame; the hover's hit as the emphasis and the hit-test read it; the session as an
// export draws it.

type Mutable<T> = { -readonly [K in keyof T]: T[K] };

export class ViewerEmphasis {
  readonly emphasis: Mutable<Emphasis> = {};
  /** The projections onto the view's sections, cached (`Projector`): a redraw makes none. */
  readonly projector: Projector;
  /** The emphasis's projected record, reused: its fields are the cached projections. */
  private readonly _projected: { hover?: Location; selection?: Location; hoverSection?: string } = {};
  private readonly _cursor: {
    system: string;
    position: number;
    section: string | undefined;
    ySystem: string | undefined;
    yPosition: number | undefined;
    part: 'cell' | 'header' | undefined;
    track: string | undefined;
  } = {
    system: '',
    position: 0,
    section: undefined,
    ySystem: undefined,
    yPosition: undefined,
    part: undefined,
    track: undefined,
  };
  /** The highlights as drawn, made again only when a channel, the tokens, the view or the model changed. */
  private _highlights: readonly EmphasisHighlight[] = [];
  private readonly _highlightsFor: unknown[] = [];
  private readonly _hit = {} as Mutable<RenderHit>;
  /** Section ids by address (the last step), kept: a hover makes no string. */
  private readonly _sectionIds = new Map<string, string>();
  private _hoverVersion = -1;

  constructor(private readonly viewer: SeqViewer) {
    this.projector = new Projector(viewer.spec, viewer.view, viewer.store.regions);
  }

  /**
   * Sets the emphasis the passes draw. Reused records, mutated when what they show changes; the projections are cached
   * (`Projector`, `HoverModel.projectedAll`), so a redraw makes none.
   */
  update(): void {
    const viewer = this.viewer;
    const e = this.emphasis;
    const p = this._projected;
    const sel = viewer.selection.value.value;
    e.elements = sel.elements.keys.length ? sel.elements : undefined;
    if (isNothing(sel)) e.selection = p.selection = undefined;
    else {
      e.selection = sel.shown;
      p.selection = this.projector.of(sel.shown);
    }
    // the cursor shows once the keyboard moved it (focus-visible: a click places it unseen)
    const c = viewer.selection.cursor.value;
    if (c?.keyboard) {
      const k = this._cursor;
      k.system = c.system;
      k.position = c.position;
      // on a section of two axes its row and its part, drawn by that section's layer (the others mark its column)
      k.section = c.section;
      k.ySystem = c.ySystem;
      k.yPosition = c.yPosition;
      k.part = c.part;
      // on a track's header: its actions drawn
      k.track = c.track;
      e.cursor = k;
    } else e.cursor = undefined;
    // the header whose menu is open keeps its actions shown (the pointer leaves the line on its way to the menu)
    const menu = (
      viewer.behaviors.get('headers') as
        | { menu?: { value: { menu: { target: { section: string; track?: string } } } | undefined } }
        | undefined
    )?.menu?.value?.menu.target;
    e.menu = menu ? (menu.track ?? menu.section) : undefined;
    this.updateHighlights();
    e.projected = p;
    p.hover = viewer.hover.projectedAll();
    p.hoverSection = viewer.hover.current.section;
    // a window's part under the pointer (a rail's, a ring's) lights it; a window being dragged grows its rail
    e.drag = viewer.input?.dragHit;
    const wh = viewer.input?.windowHit;
    if (wh) {
      e.hover = wh;
      this._hoverVersion = -1;
      return;
    }
    const v = viewer.hover.version.value;
    if (v === this._hoverVersion) return;
    this._hoverVersion = v;
    e.hover = this.hoverHit();
  }

  /** The hover model's hit as the emphasis (and the hits) read it, in the reused record; undefined with no hover. */
  hoverHit(): RenderHit | undefined {
    const r = this.viewer.hover.current;
    if (r.owner === undefined || r.layer === undefined || r.system === undefined) return undefined;
    const h = this._hit;
    h.layer = r.layer;
    h.section = r.section === undefined ? '' : this._sectionId(r.section);
    h.system = r.system;
    (h as { position: number | undefined }).position = Number.isNaN(r.position) ? undefined : r.position;
    (h as { element: unknown }).element = r.element;
    (h as { extent: unknown }).extent = r.extent;
    // a cell of two axes (a grid's, a rows cell): its second position, what the crosshair draws
    (h as { ySystem: string | undefined }).ySystem = r.ySystem;
    (h as { yPosition: number | undefined }).yPosition = Number.isNaN(r.yPosition) ? undefined : r.yPosition;
    // the clade hover of a rows layer reads it (_spec/README.md "Render", The rows layer, Hover)
    (h as { target: unknown }).target = r.target;
    // the row column cell a header hover is on: its cut value read whole in the tip pass (Emphasis along y)
    (h as { state: string | undefined }).state = r.state;
    return h;
  }

  /**
   * The highlight channels as drawn: each its Location on every section it shows on and its element keys, in its
   * channel's color (the export's too, `svgEmphasis`).
   */
  updateHighlights(): void {
    const viewer = this.viewer;
    const k = this._highlightsFor;
    const version = viewer.highlight.version.value;
    const view = viewer.view.resolved.value;
    const model = viewer.view.model.value;
    if (k[0] === version && k[1] === viewer.tokens.version && k[2] === view && k[3] === model) {
      this.emphasis.highlights = this._highlights.length ? this._highlights : undefined;
      return;
    }
    k[0] = version;
    k[1] = viewer.tokens.version;
    k[2] = view;
    k[3] = model;
    const out: EmphasisHighlight[] = [];
    for (const name of viewer.highlight.names.value) {
      const h = viewer.highlight.get(name);
      if (!h || isNothing(h)) continue;
      const location = this.projector.of(h.shown, { own: true }) ?? h.shown;
      // its element keys too: a layer of keyed items draws those it holds (a pair's arc, EmphasisHighlight.elements)
      out.push({
        location,
        color: highlightColorOf(viewer, name),
        ...(h.elements.keys.length ? { elements: h.elements } : {}),
        ...(h.style ? { style: h.style } : {}),
      });
    }
    this._highlights = out;
    this.emphasis.highlights = out.length ? out : undefined;
  }

  /** The session as an export draws it: the selection (its address, else its Location) and the highlight channels. */
  svgEmphasis(): Pick<RenderSvgOptions, 'selection' | 'emphasis'> {
    this.updateHighlights();
    const sel = this.viewer.selection.value.value;
    const highlights = this.emphasis.highlights;
    const location = sel.address === undefined && !isNothing(sel) ? sel.shown : undefined;
    return {
      ...(sel.address !== undefined ? { selection: sel.address } : {}),
      ...(location || highlights
        ? { emphasis: { ...(location ? { selection: location } : {}), ...(highlights ? { highlights } : {}) } }
        : {}),
    };
  }

  /** A section's id, its address's last step (kept: a hover allocates nothing once warm). */
  private _sectionId(address: string): string {
    let id = this._sectionIds.get(address);
    if (id === undefined) this._sectionIds.set(address, (id = address.slice(address.lastIndexOf('/') + 1)));
    return id;
  }
}

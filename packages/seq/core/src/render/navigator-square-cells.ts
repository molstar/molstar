/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node, Value } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import type { Overlay, StaticModel } from '@molstar/seqx-core/resolve';
import { type LayerKinds, type ResolvedNavigator, type ResolvedView, resolveViewNode } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import type { PassWriter } from './draw-lists.js';
import type { LayerType, SquareCells, SquareLattice, StoreReader } from './layer-type.js';
import type { NavigatorPicture } from './navigator-picture.js';
import { Renderer, type RendererOptions } from './renderer.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';

// A square's picture (_spec/README.md "Render", Navigators, The square, A picture's cells): when the navigator's
// `content` names a registered picture whose track's layer type has `square`, the square lays that track out over the
// assembly's blocks as a hidden section of its own (a renderer of its own over the viewer's store: the layer reads its
// data as on any row, its wants the navigator's) and draws the type's cells over its lattice. Without one, or without a
// store, the lattice alone (the `axis` picture's square).

/** What a square's picture is given beyond the square's own options: the composition's and the store's. */
export interface SquarePictureOptions {
  readonly data?: StoreReader;
  readonly text?: TextMetrics;
  readonly types?: RendererOptions['types'];
  readonly layouts?: RendererOptions['layouts'];
  readonly pictures?: Readonly<Record<string, NavigatorPicture>>;
  readonly kinds?: LayerKinds;
}

/** The section's id inside the square's hidden view. */
const HiddenSection = 'cells';
/** The hidden section's height, CSS px: its layer draws nothing there, so any. */
const HiddenHeight = 32;

/** A square's picture: its hidden renderer and the cells' layer type, laid out when the square asks. */
export class SquarePictureCells {
  private readonly renderer: Renderer;
  private readonly cells: SquareCells<unknown>;
  private lay: unknown = undefined;
  private settled: boolean | undefined = undefined;
  private waiting = true;
  /** The hidden layout's wants (the navigator's, on a frame that laid it out again). */
  wants: readonly Want[] = [];
  private readonly input = {
    size: { width: 1, height: HiddenHeight },
    dpr: 1,
    settled: false,
    overlay: undefined as Overlay | undefined,
  };

  private constructor(
    renderer: Renderer,
    private readonly type: LayerType<any, any>,
  ) {
    this.renderer = renderer;
    this.cells = type.square!;
  }

  /**
   * The picture a square draws, or undefined: no picture named, no store, no track, or its layer's type has no cells.
   * `of` is the lattice's assembly (what the hidden section shows: every block).
   */
  static of(
    model: StaticModel,
    view: ResolvedView,
    nav: ResolvedNavigator,
    of: Value,
    tokens: SeqTokens,
    config: RenderConfig,
    o: SquarePictureOptions,
  ): SquarePictureCells | undefined {
    const entry = nav.content ? o.pictures?.[nav.content] : undefined;
    if (!entry || !o.data || !o.text || !o.types || !o.layouts) return undefined;
    const bound = view.sections.find((s) => s.address === nav.of[0]);
    if (!bound) return undefined;
    let boundNode: Node | undefined = model.spec;
    for (const id of bound.address.split('/').filter(Boolean))
      boundNode = boundNode?.children?.find((c) => c.id === id);
    const picture = entry.picture({
      model,
      view,
      navigator: nav,
      bound,
      boundNode,
      height: 0,
      zoom: false,
      config,
      left: new Set(),
    });
    const tracks = picture?.tracks;
    if (!tracks?.length) return undefined;
    const node: Node = {
      kind: 'view',
      id: nav.id,
      params: { layout: 'tracks' },
      children: [{ kind: 'section', id: HiddenSection, params: { of }, children: tracks as Node[] }],
    };
    const hidden = resolveViewNode(model, node, nav.address, o.kinds);
    const layer = hidden.sections[0]?.tracks[0]?.layers[0];
    const type = layer ? o.types[layer.type] : undefined;
    if (!type?.square) return undefined;
    const renderer = new Renderer({
      model,
      view: hidden,
      data: o.data,
      tokens,
      text: o.text,
      config,
      types: o.types,
      layouts: o.layouts,
      picture: true,
    });
    return new SquarePictureCells(renderer, type);
  }

  setTokens(tokens: SeqTokens): void {
    this.renderer.tokens = tokens;
  }

  /** Whether `layout` would lay the hidden section out now (the base drawn again, the store's settled state moved, waiting). */
  due(settled: boolean, rebuilt: boolean): boolean {
    return rebuilt || this.waiting || settled !== this.settled;
  }

  /**
   * Lays the hidden section out when the square's base is drawn again, the store's settled state moved, or the cells
   * still wait on data (never on a window's move: nothing allocated there). True when the cells changed.
   */
  layout(width: number, dpr: number, settled: boolean, overlay: Overlay | undefined, rebuilt: boolean): boolean {
    if (!this.due(settled, rebuilt)) return false;
    this.settled = settled;
    const fi = this.input;
    fi.size.width = Math.max(1, width);
    fi.dpr = dpr;
    fi.settled = settled;
    fi.overlay = overlay;
    const f = this.renderer.frame(fi as never);
    if (!f.rebuilt) return false;
    this.wants = f.wants;
    const s = f.layout.sections[0];
    let lay: unknown;
    for (const l of s?.sectionLays ?? []) if (l.type === this.type) lay = l.lay;
    for (const band of s?.lays ?? []) for (const l of band) if (!lay && l.type === this.type) lay = l.lay;
    this.lay = lay;
    this.waiting = lay === undefined || !!this.cells.pending?.(lay);
    return true;
  }

  /** Draws the cells (the square's base pass). */
  draw(out: PassWriter, lattice: SquareLattice): void {
    if (this.lay !== undefined) this.cells.draw(this.lay, out, lattice);
  }

  /** What the square says while its cells have nothing to show (the hatch's words), else undefined. */
  note(): string | undefined {
    return this.lay !== undefined ? this.cells.note?.(this.lay) : undefined;
  }

  /** What a cell says. */
  text(i: number, j: number, lattice: SquareLattice): string | undefined {
    return this.lay !== undefined ? this.cells.text?.(this.lay, i, j, lattice) : undefined;
  }
}

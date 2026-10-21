/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Mat4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { Canvas3D } from '@molstar/graphics/canvas3d/canvas3d';
import type { Geometry, GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Representation, RepresentationProvider } from '@molstar/graphics/repr/representation';
import { ShapeRepresentation, type ShapeGetter } from '@molstar/graphics/repr/shape/representation';
import { Theme, type ThemeDataContext, type ThemeRegistryContext } from '@molstar/graphics/theme/theme';

/** Any representation, whatever its params and state (`AnyRepresentation` fixes them to the base ones) */
type AnyRepresentation = Representation<any, any, any>;

/** What a visual set needs of the canvas: its public scene API (no Canvas3D internals) */
export type VisualSetCanvas = Pick<Canvas3D, 'graphics' | 'add' | 'remove' | 'update' | 'requestDraw'>;

export interface VisualSetProviderOptions<P extends PD.Params> {
  /** The theme registries (a plugin's `representation.structure.themes`) */
  readonly themes: ThemeRegistryContext;
  /** The theme's data: `{ structure }`, `{ volume }` */
  readonly themeData: ThemeDataContext;
  /** For the custom properties the representation and its themes need */
  readonly properties: CustomProperty.Context;
  readonly props?: Partial<PD.Values<P>>;
  /** The provider's default themes when not given */
  readonly colorTheme?: { name: string; params?: {} };
  readonly sizeTheme?: { name: string; params?: {} };
  readonly state?: Partial<Representation.State>;
}

interface Entry {
  /** What was attached for it (custom properties, theme dependencies), released with it */
  readonly releases: (() => void)[];
}

/**
 * Representations drawn in one canvas outside the state tree: gizmos, guides, shadow molecules an editor previews. The
 * set adds them to the canvas, tells its own from others' (`has`, `owns` for an interaction event's loci) and removes
 * and destroys them, one or all at once.
 *
 * They stay in the canvas's scene as any representation (picked, marked, in its bounds) and are re-created by the
 * canvas when its graphics backend changes. They are auto-synced: an update (`update`, or the representation's own)
 * reaches the canvas without a commit of the state.
 */
export class VisualSet {
  private readonly entries = new Map<AnyRepresentation, Entry>();
  private disposed = false;

  constructor(private readonly canvas: VisualSetCanvas) {}

  get size() {
    return this.entries.size;
  }

  get isDisposed() {
    return this.disposed;
  }

  get representations(): readonly AnyRepresentation[] {
    return [...this.entries.keys()];
  }

  has(repr: AnyRepresentation | undefined): boolean {
    return !!repr && this.entries.has(repr);
  }

  /** Whether an interaction event's loci (hover, click, drag) is on one of the set's representations */
  owns(loci: Representation.Loci | undefined): boolean {
    return this.has(loci?.repr);
  }

  /** Adds a representation that is already built (and owns it from now on: `remove` and `clear` destroy it) */
  add<R extends AnyRepresentation>(repr: R, state?: Partial<Representation.State>): R {
    return this.track(repr, { releases: [] }, state);
  }

  /**
   * Builds a representation with the canvas's graphics, creates it from `props` and `data` and adds it. When the set is
   * cleared or disposed before that finishes, the representation is destroyed and not added (the promise still
   * resolves to it).
   */
  async create<R extends AnyRepresentation>(
    build: (graphics: GraphicsContext) => R,
    props?: Partial<PD.Values<R['params']>>,
    data?: unknown,
    state?: Partial<Representation.State>,
  ): Promise<R> {
    this.assertLive();
    const repr = build(this.canvas.graphics);
    return this.finish(repr, { releases: [] }, state, () => repr.createOrUpdate(props, data).run());
  }

  /** A shape (gizmos, guides): `ShapeRepresentation` of `getShape` over `data` */
  addShape<D, G extends Geometry, P extends Geometry.Params<G>>(
    getShape: ShapeGetter<D, G, P>,
    geometryUtils: GeometryUtils<G>,
    data: D,
    props?: Partial<PD.Values<P>>,
    state?: Partial<Representation.State>,
  ): Promise<ShapeRepresentation<D, G, P>> {
    return this.create(
      (graphics) => ShapeRepresentation(getShape, geometryUtils, { graphics }),
      props as any,
      data,
      state,
    );
  }

  /**
   * A representation of a provider (a structure's cartoon, a volume's isosurface) over `data`, themed as the state
   * tree's would be: its custom properties and its themes' dependencies attached, released when it is removed
   */
  async fromProvider<D, P extends PD.Params, S extends Representation.State>(
    provider: RepresentationProvider<D, P, S>,
    data: D,
    options: VisualSetProviderOptions<P>,
  ): Promise<Representation<D, P, S>> {
    this.assertLive();
    const { themes, themeData, properties } = options;
    const themeProps = {
      colorTheme: options.colorTheme ?? {
        name: provider.defaultColorTheme.name,
        params: provider.defaultColorTheme.props ?? {},
      },
      sizeTheme: options.sizeTheme ?? {
        name: provider.defaultSizeTheme.name,
        params: provider.defaultSizeTheme.props ?? {},
      },
    };
    const repr = provider.factory({ ...themes, graphics: this.canvas.graphics }, provider.getParams);
    const entry: Entry = { releases: [] };
    return this.finish(repr, entry, options.state, async () => {
      const custom = provider.ensureCustomProperties;
      if (custom) {
        await custom.attach(properties, data);
        entry.releases.push(() => custom.detach(data));
      }
      await Theme.ensureDependencies(properties, themes, themeData, themeProps);
      entry.releases.push(() => Theme.releaseDependencies(themes, themeData, themeProps));
      repr.setTheme(Theme.create(themes, themeData, themeProps));
      const props = { ...PD.getDefaultValues(provider.getParams(themes, data)), ...options.props } as PD.Values<P>;
      await repr.createOrUpdate(props, provider.getData?.(data, props) ?? data).run();
    });
  }

  /** Updates one of the set's representations; the canvas follows (auto-sync) */
  async update<R extends AnyRepresentation>(repr: R, props?: Partial<PD.Values<R['params']>>, data?: unknown) {
    if (!this.has(repr)) return;
    await repr.createOrUpdate(props, data).run();
    this.canvas.requestDraw();
  }

  setState(repr: AnyRepresentation, state: Partial<Representation.State>) {
    if (!this.has(repr)) return;
    repr.setState(state);
    // the scene's bounds follow (a moved gizmo: the camera's near and far planes, culling, picking)
    this.canvas.update(repr);
    this.canvas.requestDraw();
  }

  /** Moves a representation as a whole (a gizmo's drag), without rebuilding it */
  setTransform(repr: AnyRepresentation, transform: Mat4) {
    this.setState(repr, { transform });
  }

  remove(repr: AnyRepresentation) {
    const entry = this.entries.get(repr);
    if (!entry) return;
    this.entries.delete(repr);
    this.release(repr, entry);
    this.canvas.requestDraw();
  }

  /** Removes and destroys every representation; the set can be used again */
  clear() {
    this.generation++;
    if (this.entries.size === 0) return;
    const entries = [...this.entries];
    this.entries.clear();
    for (const [repr, entry] of entries) this.release(repr, entry);
    this.canvas.requestDraw();
  }

  /** `clear`, and no representation can be added anymore */
  dispose() {
    if (this.disposed) return;
    this.clear();
    this.disposed = true;
  }

  /** Bumped by `clear`: a `create` that started before it does not add what it built */
  private generation = 0;

  private async finish<R extends AnyRepresentation>(
    repr: R,
    entry: Entry,
    state: Partial<Representation.State> | undefined,
    run: () => Promise<void>,
  ): Promise<R> {
    const generation = this.generation;
    try {
      await run();
    } catch (e) {
      this.release(repr, entry, false);
      throw e;
    }
    if (this.disposed || generation !== this.generation) {
      this.release(repr, entry, false);
      return repr;
    }
    return this.track(repr, entry, state);
  }

  private track<R extends AnyRepresentation>(repr: R, entry: Entry, state?: Partial<Representation.State>): R {
    this.assertLive();
    if (this.entries.has(repr)) return repr;
    // auto-synced by the canvas (`Canvas3D.add` follows `repr.updated`): no state commit adds it
    repr.setState({ ...state, syncManually: false });
    this.entries.set(repr, entry);
    this.canvas.add(repr);
    this.canvas.requestDraw();
    return repr;
  }

  private release(repr: AnyRepresentation, entry: Entry, added = true) {
    if (added) this.canvas.remove(repr);
    repr.destroy();
    for (const release of entry.releases) release();
  }

  private assertLive() {
    if (this.disposed) throw new Error('VisualSet: disposed');
  }
}

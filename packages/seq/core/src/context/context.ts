/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue, Node } from '@molstar/core/repr/tree';
import { type AssetManager, DefaultAssetManager } from '@molstar/core/util/assets';
import type { AssetReader, FormatRuntime } from '@molstar/seqx-core/data';
import { runtimesByFormat } from '@molstar/seqx-core/data/channels';
import { checkSeqSpec, type SpecBounds, specBounds } from '@molstar/seqx-core/lang';
import {
  type QueryExtensionEntry,
  type QueryExtensions,
  type QueryFormEntry,
  QueryFormKinds,
  type QueryOperatorEntry,
  type QueryPropertyEntry,
  queryExtensions,
} from '@molstar/seqx-core/query';
import {
  type RenderConfig,
  type RenderConfigOverrides,
  exportSvg,
  type LayerType,
  type NavigatorPicture,
  type RenderSvgOptions,
  type RenderSvgResult,
  renderConfig,
  type SectionLayoutType,
  type SeqTokens,
  type SizedTrack,
  type TextMetrics,
  type TrackSizing,
  type TrackSizingHook,
  builtInTokens,
  bundledTextMetrics,
  createTextMetrics,
  NoAssets,
  refuseMuteTypes,
} from '@molstar/seqx-core/render';
import { type ResolveConfig, resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import type { SchemeEntry, SchemeRegistry } from '@molstar/seqx-core/theme';
import {
  BuiltInLayerKinds,
  type LayerKinds,
  type ResolvedSection,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import { type AssetPolicy, HeadlessAssetPolicy, policyAssets } from './assets.js';
import { checkExtensions, type ExtensionTables, type SeqExtension } from './extension.js';

// The headless composition (_spec/README.md "Context"): what a render needs, listed in one place: the registry (the
// format runtimes a spec's sources are read with), the configuration and the ports. No global registry, nothing
// registered on import; the viewer composes the same pieces.

/**
 * What a registry holds in this version (_spec/README.md "Context", Registry entries): a format runtime, read by the
 * sources whose `format` it names; a layer type, drawn for the layers of its name (the canvas and the SVG alike); a
 * section layout, placing the sections that name it; a scheme the themes name; a query form of a role and a query
 * property, namespaced (`<id>:<name>`), what every query of the composition reads (`SeqContext.query`).
 */
export type SeqRegistryEntry =
  | { readonly kind: 'format-runtime'; readonly runtime: FormatRuntime }
  | { readonly kind: 'layer-type'; readonly type: LayerType<any, any> }
  | { readonly kind: 'layout'; readonly layout: SectionLayoutType<any> }
  | { readonly kind: 'track-sizing'; readonly sizing: TrackSizing }
  | { readonly kind: 'scheme'; readonly scheme: SchemeEntry }
  | { readonly kind: 'navigator-picture'; readonly picture: NavigatorPicture }
  | QueryFormEntry
  | QueryOperatorEntry
  | QueryPropertyEntry;

/** The registry's query forms and properties, composed and checked (names namespaced, no built-in shadowed). */
export function queryExtensionsOf(registry: readonly SeqRegistryEntry[]): QueryExtensions {
  return queryExtensions(
    registry.filter((e): e is QueryExtensionEntry => e.kind === 'query-form' || e.kind === 'query-property'),
    QueryFormKinds,
  );
}

/** A setting: its key, and the type of its value. Items are pairs `[item, value]` in a spec's `config`. */
export interface SeqConfigItem<T = unknown> {
  readonly key: string;
  /** Never set: carries the value's type. */
  readonly type?: T;
}

export function configItem<T>(key: string): SeqConfigItem<T> {
  return { key };
}

/** The core's settings. */
export const SeqCoreConfig = {
  /** The built-in token table: `light` by default. */
  theme: configItem<'light' | 'dark'>('theme'),
  /** The layout's numbers over `DefaultRenderConfig`. */
  render: configItem<RenderConfigOverrides>('render'),
  /** The resolve passes' bounds over `DefaultResolveConfig`. */
  resolve: configItem<Partial<ResolveConfig>>('resolve'),
  /** The bounds on what a spec sizes over `DefaultSpecBounds` (size, depth, query depth, column length, issues kept). */
  bounds: configItem<Partial<SpecBounds>>('bounds'),
} as const;

/** A port's adapter and whether the context that asked for it disposes it (true unless the host keeps it). */
export interface PortResult {
  /** True by default: the adapter goes with its owner; a host that passes its own sets false and keeps it. */
  readonly disposeWithViewer?: boolean;
}

/** The default adapter: the application's `AssetManager` behind a policy. */
export interface ManagedAssetsPort extends PortResult {
  readonly manager: AssetManager;
  /** What may be read through the manager; the headless default (nothing) when absent. */
  readonly policy?: AssetPolicy;
}

/** Another adapter, a second one of the port: a reader of its own (the test adapter that refuses the network). */
export interface ReaderAssetsPort extends PortResult {
  readonly reader: AssetReader;
  dispose?(): void;
}

export type AssetsPort = ManagedAssetsPort | ReaderAssetsPort;

export interface TextMetricsPort extends PortResult {
  /** The width of `text` set in `font`, CSS px. */
  measure(font: string, text: string): number;
  /**
   * The cap height of `font`, CSS px (a canvas's `actualBoundingBoxAscent` of `H`): what makes a letter's size match a
   * label's (`letters.capRatio`). Absent: the faces are taken to have the same cap height per em.
   */
  capHeight?(font: string): number;
  dispose?(): void;
}

/**
 * The host's last word on track heights (_spec/README.md "Render", Track heights): each section's tracks with their
 * content heights and the named sizing's answer, returning the heights (an application fitting the viewer into a fixed
 * panel or its own row grid).
 */
export interface TrackSizingPort extends PortResult {
  heights(section: ResolvedSection, tracks: readonly SizedTrack[], named: readonly number[]): readonly number[];
  dispose?(): void;
}

/** The ports of a composition, each a factory called once. Unset: the default adapter, owned by the owner. */
export interface SeqPorts {
  readonly assets?: () => AssetsPort;
  readonly textMetrics?: () => TextMetricsPort;
  readonly trackSizing?: () => TrackSizingPort;
}

export interface SeqContextSpec {
  readonly registry?: readonly SeqRegistryEntry[];
  /**
   * Extensions (_spec/README.md "Context", Extensions): each one's entries registered with the registry's, the
   * descriptors whole (checked: namespaces, versions, conflicts; nothing registered on an error).
   */
  readonly extensions?: readonly SeqExtension[];
  /** Settings as `[item, value]` pairs. */
  readonly config?: readonly (readonly [SeqConfigItem<any>, unknown])[];
  readonly ports?: SeqPorts;
}

/**
 * The registered runtimes, keyed by the formats they read (one adapter a format: `withAdapter` replaces one): the same
 * object again is no change, a format two of them read is an error listing every conflict. Not by name: a runtime
 * narrowed by `withAdapter` keeps its name beside an adapter that may share it.
 */
export function runtimesOf(registry: readonly SeqRegistryEntry[]): FormatRuntime[] {
  const runtimes = [...new Set(registry.flatMap((e) => (e.kind === 'format-runtime' ? [e.runtime] : [])))];
  const { conflicts } = runtimesByFormat(runtimes);
  if (conflicts.size)
    throw new Error(
      `The registry holds several format runtimes for one format: ${[...conflicts].map(([f, names]) => `${f} (${names.join(', ')})`).join('; ')}`,
    );
  return runtimes;
}

/** The registered layer types by name, as `runtimesOf`: the same object again is no change, another an error. */
function layerTypesOf(registry: readonly SeqRegistryEntry[]): Readonly<Record<string, LayerType<any, any>>> {
  const types = byName(registry, 'layer-type', 'layer types', (e) => e.type);
  refuseMuteTypes(types.values());
  return Object.fromEntries(types);
}

/** The registered track sizings by name, as `runtimesOf` (the built-ins need no entry). */
function sizingsOf(registry: readonly SeqRegistryEntry[]): Readonly<Record<string, TrackSizing>> {
  return Object.fromEntries(byName(registry, 'track-sizing', 'track sizings', (e) => e.sizing));
}

/** The registered section layouts by name, as `runtimesOf`. */
function layoutsOf(registry: readonly SeqRegistryEntry[]): Readonly<Record<string, SectionLayoutType<any>>> {
  return Object.fromEntries(byName(registry, 'layout', 'layouts', (e) => e.layout));
}

/** One kind's entries by their provider's name; every conflict listed in one error. */
function byName<K extends SeqRegistryEntry['kind'], T extends { readonly name: string }>(
  registry: readonly SeqRegistryEntry[],
  kind: K,
  what: string,
  of: (e: Extract<SeqRegistryEntry, { kind: K }>) => T,
): Map<string, T> {
  const out = new Map<string, T>();
  const conflicts: string[] = [];
  for (const e of registry) {
    if (e.kind !== kind) continue;
    const p = of(e as Extract<SeqRegistryEntry, { kind: K }>);
    const had = out.get(p.name);
    if (had && had !== p) conflicts.push(p.name);
    else out.set(p.name, p);
  }
  if (conflicts.length) throw new Error(`The registry holds different ${what} under one name: ${conflicts.join(', ')}`);
  return out;
}

/** The value of a config item in the pairs (the last), else `undefined`. */
export function configValue<T>(
  config: readonly (readonly [SeqConfigItem<any>, unknown])[] | undefined,
  item: SeqConfigItem<T>,
): T | undefined {
  let out: T | undefined;
  for (const [i, v] of config ?? []) if (i.key === item.key) out = v as T;
  return out;
}

/** The registry's tables, composed from its entries and its extensions': what `_extend` swaps whole. */
type Tables = Pick<
  SeqContext,
  'runtimes' | 'types' | 'layouts' | 'sizings' | 'pictures' | 'query' | 'schemes' | 'kinds' | 'extensions'
>;

/**
 * The tables of a registry and its checked extensions: every conflict across their entries found first, one error,
 * nothing changed; then what each extension adds (`defined.compose`: its layer kinds, its schemes).
 */
function tablesOf(registry: readonly SeqRegistryEntry[], extensions: readonly SeqExtension[]): Tables {
  const all = [...registry, ...extensions.flatMap((e) => e.registry)];
  let added: ExtensionTables = { kinds: BuiltInLayerKinds };
  for (const e of extensions) added = e.defined!.compose(added);
  // the registered navigator pictures by name, as `runtimesOf`
  const pictures = Object.fromEntries(byName(all, 'navigator-picture', 'navigator pictures', (e) => e.picture));
  return {
    runtimes: runtimesOf(all),
    types: layerTypesOf(all),
    layouts: layoutsOf(all),
    sizings: sizingsOf(all),
    pictures,
    query: queryExtensionsOf(all),
    schemes: undefined,
    ...added,
    // the pictures a navigator's `content` may name: the registry's (the axis always drawn)
    kinds: { ...added.kinds, pictures: Object.keys(pictures) },
    extensions,
  };
}

export interface SeqContext {
  readonly runtimes: readonly FormatRuntime[];
  /** The registered layer types by name: what the canvas and `renderSvg` draw (a layer of another type is said). */
  readonly types: Readonly<Record<string, LayerType<any, any>>>;
  /** The registered section layouts by name (a section naming another is said, and laid out as tracks). */
  readonly layouts: Readonly<Record<string, SectionLayoutType<any>>>;
  /** The registered track sizings by name (a section's `heights` naming another is said; its tracks fit their content). */
  readonly sizings: Readonly<Record<string, TrackSizing>>;
  /** The registered navigator pictures by name (a `content` naming another is said, and drawn as the axis). */
  readonly pictures: Readonly<Record<string, NavigatorPicture>>;
  /** The registered query forms and properties: what `evaluateQuery` is given as `extensions`. */
  readonly query: QueryExtensions;
  /** The schemes the themes read, when the registry adds any to the built-ins (undefined: the built-ins). */
  readonly schemes?: SchemeRegistry;
  /** The layer kinds a view is resolved with: the built-ins', and the extensions' (_spec/README.md "Resolved view"). */
  readonly kinds: LayerKinds;
  /** The extensions composed, each once. */
  readonly extensions: readonly SeqExtension[];
  /** The tables' version: moves each time `_extend` swaps them (an extension registered late, or its undo). */
  readonly version: number;
  /** The `trackSizing` port's answer, when one was bound. */
  readonly trackSizing?: TrackSizingHook;
  readonly renderConfig: RenderConfig;
  readonly resolveConfig: Partial<ResolveConfig>;
  /** The bounds the spec's check holds (`SeqCoreConfig.bounds` over the defaults). */
  readonly bounds: SpecBounds;
  readonly tokens: SeqTokens;
  /** The assets port: what a source's location is read through. */
  readonly assets: AssetReader;
  /** The text metrics port. */
  readonly text: TextMetrics;
  /** The spec's issues (the language's checks). */
  check(spec: Node): Issue[];
  /** The static model of a spec with the registered runtimes. */
  resolve(spec: Node): StaticModel;
  /** A view of a model: the address named, else the first. */
  view(model: StaticModel, address?: string): ResolvedView | undefined;
  /** A spec as an SVG string once its data has settled, through this context's ports. */
  renderSvg(spec: Node, options: Pick<RenderSvgOptions, 'size'> & Partial<RenderSvgOptions>): Promise<RenderSvgResult>;
  /**
   * The extensions composed made these, checked as at composition (nothing changed on an error), the tables and their
   * version moved: what `SeqViewer.extend` and its undo call (a headless host composes another context).
   * @internal
   */
  _extend(extensions: readonly SeqExtension[]): void;
  /** Disposes the ports this context made or was told to dispose; a host's are left. */
  retire(): void;
}

/**
 * The config with `letters.capRatio` measured: the label face's cap height over the mono face's per em, so a letter
 * matches a label's cap height (the host's own `capRatio` is kept). Unchanged where the port measures no cap heights.
 */
function withCapRatio(config: RenderConfig, port: TextMetricsPort): RenderConfig {
  if (!port.capHeight || config.letters.capRatio !== 1) return config;
  const label = port.capHeight('100px sans');
  const mono = port.capHeight('100px mono');
  if (!(label > 0) || !(mono > 0)) return config;
  return { ...config, letters: { ...config.letters, capRatio: label / mono } };
}

/**
 * The headless context of a composition as listed (_spec/README.md "Context"): the registry, the configuration and the
 * ports a render needs, and the resolve passes and `renderSvg` through them, with nothing the registry does not list (no
 * layer type, no layout, no navigators: `createSeqContext` lists the built-ins). Ports unset take their defaults (no
 * reads, the bundled text metrics), owned and disposed by `retire`.
 */
export function composeSeqContext(spec: SeqContextSpec = {}): SeqContext {
  const registry = spec.registry ?? [];
  let t = tablesOf(registry, checkExtensions(spec.extensions ?? []));
  const theme = configValue(spec.config, SeqCoreConfig.theme) ?? 'light';
  let renderCfg = renderConfig(configValue(spec.config, SeqCoreConfig.render));
  const resolveCfg = configValue(spec.config, SeqCoreConfig.resolve) ?? {};
  const bounds = specBounds(configValue(spec.config, SeqCoreConfig.bounds));
  const disposers: (() => void)[] = [];

  // a factory that throws leaves nothing made behind: the ports already made are disposed, last first
  let assets: AssetReader = NoAssets;
  let text: TextMetrics = undefined as never;
  let trackSizing: TrackSizingHook | undefined;
  try {
    if (spec.ports?.assets) {
      const port = spec.ports.assets();
      if ('reader' in port) {
        assets = port.reader;
        if (port.disposeWithViewer !== false && port.dispose) disposers.push(() => port.dispose!());
      } else {
        assets = policyAssets(port.manager, port.policy ?? HeadlessAssetPolicy);
        if (port.disposeWithViewer !== false) disposers.push(() => port.manager.dispose());
      }
    }

    if (spec.ports?.textMetrics) {
      const port = spec.ports.textMetrics();
      text = createTextMetrics((font) => (s) => port.measure(font, s));
      renderCfg = withCapRatio(renderCfg, port);
      if (port.disposeWithViewer !== false && port.dispose) disposers.push(() => port.dispose!());
    } else text = bundledTextMetrics();

    if (spec.ports?.trackSizing) {
      const port = spec.ports.trackSizing();
      trackSizing = (section, tracks, named) => port.heights(section, tracks, named);
      if (port.disposeWithViewer !== false && port.dispose) disposers.push(() => port.dispose!());
    }
  } catch (e) {
    for (const d of disposers.reverse()) {
      try {
        d();
      } catch {
        // the factory's failure is the one reported
      }
    }
    throw e;
  }

  const tokens = builtInTokens(theme);
  let retired = false;
  const context: SeqContext = {
    ...t,
    version: 0,
    ...(trackSizing ? { trackSizing } : {}),
    renderConfig: renderCfg,
    resolveConfig: resolveCfg,
    bounds,
    tokens,
    assets,
    text,
    check: (s) => [...checkSeqSpec(s, bounds), ...t.extensions.flatMap((e) => e.check?.(s) ?? [])],
    resolve: (s) => resolveSpec(s, { runtimes: t.runtimes, config: resolveCfg }),
    view: (model, address) => (address ? resolveView(model, address, t.kinds) : resolveViews(model, t.kinds)[0]),
    renderSvg: (s, options) =>
      exportSvg(s, {
        runtimes: t.runtimes,
        types: t.types,
        layouts: t.layouts,
        sizings: t.sizings,
        pictures: t.pictures,
        kinds: t.kinds,
        ...(t.schemes ? { registry: { schemes: t.schemes } } : {}),
        ...(trackSizing ? { trackSizing } : {}),
        theme,
        tokens,
        bounds,
        text,
        assets,
        config: renderCfg,
        resolve: resolveCfg,
        ...options,
      }),
    _extend(next) {
      t = tablesOf(registry, checkExtensions(next));
      Object.assign(context, t, { version: context.version + 1 });
    },
    retire() {
      if (retired) return;
      retired = true;
      for (const d of disposers) d();
      disposers.length = 0;
    },
  };
  return context;
}

/** A default assets port: a new `DefaultAssetManager` under a policy, disposed with its owner. */
export function defaultAssetsPort(policy: AssetPolicy = HeadlessAssetPolicy): ManagedAssetsPort {
  return { manager: new DefaultAssetManager(), policy };
}

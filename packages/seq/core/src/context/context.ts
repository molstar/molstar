/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue, Node } from '@molstar/core/repr/tree';
import { type AssetManager, DefaultAssetManager } from '@molstar/core/util/assets';
import type { AssetReader, FormatRuntime } from '@molstar/seqx-core/data';
import { checkSeqSpec, type SpecBounds, specBounds } from '@molstar/seqx-core/lang';
import {
  type RenderConfig,
  type RenderConfigOverrides,
  type RenderSvgOptions,
  type RenderSvgResult,
  renderConfig,
  renderSvg,
  type SeqTokens,
  type TextMetrics,
  builtInTokens,
  bundledTextMetrics,
  createTextMetrics,
  NoAssets,
} from '@molstar/seqx-core/render';
import { type ResolveConfig, resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import { type ResolvedView, resolveView, resolveViews } from '@molstar/seqx-core/view';
import { type AssetPolicy, HeadlessAssetPolicy, policyAssets } from './assets.js';

// The headless composition (_spec/README.md "Context"): what a render needs, listed in one place: the registry (the
// format runtimes a spec's sources are read with), the configuration and the ports. No global registry, nothing
// registered on import; the viewer composes the same pieces.

/** What a registry holds in this version: a format runtime, read by the sources whose `format` it names. */
export type SeqRegistryEntry = { readonly kind: 'format-runtime'; readonly runtime: FormatRuntime };

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
  dispose?(): void;
}

/** A structure runtime binds structure systems to a structure viewer: no adapter in this version, unset. */
export interface StructurePort extends PortResult {
  readonly runtime: unknown;
  dispose?(): void;
}

/** The ports of a composition, each a factory called once. Unset: the default adapter, owned by the owner. */
export interface SeqPorts {
  readonly assets?: () => AssetsPort;
  readonly textMetrics?: () => TextMetricsPort;
  readonly structure?: () => StructurePort;
}

export interface SeqContextSpec {
  readonly registry?: readonly SeqRegistryEntry[];
  /** Settings as `[item, value]` pairs. */
  readonly config?: readonly (readonly [SeqConfigItem<any>, unknown])[];
  readonly ports?: SeqPorts;
}

/** The registered runtimes by name: the same object again is no change, another under one name is an error (all listed). */
export function runtimesOf(registry: readonly SeqRegistryEntry[]): FormatRuntime[] {
  const byName = new Map<string, FormatRuntime>();
  const conflicts: string[] = [];
  for (const e of registry) {
    if (e.kind !== 'format-runtime') continue;
    const had = byName.get(e.runtime.name);
    if (had && had !== e.runtime) conflicts.push(e.runtime.name);
    else byName.set(e.runtime.name, e.runtime);
  }
  if (conflicts.length)
    throw new Error(`The registry holds different format runtimes under one name: ${conflicts.join(', ')}`);
  return [...byName.values()];
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

export interface SeqContext {
  readonly runtimes: readonly FormatRuntime[];
  readonly renderConfig: RenderConfig;
  readonly resolveConfig: Partial<ResolveConfig>;
  /** The bounds the spec's check holds (`SeqCoreConfig.bounds` over the defaults). */
  readonly bounds: SpecBounds;
  readonly tokens: SeqTokens;
  /** The assets port: what a source's location is read through. */
  readonly assets: AssetReader;
  /** The text metrics port. */
  readonly text: TextMetrics;
  /** The structure port's runtime, when one was bound. */
  readonly structure: unknown;
  /** The spec's issues (the language's checks). */
  check(spec: Node): Issue[];
  /** The static model of a spec with the registered runtimes. */
  resolve(spec: Node): StaticModel;
  /** A view of a model: the address named, else the first. */
  view(model: StaticModel, address?: string): ResolvedView | undefined;
  /** A spec as an SVG string once its data has settled, through this context's ports. */
  renderSvg(spec: Node, options: Pick<RenderSvgOptions, 'size'> & Partial<RenderSvgOptions>): Promise<RenderSvgResult>;
  /** Disposes the ports this context made or was told to dispose; a host's are left. */
  retire(): void;
}

/**
 * The headless context of a composition: the registry, the configuration and the ports a render needs, and the resolve
 * passes and `renderSvg` through them. Ports unset take their defaults (no reads, the bundled text metrics), owned and
 * disposed by `retire`.
 */
export function createSeqContext(spec: SeqContextSpec = {}): SeqContext {
  const runtimes = runtimesOf(spec.registry ?? []);
  const theme = configValue(spec.config, SeqCoreConfig.theme) ?? 'light';
  const renderCfg = renderConfig(configValue(spec.config, SeqCoreConfig.render));
  const resolveCfg = configValue(spec.config, SeqCoreConfig.resolve) ?? {};
  const bounds = specBounds(configValue(spec.config, SeqCoreConfig.bounds));
  const disposers: (() => void)[] = [];

  // a factory that throws leaves nothing made behind: the ports already made are disposed, last first
  let assets: AssetReader = NoAssets;
  let structure: unknown;
  let text: TextMetrics = undefined as never;
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
      if (port.disposeWithViewer !== false && port.dispose) disposers.push(() => port.dispose!());
    } else text = bundledTextMetrics();

    if (spec.ports?.structure) {
      const port = spec.ports.structure();
      structure = port.runtime;
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
    runtimes,
    renderConfig: renderCfg,
    resolveConfig: resolveCfg,
    bounds,
    tokens,
    assets,
    text,
    structure,
    check: (s) => checkSeqSpec(s, bounds),
    resolve: (s) => resolveSpec(s, { runtimes, config: resolveCfg }),
    view: (model, address) => (address ? resolveView(model, address) : resolveViews(model)[0]),
    renderSvg: (s, options) =>
      renderSvg(s, {
        runtimes,
        theme,
        tokens,
        text,
        assets,
        config: renderCfg,
        resolve: resolveCfg,
        ...options,
      }),
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

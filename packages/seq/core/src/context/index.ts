/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/context: the headless composition (_spec/README.md "Context"): the registry, the configuration and
// the ports a render needs, the policy wrapper of the assets port, and `createSeqContext`. Named re-exports only.

export { withAdapter } from './adapters.js';
export {
  type AssetPolicy,
  AssetRefused,
  BrowserAssetPolicy,
  checkLocation,
  HeadlessAssetPolicy,
  headlessAssetPolicy,
  policyAssets,
  rootUrl,
  utf8Length,
} from './assets.js';
export {
  type AssetsPort,
  configItem,
  composeSeqContext,
  configValue,
  defaultAssetsPort,
  type ManagedAssetsPort,
  type PortResult,
  queryExtensionsOf,
  type ReaderAssetsPort,
  runtimesOf,
  type SeqConfigItem,
  type SeqContext,
  type SeqContextSpec,
  SeqCoreConfig,
  type SeqPorts,
  type SeqRegistryEntry,
  type TextMetricsPort,
  type TrackSizingPort,
} from './context.js';
export type { ExtensionTables, SeqExtension } from './extension.js';

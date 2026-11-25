/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/extension (this file, through the package's `./*` export): the one entry an extension imports
// (_spec/README.md "Context", Extensions): the descriptor and `defineExtension`, the layer-type, layout and navigator-picture
// contracts and what they are given (the window's geometry and drawing for one that draws a window of its own), the
// format runtime's, the scheme's. What it re-exports is the contract a review reads before
// a break; everything else is internal and may change in a 0.x minor. Named re-exports only.

export {
  configItem,
  settingOf,
  type SeqConfigItem,
  type SeqExtension,
  type SeqRegistryEntry,
} from '@molstar/seqx-core/context';
export {
  coarsenBins,
  type FormatRuntime,
  type PairsBins,
  type PairsList,
  type PairsTable,
} from '@molstar/seqx-core/data';
export { defineExtension, requiresCovers } from '@molstar/seqx-core/context/define';
export {
  type Address,
  type NavigatorParams,
  type Params,
  type RefValue,
  type SectionBuilder,
  theme,
  type TrackBuilder,
  type ViewBuilder,
} from '@molstar/seqx-core/lang';
export { type Key, type Keys, type Location, rangesLocation } from '@molstar/seqx-core/query';
export type { Overlay, PairsData, StaticModel } from '@molstar/seqx-core/resolve';
export {
  type Box,
  channelLegends,
  type ContentContext,
  type Emphasis,
  type Hit,
  type HitMark,
  type ItemDescription,
  type LayerPlace,
  type LayerType,
  type SquareCells,
  type SquareLattice,
  drawNavigatorWindow,
  floorExtent,
  type LayoutContext,
  type NavigatorPicture,
  type NavigatorPictureContext,
  type NavigatorPictureTracks,
  navigatorPart,
  NavWindowFlags,
  type PairsLevels,
  pairsLevels,
  type PairsRead,
  PairsReader,
  type PairsReadOptions,
  type PairsSource,
  type PassName,
  type PassWriter,
  type Placement,
  PlacePoint,
  PolarWarp,
  TakenLabels,
  type TakenOwner,
  type SectionLayoutType,
  type StoreReader,
  sectionId,
  type ThemeContext,
  TranslateWarp,
} from '@molstar/seqx-core/render';
export type { Legend, SchemeEntry, SchemeSubject } from '@molstar/seqx-core/theme';
export type { NavigatorWindow } from '@molstar/seqx-core/render/navigator';
export {
  type BlockPlace,
  type LayerShape,
  placesOf,
  type ResolvedLayer,
  type ResolvedSection,
} from '@molstar/seqx-core/view';

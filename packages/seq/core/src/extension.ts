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
  type SeqConfigItem,
  type SeqExtension,
  type SeqRegistryEntry,
} from '@molstar/seqx-core/context';
export type { FormatRuntime } from '@molstar/seqx-core/data';
export { defineExtension, requiresCovers } from '@molstar/seqx-core/context/define';
export type { Key, Keys, Location } from '@molstar/seqx-core/query';
export type { Overlay, StaticModel } from '@molstar/seqx-core/resolve';
export {
  type Box,
  type ContentContext,
  type Emphasis,
  type Hit,
  type HitMark,
  type ItemDescription,
  type LayerType,
  drawNavigatorWindow,
  floorExtent,
  type LayoutContext,
  type NavigatorPicture,
  type NavigatorPictureContext,
  type NavigatorPictureTracks,
  navigatorPart,
  NavWindowFlags,
  type PassName,
  type PassWriter,
  type Placement,
  PolarWarp,
  type SectionLayoutType,
  type StoreReader,
  sectionId,
  type ThemeContext,
  TranslateWarp,
} from '@molstar/seqx-core/render';
export type { Legend, SchemeEntry, SchemeSubject } from '@molstar/seqx-core/theme';
export type { NavigatorWindow } from '@molstar/seqx-core/render/navigator';
export type { LayerShape, ResolvedLayer, ResolvedSection } from '@molstar/seqx-core/view';

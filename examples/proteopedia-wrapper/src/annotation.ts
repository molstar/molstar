/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomElementProperty } from '@molstar/plugin/properties/custom-element-property';
import { Model, type ElementIndex, type ResidueIndex } from '@molstar/model/model/structure';
import { Color } from '@molstar/core/util/color';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { Asset } from '@molstar/core/util/assets';

const EvolutionaryConservationPalette: Color[] = [
  [255, 255, 129], // insufficient
  [160, 37, 96], // 9
  [240, 125, 171],
  [250, 201, 222],
  [252, 237, 244],
  [255, 255, 255],
  [234, 255, 255],
  [215, 255, 255],
  [140, 255, 255],
  [16, 200, 209], // 1
]
  .reverse()
  .map(([r, g, b]) => Color.fromRgb(r, g, b));
const EvolutionaryConservationDefaultColor = Color(0x999999);

export const EvolutionaryConservation = CustomElementProperty.create<number>({
  name: 'proteopedia-wrapper-evolutionary-conservation',
  label: 'Evolutionary Conservation',
  type: 'static',
  async getData(model: Model, ctx: CustomProperty.Context) {
    const id = model.entryId.toLowerCase();
    const url = Asset.getUrlAsset(ctx.assetManager, `https://proteopedia.org/cgi-bin/cnsrf?${id}`);
    const json = await ctx.assetManager.resolve(url, 'json').runInContext(ctx.runtime);
    const annotations = json.data?.residueAnnotations || [];

    const conservationMap = new Map<string, number>();

    for (const e of annotations) {
      for (const r of e.ids) {
        conservationMap.set(r, e.annotation);
      }
    }

    const map = new Map<ElementIndex, number>();

    const { _rowCount: residueCount } = model.atomicHierarchy.residues;
    const { offsets: residueOffsets } = model.atomicHierarchy.residueAtomSegments;
    const chainIndex = model.atomicHierarchy.chainAtomSegments.index;

    for (let rI = 0 as ResidueIndex; rI < residueCount; rI++) {
      const cI = chainIndex[residueOffsets[rI]];
      const key = `${model.atomicHierarchy.chains.auth_asym_id.value(cI)} ${model.atomicHierarchy.residues.auth_seq_id.value(rI)}`;
      if (!conservationMap.has(key)) continue;
      const ann = conservationMap.get(key)!;
      for (let aI = residueOffsets[rI]; aI < residueOffsets[rI + 1]; aI++) {
        map.set(aI, ann);
      }
    }

    return { value: map, assets: [json] };
  },
  coloring: {
    getColor(e: number) {
      if (e < 1 || e > 10) return EvolutionaryConservationDefaultColor;
      return EvolutionaryConservationPalette[e - 1];
    },
    defaultColor: EvolutionaryConservationDefaultColor,
  },
  getLabel(e) {
    if (e === 10) return `Evolutionary Conservation: Insufficient Data`;
    return e ? `Evolutionary Conservation: ${e}` : void 0;
  },
});

/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '@molstar/model/model/structure';
import type { VolumeServerInfo } from './model.js';
import { MmcifFormat } from '@molstar/model/formats/structure/mmcif-format';

export function getStreamingMethod(s?: Structure, defaultKind: VolumeServerInfo.Kind = 'x-ray'): VolumeServerInfo.Kind {
  if (!s) return defaultKind;

  const model = s.models[0];
  if (!MmcifFormat.is(model.sourceData)) return defaultKind;

  // Prefer EMDB entries over structure-factors (SF) e.g. for 'ELECTRON CRYSTALLOGRAPHY' entries
  // like 6AXZ or 6KJ3 for which EMDB entries are available but map calculation from SF is hard.
  if (Model.hasEmMap(model)) return 'em';
  if (Model.hasXrayMap(model)) return 'x-ray';

  // Fallbacks based on experimental method
  if (Model.isFromEm(model)) return 'em';
  if (Model.isFromXray(model)) return 'x-ray';
  return defaultKind;
}

/** Returns EMD ID when available, otherwise falls back to PDB ID */
export function getEmIds(model: Model): string[] {
  const ids: string[] = [];
  if (!MmcifFormat.is(model.sourceData)) return [model.entryId];

  const { db_id, db_name, content_type } = model.sourceData.data.db.pdbx_database_related;
  if (!db_name.isDefined) return [model.entryId];

  for (let i = 0, il = db_name.rowCount; i < il; ++i) {
    if (db_name.value(i).toUpperCase() === 'EMDB' && content_type.value(i) === 'associated EM volume') {
      ids.push(db_id.value(i));
    }
  }

  return ids;
}

export function getXrayIds(model: Model): string[] {
  return [model.entryId];
}

export function getIds(method: VolumeServerInfo.Kind, s?: Structure): string[] {
  if (!s || !s.models.length) return [];
  const model = s.models[0];
  switch (method) {
    case 'em':
      return getEmIds(model);
    case 'x-ray':
      return getXrayIds(model);
  }
}

/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { RuntimeContext } from '@molstar/core/task';
import { PluginConfig } from '@molstar/plugin/config';

/**
 * EMDB/PDBe lookups shared by volume streaming and the volume format providers. Kept out of the volume streaming
 * behavior modules so base modules can use them without loading the behavior.
 */

export async function getContourLevel(
  provider: 'emdb' | 'pdbe',
  plugin: PluginContext,
  taskCtx: RuntimeContext,
  emdbId: string,
) {
  switch (provider) {
    case 'emdb':
      return getContourLevelEmdb(plugin, taskCtx, emdbId);
    case 'pdbe':
      return getContourLevelPdbe(plugin, taskCtx, emdbId);
  }
}

export async function getContourLevelEmdb(plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
  const emdbHeaderServer = plugin.config.get(PluginConfig.VolumeStreaming.EmdbHeaderServer);
  const header = await plugin
    .fetch({ url: `${emdbHeaderServer}/${emdbId.toUpperCase()}/header/${emdbId.toLowerCase()}.xml`, type: 'xml' })
    .runInContext(taskCtx);
  const map = header.getElementsByTagName('map')[0];
  const contours = map.getElementsByTagName('contour');

  let primaryContour = contours[0];
  for (let i = 1; i < contours.length; i++) {
    if (contours[i].getAttribute('primary') === 'true') {
      primaryContour = contours[i];
      break;
    }
  }
  const contourLevel = parseFloat(primaryContour.getElementsByTagName('level')[0].textContent!);
  return contourLevel;
}

export async function getContourLevelPdbe(plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
  // TODO: parametrize URL in plugin settings?
  emdbId = emdbId.toUpperCase();
  const header = await plugin
    .fetch({ url: `https://www.ebi.ac.uk/emdb/api/entry/map/${emdbId}`, type: 'json' })
    .runInContext(taskCtx);
  const contours = header?.map?.contour_list?.contour;

  if (!contours || contours.length === 0) {
    // try fallback to the old API
    return getContourLevelPdbeLegacy(plugin, taskCtx, emdbId);
  }

  return contours.find((c: any) => c.primary)?.level ?? contours[0].level;
}

async function getContourLevelPdbeLegacy(plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
  // TODO: parametrize URL in plugin settings?
  emdbId = emdbId.toUpperCase();
  const header = await plugin
    .fetch({ url: `https://www.ebi.ac.uk/pdbe/api/emdb/entry/map/${emdbId}`, type: 'json' })
    .runInContext(taskCtx);
  const emdbEntry = header?.[emdbId];
  let contourLevel: number | undefined = void 0;
  if (emdbEntry?.[0]?.map?.contour_level?.value !== void 0) {
    contourLevel = +emdbEntry[0].map.contour_level.value;
  }

  return contourLevel;
}

export async function getEmdbIds(plugin: PluginContext, taskCtx: RuntimeContext, pdbId: string) {
  // TODO: parametrize to a differnt URL? in plugin settings perhaps
  const summary = await plugin
    .fetch({ url: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/${pdbId}`, type: 'json' })
    .runInContext(taskCtx);

  const summaryEntry = summary?.[pdbId];
  const emdbIds: string[] = [];
  if (summaryEntry?.[0]?.related_structures) {
    const emdb = summaryEntry[0].related_structures.filter(
      (s: any) => s.resource === 'EMDB' && s.relationship === 'associated EM volume',
    );
    if (!emdb.length) {
      throw new Error(`No related EMDB entry found for '${pdbId}'.`);
    }
    emdbIds.push(...emdb.map((e: { accession: string }) => e.accession));
  } else {
    throw new Error(`No related EMDB entry found for '${pdbId}'.`);
  }

  return emdbIds;
}

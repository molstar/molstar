/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-molstar-poc/plugin`: a Mol* plugin (a canvas, no Mol* UI) as the link's 3D view
// (seqspec:libraries#the-mol-adapter), and the structure an MVS load left in it. Loci are carried onto that structure
// (`remap`): Mol*'s hover reports a component's substructure, whose positions differ from the root's.

import { MolViewSpecBehavior } from "molstar/lib/extensions/mvs/behavior";
import { loadMVS } from "molstar/lib/extensions/mvs/load";
import type { MVSData_State } from "molstar/lib/extensions/mvs/mvs-data";
import { type Structure, StructureElement } from "molstar/lib/mol-model/structure";
import { PluginContext } from "molstar/lib/mol-plugin/context";
import { DefaultPluginSpec } from "molstar/lib/mol-plugin/spec";
import { Asset } from "molstar/lib/mol-util/assets";
import { Color } from "molstar/lib/mol-util/color/color";
import type { StructureView3D } from "./viewer.ts";

export type { PluginContext };

/**
 * A Mol* plugin mounted on `host` as a canvas only (no Mol* UI): MVS's behavior, no axes, clicks selecting, and
 * `files` registered as Mol*'s assets at their URLs (an MVS download of such a URL reads them: offline fixtures).
 */
export async function mountPlugin(
  host: HTMLElement,
  o: { files?: ReadonlyMap<string, Uint8Array>; background?: number } = {},
): Promise<PluginContext> {
  const spec = DefaultPluginSpec();
  spec.behaviors = [...spec.behaviors, MolViewSpecBehavior];
  spec.canvas3d = {
    ...spec.canvas3d,
    camera: { ...spec.canvas3d?.camera, helper: { axes: { name: "off", params: {} } } },
  };
  const plugin = new PluginContext(spec);
  await plugin.init();
  await plugin.mountAsync(host);
  for (const [url, bytes] of o.files ?? []) {
    const name = url.split("/").pop() ?? url;
    plugin.managers.asset.set(Asset.Url(url), new File([bytes as Uint8Array<ArrayBuffer>], name), { isStatic: true });
  }
  plugin.behaviors.interaction.selectionMode.next(true);
  if (o.background !== undefined) setBackground(plugin, o.background);
  return plugin;
}

/** The canvas's background, 0xRRGGBB (a host's light or dark surface). */
export function setBackground(plugin: PluginContext, color: number): void {
  plugin.canvas3d?.setProps({ renderer: { backgroundColor: Color(color) } });
}

/**
 * Loads an MVS state (replacing what was shown) and gives the structure it made, for the link; null: none.
 * `keepCamera`: the same structure recolored keeps the view.
 */
export async function showMvs(
  plugin: PluginContext,
  state: MVSData_State,
  o: { keepCamera?: boolean } = {},
): Promise<Structure | null> {
  await loadMVS(plugin, state, { doNotReportErrors: false, keepCamera: o.keepCamera ?? false });
  return loadedStructure(plugin);
}

/** The first structure an MVS load made in the plugin (its model or assembly); null: none yet. */
export function loadedStructure(plugin: PluginContext): Structure | null {
  return plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data ?? null;
}

/** Mol*'s selection on `structure`, carried onto it; null: none (SeqQL's `structure` input, query.ts). */
export function selectedLoci(plugin: PluginContext, structure: Structure): StructureElement.Loci | null {
  const loci = plugin.managers.structure.selection.getLoci(structure);
  return StructureElement.Loci.is(loci) && loci.elements.length ? StructureElement.Loci.remap(loci, structure) : null;
}

/** The plugin's highlight, selection, camera and events over `structure`. */
export function pluginView3D(plugin: PluginContext, structure: Structure): StructureView3D {
  const onto = (l: unknown): StructureElement.Loci | null =>
    StructureElement.Loci.is(l) && l.structure.root === structure.root
      ? StructureElement.Loci.remap(l, structure)
      : null;
  const interactivity = plugin.managers.interactivity;
  return {
    highlight: (loci) => {
      if (loci?.elements.length) interactivity.lociHighlights.highlightOnly({ loci }, false);
      else interactivity.lociHighlights.clearHighlights();
    },
    select: (loci) => {
      if (loci?.elements.length) interactivity.lociSelects.selectOnly({ loci }, false);
      else interactivity.lociSelects.deselectAll();
    },
    focus: (loci) => plugin.managers.camera.focusLoci(loci),
    onHover: (fn) => {
      const sub = plugin.behaviors.interaction.hover.subscribe((e) => fn(onto(e.current.loci)));
      return () => sub.unsubscribe();
    },
    onSelect: (fn) => {
      const sub = plugin.managers.structure.selection.events.changed.subscribe(() =>
        fn(selectedLoci(plugin, structure)),
      );
      return () => sub.unsubscribe();
    },
  };
}

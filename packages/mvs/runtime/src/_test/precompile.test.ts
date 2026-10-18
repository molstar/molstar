/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { setFSModule } from '@molstar/core/util/data-source';
import {
  getRenderObjectSignature,
  getSignatureKey,
  mergeSignatures,
  type ObjectSignature,
} from '@molstar/graphics/gl/precompile';
import type { Representation } from '@molstar/graphics/repr/representation';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginContext } from '@molstar/plugin/context';
import { PluginSpec } from '@molstar/plugin/spec';
import { getSnapshotsPrecompileList } from '@molstar/plugin/state/precompile';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';

setFSModule(fs);

/**
 * `data/examples/mvs/kinase-story.mvsj` with every structure read from the local `1cbs_full.bcif` (selections that find
 * nothing there draw nothing) and without labels: text needs the canvas API, which Node does not have
 */
function createLocalStory() {
  const mvsj = JSON.parse(fs.readFileSync('data/examples/mvs/kinase-story.mvsj', 'utf8'));
  const visit = (node: any) => {
    if (node.kind === 'download') node.params.url = `file://${process.cwd()}/data/examples/1cbs_full.bcif`;
    node.children = node.children?.filter((c: any) => c.kind !== 'label');
    node.children?.forEach(visit);
  };
  for (const snapshot of mvsj.snapshots) visit(snapshot.root);
  return MVSData.fromMVSJ(JSON.stringify(mvsj));
}

/** The signatures of what the representations in the plugin's state drew */
function getDrawnSignatures(plugin: PluginContext) {
  const signatures: ObjectSignature[] = [];
  for (const cell of plugin.state.data.cells.values()) {
    const repr = (cell.obj?.data as { repr?: Representation<any> } | undefined)?.repr;
    for (const object of repr?.renderObjects ?? []) {
      const signature = getRenderObjectSignature(object);
      if (signature) signatures.push(signature);
    }
  }
  return signatures;
}

/** The ways `r` drew (opaque, transparent, emissive) that its derived signature `d` doesn't allow */
function getMissingFlags(r: ObjectSignature, d: ObjectSignature) {
  return (['opaque', 'transparent', 'emissive'] as const).filter((f) => r[f] && !d[f]);
}

/** The defines (and other fields) in which `r` differs from the closest of `derived`, to name what a derivation misses */
function describeMiss(r: ObjectSignature, derived: readonly ObjectSignature[]) {
  let best: string[] | undefined;
  for (const d of derived) {
    if (d.type !== r.type) continue;
    const diff = Object.keys({ ...r.defines, ...d.defines })
      .filter((k) => r.defines[k] !== d.defines[k])
      .map((k) => `${k}: ${r.defines[k]} (derived ${d.defines[k]})`);
    if (r.doubleSided !== d.doubleSided) diff.push(`doubleSided: ${r.doubleSided}`);
    if (r.merged !== d.merged) diff.push(`merged: ${r.merged}`);
    if (r.linearTextures.join() !== d.linearTextures.join()) diff.push(`linear: ${r.linearTextures}`);
    if (!best || diff.length < best.length) best = diff;
  }
  return `${r.type}${r.merged ? ' merged' : ''}: ${best ? best.join('; ') : 'no derived signature of the kind'}`;
}

describe('precompile lists of MVS stories', () => {
  it('derived from the snapshots, cover what playing them draws', async () => {
    const plugin = new PluginContext({ registry: [MVSRuntimeRegistry], behaviors: [PluginSpec.Behavior(MolViewSpec)] });
    await plugin.init();
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    // primitive labels need the canvas API
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    try {
      // `keepCamera`: the story has camera nodes and the plugin has no canvas
      await loadMVS(plugin, createLocalStory(), { keepCamera: true });
      const entries = plugin.managers.snapshot.state.entries.toArray();
      expect(entries.length).toBe(11);

      // derived before any snapshot but the first is applied, as when loading
      const { list, props } = getSnapshotsPrecompileList(
        plugin,
        entries.map((e) => e.snapshot),
      );
      // the snapshots' props, each distinct once: the story's are all the same
      expect(props.length).toBe(1);
      const derived = new Set(list.objects.map(getSignatureKey));

      const drawn = new Map<string, ObjectSignature>();
      for (const entry of entries) {
        await PluginCommands.State.Snapshots.Apply(plugin, { id: entry.snapshot.id });
        for (const s of getDrawnSignatures(plugin)) {
          const key = getSignatureKey(s);
          const known = drawn.get(key);
          drawn.set(key, known ? mergeSignatures(known, s) : s);
        }
      }
      expect(drawn.size).toBeGreaterThan(3);
      const missing = [...drawn.entries()]
        .filter(([k]) => !derived.has(k))
        .map(([, s]) => describeMiss(s, list.objects));
      expect(missing).toEqual([]);
      const byKey = new Map(list.objects.map((d) => [getSignatureKey(d), d]));
      const flags = [...drawn.values()]
        .map((r) => [r.type, getMissingFlags(r, byKey.get(getSignatureKey(r))!)] as const)
        .filter(([, f]) => f.length);
      expect(flags).toEqual([]);
      // what data decides is over-approximated (merged objects, double-sidedness by quality, the multilayer theme's
      // granularity by selections), and text derives but doesn't draw here
      expect(list.objects.length, `drawn ${drawn.size}, derived ${list.objects.length}`).toBeLessThan(drawn.size * 5);
    } finally {
      warn.mockRestore();
      error.mockRestore();
      plugin.dispose();
    }
  });
});

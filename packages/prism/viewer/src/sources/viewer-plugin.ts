/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasBackgroundProps, observableState, parseColor, type ReadonlyState } from '@molstar/prism-ui';
import type { Structure, Trajectory } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureSource, structureMvs } from '@molstar/seq-molstar-poc';
import {
  mountPlugin,
  pluginView3D,
  setBackground as setPluginBackground,
  showMvs,
} from '@molstar/seq-molstar-poc/plugin';
import type { StructureView3D } from '@molstar/seq-molstar-poc/viewer';
import type { HostStructure } from '../sequence/index.js';
import { hoverLines } from './hover-lines.js';

// The viewer's plugin source (contexts.md §7: only `sources/` and `actions/` import `@molstar/plugin`): a real Mol*
// plugin on the canvas host, composed as molstar-poc's `pluginSpec()` is (highlight, select, camera focus, loci labels,
// MolViewSpec), which the viewer's models follow through plain states: the canvas props (the ink), the 3D hover label
// (the shell's hover label), the structure shown (the sequence model's input). Nothing here knows the UI.

export interface ViewerPlugin {
  /** The parts of the canvas props that decide the background: what `canvasInk` follows. */
  readonly canvasProps: ReadonlyState<CanvasBackgroundProps>;
  /** The structure shown, with the host's parse of it; undefined until one is loaded. */
  readonly structure: ReadonlyState<HostStructure | undefined>;
  /** The 3D hover label as lines (`plugin.behaviors.labels.highlight`); empty: nothing hovered. */
  readonly hoverLines: ReadonlyState<readonly string[]>;
  /** The reason the structure could not be loaded; undefined: none. */
  readonly error: ReadonlyState<string | undefined>;
  /** The canvas background, a CSS color (`#rrggbb`). */
  setBackground(color: string): void;
  /** Shows a structure file through MolViewSpec with its default representation. */
  load(o: { url: string; format: 'bcif' | 'mmcif' | 'pdb'; title: string }): Promise<void>;
  /** The plugin's highlight, selection, camera and events over `structure`, for the link. */
  view3d(structure: Structure): StructureView3D;
  /** Resolves when the canvas is mounted. */
  readonly ready: Promise<void>;
  dispose(): void;
}

const same = (a: readonly string[], b: readonly string[]) => a.length === b.length && a.every((x, i) => x === b[i]);

/** The first trajectory of the structure hierarchy: what the structure was made from (shared with the sequence). */
function loadedTrajectory(plugin: PluginContext): Trajectory | undefined {
  return plugin.managers.structure.hierarchy.current.structures[0]?.model?.trajectory?.cell.obj?.data;
}

/** A Mol* plugin mounted on `host` (a positioned, sized element), disposed with `dispose()`. */
export function createViewerPlugin(host: HTMLElement, o: { background: string }): ViewerPlugin {
  let background = o.background;
  let plugin: PluginContext | undefined;
  let disposed = false;
  let loads = 0;
  const canvasProps = observableState<CanvasBackgroundProps>({
    renderer: { backgroundColor: parseColor(background) ?? 0xffffff },
  });
  const structure = observableState<HostStructure | undefined>(undefined);
  const lines = observableState<readonly string[]>([]);
  const error = observableState<string | undefined>(undefined);
  const subs: { unsubscribe(): void }[] = [];

  /** The canvas props as the plugin has them now (a load may set the background; so may Mol*'s own settings). */
  const sync = () => {
    const props = plugin?.canvas3d?.props;
    if (!props) return;
    canvasProps.set({
      renderer: { backgroundColor: props.renderer.backgroundColor },
      transparentBackground: props.transparentBackground,
      postprocessing: { background: props.postprocessing.background as never },
    });
  };

  const ready = mountPlugin(host, { background: parseColor(background) ?? 0xffffff }).then((p) => {
    if (disposed) {
      p.dispose();
      return;
    }
    plugin = p;
    subs.push(
      p.events.canvas3d.settingsUpdated.subscribe(sync),
      p.behaviors.labels.highlight.subscribe((e) => {
        const next = hoverLines(e.labels);
        if (!same(next, lines.value)) lines.set(next);
      }),
    );
    // A background chosen while the canvas was mounting.
    const wanted = parseColor(background);
    if (wanted !== undefined && wanted !== p.canvas3d?.props.renderer.backgroundColor) setPluginBackground(p, wanted);
    else sync();
  });

  return {
    canvasProps: canvasProps.asReadonly,
    structure: structure.asReadonly,
    hoverLines: lines.asReadonly,
    error: error.asReadonly,
    ready,
    setBackground(color) {
      background = color;
      const value = parseColor(color);
      if (plugin && value !== undefined) setPluginBackground(plugin, value);
    },
    async load({ url, format, title }) {
      await ready;
      if (!plugin || disposed) return;
      const load = ++loads;
      error.set(undefined);
      try {
        const shown = await showMvs(plugin, structureMvs({ url, format, background }));
        if (disposed || load !== loads) return;
        // MolViewSpec sets the canvas as its state says (the background we gave it); say it to the followers.
        sync();
        const trajectory = loadedTrajectory(plugin);
        if (!shown || !trajectory) throw new Error('Mol* showed no structure');
        structure.set({ key: `${load}`, title, source: StructureSource.of(trajectory, format), structure: shown });
      } catch (e) {
        if (!disposed) error.set(e instanceof Error ? e.message : String(e));
      }
    },
    view3d: (s) => {
      if (!plugin) throw new Error('the plugin is not mounted');
      return pluginView3D(plugin, s);
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      for (const s of subs) s.unsubscribe();
      plugin?.dispose();
      plugin = undefined;
    },
  };
}

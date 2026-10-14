/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasBackgroundProps, observableState, parseColor, type ReadonlyState } from '@molstar/prism-ui';
import type { Structure, StructureElement, Trajectory } from '@molstar/model/model/structure';
import { SecondaryStructureColorThemeProvider } from '@molstar/graphics/theme/color/secondary-structure';
import { PLDDTConfidenceColorThemeProvider } from '@molstar/model-archive-extension/quality-assessment/color/plddt';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureSource, structureMvs } from '@molstar/seq-molstar-poc';
import {
  mountPlugin,
  pluginView3D,
  setBackground as setPluginBackground,
  showMvs,
} from '@molstar/seq-molstar-poc/plugin';
import type { HostActions, HostStructure, SequenceView3D, StructureColoring } from '../sequence/index.js';
import { type FramingCamera, paddedFocus } from './framing.js';
import { hoverLines } from './hover-lines.js';
import { selectionSummary } from './selection.js';

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
  /**
   * The structure selection (`managers.structure.selection`), for the selection bar: its summary in Mol*'s words, empty
   * when nothing is selected; `clear` deselects (the link mirrors it into the sequence); `focus` moves the camera to it.
   */
  readonly selection: {
    readonly summary: ReadonlyState<string | undefined>;
    clear(): void;
    focus(): void;
  };
  /** The reason the structure could not be loaded; undefined: none. */
  readonly error: ReadonlyState<string | undefined>;
  /** The canvas background, a CSS color (`#rrggbb`). */
  setBackground(color: string): void;
  /**
   * The height of the band over the canvas' top, px: the camera's default framing (the first one, a reset) fits the
   * structure below it. A camera the person has moved is left where it is; one still as framed follows the change.
   */
  setTopPadding(px: number): void;
  /** What the sequence view's track menus ask of the 3D view (`HostActions`): the polymer's color theme. */
  readonly actions: HostActions;
  /** Whether the canvas is transparent (the page shows through); the ink then falls back to light. */
  setTransparent(on: boolean): void;
  /** Shows a structure file through MolViewSpec with its default representation. */
  load(o: {
    url: string;
    format: 'bcif' | 'mmcif' | 'pdb';
    title: string;
    /** A predicted model: its cartoon colored by pLDDT, and the sequence view gets a pLDDT track. */
    predicted?: boolean;
    /** Which ligands to show, by chemical component id; absent: all. */
    ligands?: readonly string[];
  }): Promise<void>;
  /**
   * The plugin's highlight, selection, camera and events over `structure`, for the link; its `resetCamera` is the default
   * framing below the band, animated (the sequence view asks for it when the person empties the selection there).
   */
  view3d(structure: Structure): SequenceView3D;
  /** Resolves when the canvas is mounted. */
  readonly ready: Promise<void>;
  dispose(): void;
}

/** The color theme name of a polymer drawn with the default coloring (what MolViewSpec's cartoon has). */
const DefaultColorTheme = 'uniform';

const same = (a: readonly string[], b: readonly string[]) => a.length === b.length && a.every((x, i) => x === b[i]);

/** The first trajectory of the structure hierarchy: what the structure was made from (shared with the sequence). */
function loadedTrajectory(plugin: PluginContext): Trajectory | undefined {
  return plugin.managers.structure.hierarchy.current.structures[0]?.model?.trajectory?.cell.obj?.data;
}

/** A Mol* plugin mounted on `host` (a positioned, sized element), disposed with `dispose()`. */
export function createViewerPlugin(host: HTMLElement, o: { background?: string } = {}): ViewerPlugin {
  // White by default: Mol*'s own warm default clashes with the UI's glass.
  let background = o.background ?? '#ffffff';
  let transparent = false;
  let topPadding = 0;
  let selectedLoci: () => StructureElement.Loci[] = () => [];
  // Whether the camera is still as the framing left it (a drag or the wheel ends that).
  let framed = true;
  let plugin: PluginContext | undefined;
  let disposed = false;
  let loads = 0;
  const canvasProps = observableState<CanvasBackgroundProps>({
    renderer: { backgroundColor: parseColor(background) ?? 0xffffff },
  });
  const structure = observableState<HostStructure | undefined>(undefined);
  const lines = observableState<readonly string[]>([]);
  const error = observableState<string | undefined>(undefined);
  const summary = observableState<string | undefined>(undefined);
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

  /** The default framing, below the band; `durationMs` 0 at once. */
  const frame = (durationMs: number) => {
    framed = true;
    plugin?.canvas3d?.requestCameraReset({
      durationMs,
      snapshot: (scene, camera) =>
        paddedFocus(
          camera as unknown as FramingCamera,
          scene.boundingSphereVisible.center,
          scene.boundingSphereVisible.radius,
          topPadding,
        ) as never,
    });
  };

  const ready = mountPlugin(host, {
    background: parseColor(background) ?? 0xffffff,
    // The themes the sequence's track menus name: the confidence (a predicted model) and the secondary structure.
    registry: [
      { structure: { themes: { color: [PLDDTConfidenceColorThemeProvider, SecondaryStructureColorThemeProvider] } } },
    ],
  }).then((p) => {
    if (disposed) {
      p.dispose();
      return;
    }
    plugin = p;
    const moved = () => {
      framed = false;
    };
    // One subscription: the selection changes at a click or the end of a drag, never per pointer move.
    const manager = p.managers.structure.selection;
    const selected = () => {
      const out: StructureElement.Loci[] = [];
      manager.entries.forEach((e) => out.push(e.selection));
      return out;
    };
    const update = () => {
      const next = selectionSummary(selected());
      if (next !== summary.value) summary.set(next);
    };
    subs.push(manager.events.changed.subscribe(update));
    selectedLoci = selected;
    const input = p.canvas3d?.input;
    if (input) subs.push(input.drag.subscribe(moved), input.wheel.subscribe(moved));
    subs.push(
      p.events.canvas3d.settingsUpdated.subscribe(sync),
      p.behaviors.labels.highlight.subscribe((e) => {
        const next = hoverLines(e.labels);
        if (!same(next, lines.value)) lines.set(next);
      }),
    );
    if (transparent) p.canvas3d?.setProps({ transparentBackground: true });
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
    actions: {
      coloring() {
        const p = plugin;
        if (!p) return undefined;
        // The polymer's representations (not the ligands' ball and stick): one color theme among them, or none of the three.
        let name: string | undefined;
        for (const st of p.managers.structure.hierarchy.current.structures) {
          for (const c of st.components) {
            for (const r of c.representations) {
              const params = r.cell.transform.params;
              if (params?.type?.name === 'ball-and-stick') continue;
              const theme: string | undefined = params?.colorTheme?.name;
              if (name !== undefined && theme !== name) return undefined;
              name = theme;
            }
          }
        }
        return name === 'plddt-confidence'
          ? 'plddt'
          : name === 'secondary-structure'
            ? 'secondary'
            : name === 'default' || name === DefaultColorTheme
              ? 'default'
              : undefined;
      },
      async colorStructure(coloring: StructureColoring) {
        const p = plugin;
        if (!p) return;
        const color =
          coloring === 'plddt' ? 'plddt-confidence' : coloring === 'secondary' ? 'secondary-structure' : 'default';
        const components = p.managers.structure.hierarchy.current.structures.flatMap((s) => s.components);
        // One undoable update of Mol*'s state, as a theme change of the structure tools is: the polymer's representations,
        // not the ligands' ball and stick.
        await p.dataTransaction(
          async () => {
            await p.managers.structure.component.updateRepresentationsTheme(components, (_c, r) =>
              r.cell.transform.params?.type?.name === 'ball-and-stick' ? {} : ({ color } as never),
            );
          },
          { canUndo: 'Color the structure' },
        );
      },
    },
    selection: {
      summary: summary.asReadonly,
      clear: () => plugin?.managers.interactivity.lociSelects.deselectAll(),
      focus: () => {
        const locis = selectedLoci();
        if (plugin && locis.length > 0) plugin.managers.camera.focusLoci(locis);
      },
    },
    ready,
    setBackground(color) {
      background = color;
      const value = parseColor(color);
      if (plugin && value !== undefined) setPluginBackground(plugin, value);
    },
    setTopPadding(px) {
      if (px === topPadding) return;
      topPadding = px;
      if (framed && structure.value) frame(250);
    },
    setTransparent(on) {
      transparent = on;
      if (!plugin) return;
      plugin.canvas3d?.setProps({ transparentBackground: on });
      plugin.events.canvas3d.settingsUpdated.next(undefined);
    },
    async load({ url, format, title, predicted, ligands }) {
      await ready;
      if (!plugin || disposed) return;
      const load = ++loads;
      error.set(undefined);
      try {
        const shown = await showMvs(
          plugin,
          structureMvs({
            url,
            format,
            background,
            ...(predicted ? { colorTheme: 'plddt-confidence' } : {}),
            ...(ligands ? { ligandComps: ligands } : {}),
          }),
        );
        if (disposed || load !== loads) return;
        // MolViewSpec sets the canvas as its state says (the background we gave it); say it to the followers.
        sync();
        const trajectory = loadedTrajectory(plugin);
        if (!shown || !trajectory) throw new Error('Mol* showed no structure');
        frame(0);
        structure.set({
          key: `${load}`,
          title,
          source: StructureSource.of(trajectory, format),
          structure: shown,
          ...(predicted ? { predicted: true } : {}),
          ...(ligands ? { ligands } : {}),
        });
      } catch (e) {
        if (!disposed) error.set(e instanceof Error ? e.message : String(e));
      }
    },
    view3d: (s) => {
      if (!plugin) throw new Error('the plugin is not mounted');
      return { ...pluginView3D(plugin, s), resetCamera: () => frame(250) };
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

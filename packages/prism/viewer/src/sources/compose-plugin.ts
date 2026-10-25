/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { decodeColor } from '@molstar/core/util/color/utils';
import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { FocusLoci as CameraFocusLoci, NoPrimaryFocusLociBindings } from '@molstar/plugin/behavior/dynamic/camera';
import { FocusLoci as RepresentationFocusLoci } from '@molstar/plugin/behavior/dynamic/representation';
import { PluginConfig, type PluginConfigItem } from '@molstar/plugin/config';
import { PluginSpec, type PluginRegistryEntry } from '@molstar/plugin/spec';
import { getGraphicsBackends } from '@molstar/plugin/util/graphics-backends';
import type { PrismViewerOptionValues } from '../options.js';
import type { PrismViewerPluginParts, ViewerExtension } from '../spec.js';

// The plugin spec a viewer runs on (viewer-spec.md §2.1): four layers in a fixed order, nothing
// resolved by name. Registry: the base, then each extension's entries in list order (custom formats are an extension's
// entry, so they come after the built-in formats). Behaviors: the base (narrowed by the focus option), then the
// extensions'. Config: the base, the extensions', then the options' and the viewer's fixed items; later wins. Then
// `finalize`. Registry conflicts throw at `init()` (the plugin's `register`, atomic, naming every one); a behavior twice
// throws here, naming its owners.

/** Config the viewer fixes: undo is the viewer's (sequence.md §7). */
const FixedConfig: readonly (readonly [PluginConfigItem, unknown])[] = [[PluginConfig.State.HistoryCapacity, 0]];

export interface ComposeInput {
  readonly base: PrismViewerPluginParts;
  /** The eager extensions, in order; a lazy one's plugin part registers when it loads. */
  readonly extensions: readonly ViewerExtension[];
  readonly options: PrismViewerOptionValues;
  readonly finalize?: (spec: PluginSpec) => PluginSpec;
}

export class PluginCompositionError extends Error {
  constructor(readonly conflicts: readonly string[]) {
    super(`The viewer's plugin parts conflict: ${conflicts.join('; ')}`);
    this.name = 'PluginCompositionError';
  }
}

export function composePluginSpec(input: ComposeInput): PluginSpec {
  const { base, extensions, options } = input;
  const registry: PluginRegistryEntry[] = [...base.registry];
  const behaviors: { behavior: PluginSpec.Behavior; owner: string }[] = focusBehaviors(
    base.behaviors,
    options.canvas.focusBehavior,
  ).map((behavior) => ({ behavior, owner: 'core' }));
  const config: [PluginConfigItem, unknown][] = (base.config ?? []).map(([k, v]) => [k, v]);
  let canvas3d: PartialCanvas3DProps = { ...base.canvas3d };

  for (const e of extensions) {
    if (!e.plugin) continue;
    if (e.plugin.registry) registry.push(...e.plugin.registry);
    for (const b of e.plugin.behaviors ?? []) behaviors.push({ behavior: b, owner: e.id });
    for (const [k, v] of e.plugin.config ?? []) config.push([k, v]);
    if (e.plugin.canvas3d) canvas3d = mergeCanvas(canvas3d, e.plugin.canvas3d);
  }

  const conflicts: string[] = [];
  const seen = new Map<string, string>();
  for (const { behavior, owner } of behaviors) {
    const id = behavior.transformer.id;
    const first = seen.get(id);
    if (first !== undefined) conflicts.push(`behavior '${id}' is brought by ${first} and by ${owner}`);
    else seen.set(id, owner);
  }
  if (conflicts.length) throw new PluginCompositionError(conflicts);

  config.push(...optionsConfig(options), ...FixedConfig.map(([k, v]) => [k, v] as [PluginConfigItem, unknown]));
  canvas3d = mergeCanvas(canvas3d, optionsCanvas(options));
  // Reduced motion (accessibility.md): the plugin's own camera moves are cuts too (framing on load, focus by click).
  const reduced = prefersReducedMotion();
  if (reduced) canvas3d = mergeCanvas(canvas3d, { cameraResetDurationMs: 0 });

  const spec: PluginSpec = {
    registry,
    behaviors: behaviors.map(({ behavior: b }) =>
      reduced && b.transformer === CameraFocusLoci
        ? PluginSpec.Behavior(CameraFocusLoci, { ...b.defaultParams, durationMs: 0 })
        : b,
    ),
    config: lastWins(config),
    canvas3d,
  };
  return input.finalize ? input.finalize(spec) : spec;
}

/** The camera's focus behaviors as the option says (v5's `viewportFocusBehavior`). */
function focusBehaviors(
  base: readonly PluginSpec.Behavior[],
  focus: PrismViewerOptionValues['canvas']['focusBehavior'],
): PluginSpec.Behavior[] {
  if (focus === 'default') return [...base];
  const isFocus = (b: PluginSpec.Behavior) =>
    b.transformer === CameraFocusLoci || b.transformer === RepresentationFocusLoci;
  const rest = base.filter((b) => !isFocus(b));
  if (focus === 'disabled') return rest;
  return [...rest, PluginSpec.Behavior(CameraFocusLoci, { bindings: NoPrimaryFocusLociBindings })];
}

/** The renderer and data options as config items; `auto` leaves the plugin's own detection. */
export function optionsConfig(o: PrismViewerOptionValues): [PluginConfigItem, unknown][] {
  const r = o.renderer;
  const d = o.data;
  const items: [PluginConfigItem, unknown][] = [
    [PluginConfig.General.GraphicsBackends, getGraphicsBackends(r.graphicsBackend)],
    [PluginConfig.General.PixelScale, r.pixelScale],
    [PluginConfig.General.PickScale, r.pickScale],
    [PluginConfig.General.ResolutionMode, r.resolutionMode],
    [PluginConfig.General.PowerPreference, r.powerPreference],
    [PluginConfig.General.AllowMajorPerformanceCaveat, r.allowMajorPerformanceCaveat],
    [PluginConfig.General.DisableAntialiasing, r.disableAntialiasing],
    [PluginConfig.Download.DefaultPdbProvider, d.pdbProvider],
    [PluginConfig.Download.DefaultEmdbProvider, d.emdbProvider],
    [PluginConfig.VolumeStreaming.DefaultServer, d.volumeStreamingServer],
    [PluginConfig.VolumeStreaming.Enabled, !d.volumeStreamingDisabled],
    [PluginConfig.State.DefaultServer, d.stateServer],
    [PluginConfig.State.CurrentServer, d.stateServer],
    [PluginConfig.Structure.SaccharideCompIdMapType, d.saccharideCompIdMapType],
    [PluginConfig.General.PauseLoopWhenIdle, o.ui.pauseWhenIdle],
  ];
  if (r.transparency !== 'auto') items.push([PluginConfig.General.Transparency, r.transparency]);
  return items;
}

/** The canvas options as canvas props. */
export function optionsCanvas(o: PrismViewerOptionValues): PartialCanvas3DProps {
  const props: PartialCanvas3DProps = {};
  const background = o.canvas.background ? decodeColor(o.canvas.background) : undefined;
  if (background !== undefined) props.renderer = { backgroundColor: background };
  if (o.renderer.illumination) props.illumination = { enabled: true };
  return props;
}

function mergeCanvas(a: PartialCanvas3DProps, b: PartialCanvas3DProps): PartialCanvas3DProps {
  const out: Record<string, unknown> = { ...a };
  for (const [k, v] of Object.entries(b)) {
    const prev = out[k];
    out[k] =
      prev && typeof prev === 'object' && !Array.isArray(prev) && v && typeof v === 'object' && !Array.isArray(v)
        ? { ...(prev as object), ...(v as object) }
        : v;
  }
  return out as PartialCanvas3DProps;
}

/** One value per config item, the last set; in the order each item was first set. */
function lastWins(items: readonly [PluginConfigItem, unknown][]): [PluginConfigItem, unknown][] {
  const map = new Map<PluginConfigItem, unknown>();
  for (const [k, v] of items) map.set(k, v);
  return [...map];
}

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

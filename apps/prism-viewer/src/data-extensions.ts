/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Backgrounds } from '@molstar/backgrounds-extension';
import { DnatcoNtCs } from '@molstar/dnatco-extension/behavior';
import { G3DFormat, G3dProvider } from '@molstar/g3d-extension/format';
import { PDBeStructureQualityReport } from '@molstar/pdbe-extension/structure-quality-report/behavior';
import { StructureQualityReportColorThemeProvider } from '@molstar/pdbe-extension/structure-quality-report/color';
import { PluginSpec } from '@molstar/plugin/spec';
import type { AnnotationProvider } from '@molstar/prism-viewer/api/annotations';
import type { ViewerExtension } from '@molstar/prism-viewer/spec';
import { RCSBValidationReport } from '@molstar/rcsb-extension/validation-report/behavior';
import { DensityFitColorThemeProvider } from '@molstar/rcsb-extension/validation-report/color/density-fit';
import { GeometryQualityColorThemeProvider } from '@molstar/rcsb-extension/validation-report/color/geometry-quality';
import { RandomCoilIndexColorThemeProvider } from '@molstar/rcsb-extension/validation-report/color/random-coil-index';
import { SbNcbrPartialChargesColorThemeProvider } from '@molstar/sb-ncbr-extension/partial-charges/color';
import { SbNcbrPartialCharges } from '@molstar/sb-ncbr-extension/partial-charges/behavior';
import { SbNcbrTunnels } from '@molstar/sb-ncbr-extension/tunnels/behavior';
import { wwPDBChemicalComponentDictionary } from '@molstar/wwpdb-extension/ccd/behavior';

// Mol*'s data-only extensions in the Prism viewer (plan 5.2, parity.md §6, README decision 27): what each registered
// for the v5 plugin (its behavior: custom properties, color themes, representations, presets, formats), nothing of
// their own in the API or the UI. They are reached through the generic commands and tools (`representation.color`
// with the theme's name, `load.files` for a format, `load.action` for a plugin action); what a structure carries of
// them is shown by `annotations.apply`. None of these modules imports React, so they load headless too.

/**
 * An annotation from a color theme: it applies where the theme does (the theme provider's own `isApplicable`, which
 * knows the experiment and the entry), and colors the target with the theme. The themes of the validation reports and
 * the partial charges fetch or read what they need when they are applied (`ensureCustomProperties`).
 */
function themeAnnotation(
  theme: {
    readonly name: string;
    isApplicable?(ctx: { structure?: Parameters<AnnotationProvider['isApplicable']>[0] }): boolean;
  },
  a: Pick<AnnotationProvider, 'id' | 'title' | 'description'> & { target: 'all' | 'polymer' },
): AnnotationProvider {
  return {
    id: a.id,
    title: a.title,
    ...(a.description ? { description: a.description } : {}),
    isApplicable: (s) => theme.isApplicable?.({ structure: s }) ?? true,
    apply: (ctx, structure) => ctx.api.representation.color({ target: a.target, theme: theme.name, structure }),
  };
}

const Behaviors = (...behaviors: Parameters<typeof PluginSpec.Behavior>[0][]) =>
  behaviors.map((b) => PluginSpec.Behavior(b));

export const DataOnlyExtensions: readonly ViewerExtension[] = [
  {
    id: 'dnatco',
    title: 'DNATCO',
    // NtC annotations of nucleic acids: the confal pyramids and the NtC tube (representations, color themes, presets).
    plugin: { registry: [], behaviors: Behaviors(DnatcoNtCs) },
  },
  {
    id: 'rcsb',
    title: 'RCSB validation',
    // The wwPDB validation report through RCSB: three color themes, the clashes representation, presets, labels.
    plugin: { registry: [], behaviors: Behaviors(RCSBValidationReport) },
    annotations: [
      themeAnnotation(GeometryQualityColorThemeProvider, {
        id: 'rcsb-geometry-quality',
        title: 'Geometry quality (RCSB)',
        description: "The wwPDB validation report's geometry outliers per residue, as color",
        target: 'all',
      }),
      themeAnnotation(DensityFitColorThemeProvider, {
        id: 'rcsb-density-fit',
        title: 'Density fit (RCSB)',
        description: "The wwPDB validation report's fit of the model to the density, as color",
        target: 'all',
      }),
      themeAnnotation(RandomCoilIndexColorThemeProvider, {
        id: 'rcsb-random-coil-index',
        title: 'Random coil index (RCSB)',
        description: "The wwPDB validation report's random coil index of an NMR entry, as color",
        target: 'polymer',
      }),
    ],
  },
  {
    id: 'pdbe',
    title: 'PDBe structure quality',
    plugin: { registry: [], behaviors: Behaviors(PDBeStructureQualityReport) },
    annotations: [
      themeAnnotation(StructureQualityReportColorThemeProvider, {
        id: 'pdbe-structure-quality',
        title: 'Structure quality (PDBe)',
        description: "PDBe's validation issues per residue, as color",
        target: 'all',
      }),
    ],
  },
  {
    id: 'wwpdb',
    title: 'wwPDB chemical components',
    // Presets for data loaded from the chemical component dictionary. (Its struct_conn functions are plugin helpers
    // with no behavior: the v5 `PluginExtensions.wwPDBStructConn`, not part of the viewer's operations.)
    plugin: { registry: [], behaviors: Behaviors(wwPDBChemicalComponentDictionary) },
  },
  {
    id: 'sb-ncbr',
    title: 'SB NCBR',
    // Partial charges (a model property, a color theme, a preset, labels) and the tunnels actions and preset.
    plugin: { registry: [], behaviors: Behaviors(SbNcbrPartialCharges, SbNcbrTunnels) },
    annotations: [
      themeAnnotation(SbNcbrPartialChargesColorThemeProvider, {
        id: 'partial-charges',
        title: 'Partial charges',
        description: "A structure's partial charges (SB NCBR), as color",
        target: 'all',
      }),
    ],
  },
  {
    id: 'g3d',
    title: 'G3D',
    // The G3D trajectory format (`.g3d`: `load.files` reads it), its labels and its load action.
    plugin: { registry: [{ formats: [G3dProvider] }], behaviors: Behaviors(G3DFormat) },
  },
  {
    id: 'backgrounds',
    title: 'Backgrounds',
    // The background styles the canvas offers (a config item), with their images and skybox.
    plugin: { registry: [], behaviors: Behaviors(Backgrounds) },
  },
];

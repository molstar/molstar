/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PostprocessingParams } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { PluginConfig } from '@molstar/plugin/config';
import { type ApiContext, ApiArea, ApiResult } from './area.js';

// `styles` (api.md §2): Mol*'s Quick Styles as one operation. A representation style applies a preset to every
// structure (through `structure.preset`, part of this step); a look sets how everything is lit and outlined (the
// component options' `ignoreLight`, the canvas's outline and occlusion), as the classic panel does.

export const RepresentationStyles = {
  default: { label: 'Default', preset: '' },
  cartoon: { label: 'Cartoon', preset: 'preset-structure-representation-polymer-and-ligand' },
  spacefill: { label: 'Spacefill', preset: 'preset-structure-representation-illustrative' },
  surface: { label: 'Surface', preset: 'preset-structure-representation-molecular-surface' },
} as const;

export const Looks = {
  plain: { label: 'Plain' },
  illustrative: { label: 'Illustrative' },
} as const;

export type StyleId = keyof typeof RepresentationStyles | keyof typeof Looks;

const StyleOptions = [
  ...Object.entries(RepresentationStyles).map(([id, s]) => [id, s.label] as [StyleId, string]),
  ...Object.entries(Looks).map(([id, s]) => [id, s.label] as [StyleId, string]),
];

export class StylesApi extends ApiArea {
  readonly id = 'styles';

  readonly apply = this.operation({
    id: 'styles.apply',
    title: 'Apply a style',
    params: { style: PD.Select<StyleId>('default', StyleOptions, { description: 'A quick style' }) },
    async run({ style }, ctx) {
      if (style in Looks) return applyLook(ctx, style as keyof typeof Looks);
      const preset = presetOf(ctx, style as keyof typeof RepresentationStyles);
      if (!preset) {
        return ApiResult.fail('unavailable', `The ${style} style needs a preset this viewer does not have.`);
      }
      const r = await ctx.api.structure.preset({ preset });
      if (!r.ok) return r;
      // A preset may bring its own lighting (spacefill ignores light): the look in effect is kept.
      await keepLook(ctx);
      return ApiResult.ok(
        `Applied the ${RepresentationStyles[style as keyof typeof RepresentationStyles].label} style.`,
        r.changed,
      );
    },
  });
}

/** The preset the viewer is configured to apply by default (`PluginConfig.Structure.DefaultRepresentationPreset`). */
export function defaultPreset(ctx: ApiContext): string | undefined {
  return presetOf(ctx, 'default');
}

function presetOf(ctx: ApiContext, style: keyof typeof RepresentationStyles): string | undefined {
  const builder = ctx.plugin.builders.structure.representation;
  const id =
    style === 'default'
      ? (ctx.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) ??
        'preset-structure-representation-auto')
      : RepresentationStyles[style].preset;
  return builder.resolveProvider(id) ? id : style === 'default' && builder.resolveProvider('auto') ? 'auto' : undefined;
}

/** The look in effect: illustrative when everything ignores light. */
function currentLook(ctx: ApiContext): keyof typeof Looks {
  return ctx.plugin.managers.structure.component.state.options.ignoreLight ? 'illustrative' : 'plain';
}

async function keepLook(ctx: ApiContext) {
  const look = currentLook(ctx);
  const options = ctx.plugin.managers.structure.component.state.options;
  // the look is the scene's: forced onto every representation, a preset's own lighting included (api.md §2 "Look")
  await ctx.plugin.managers.structure.component.setOptions(
    { ...options, ignoreLight: look === 'illustrative' },
    { reset: ['ignoreLight'] },
  );
}

async function applyLook(ctx: ApiContext, look: keyof typeof Looks) {
  const options = ctx.plugin.managers.structure.component.state.options;
  const illustrative = look === 'illustrative';
  const canvas = ctx.plugin.canvas3d;
  const pp = canvas?.props.postprocessing;
  const already = options.ignoreLight === illustrative && (!pp || (pp.outline.name === 'on') === illustrative);
  if (already) return ApiResult.already(`The look is already ${Looks[look].label.toLowerCase()}.`);
  // a look is the scene's: every representation takes it, those with their own light included
  await ctx.plugin.managers.structure.component.setOptions(
    { ...options, ignoreLight: illustrative },
    { reset: ['ignoreLight'] },
  );
  if (pp) {
    if (!illustrative) {
      const p = PD.getDefaultValues(PostprocessingParams);
      ctx.plugin.setCanvas3DProps({ postprocessing: { outline: p.outline, occlusion: p.occlusion, shadow: p.shadow } });
    } else {
      ctx.plugin.setCanvas3DProps({
        postprocessing: {
          outline: {
            name: 'on',
            params:
              pp.outline.name === 'on'
                ? pp.outline.params
                : { scale: 1, color: Color(0x000000), threshold: 0.33, includeTransparent: true },
          },
          occlusion: {
            name: 'on',
            params:
              pp.occlusion.name === 'on'
                ? pp.occlusion.params
                : {
                    multiScale: { name: 'off', params: {} },
                    radius: 5,
                    bias: 0.8,
                    blurKernelSize: 15,
                    blurDepthBias: 0.5,
                    samples: 32,
                    resolutionScale: 1,
                    color: Color(0x000000),
                    transparentThreshold: 0.4,
                  },
          },
          shadow: { name: 'off', params: {} },
        },
      });
    }
  }
  return ApiResult.ok(`The look is ${Looks[look].label.toLowerCase()}.`, { view: true });
}

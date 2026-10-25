/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { ObjectSignature } from '@molstar/graphics/gl/precompile';
import type { SignatureContext } from './signatures.js';
import type { Visual } from './visual.js';

/**
 * Creates the visuals of a representation, given the data their `mustRecreate` and `createOrUpdate` get (a symmetry
 * group, a structure, a volume and key, particles). One that creates the same visual whatever the data, props and
 * context takes the material id only; one that chooses is made by `selectVisual` or `chooseVisual` and names the
 * visuals it may create, as signatures are taken without data (`getFactorySignatures`)
 */
export type VisualFactory<D, P extends PD.Params> =
  | ((materialId: number) => Visual<D, P>)
  | SelectingVisualFactory<D, P>;

export interface SelectingVisualFactory<D, P extends PD.Params> {
  (materialId: number, data: D, props: PD.Values<P>, graphics?: GraphicsContext): Visual<D, P>;
  /** the visuals `props` and `graphics` allow, before any data */
  readonly signatureVisuals: (props: PD.Values<P>, graphics?: GraphicsContext) => Visual<D, P>[];
}

/** The visual `factory` creates for `data` */
export function createVisual<D, P extends PD.Params>(
  factory: VisualFactory<D, P>,
  materialId: number,
  data: D,
  props: PD.Values<P>,
  graphics?: GraphicsContext,
) {
  return 'signatureVisuals' in factory ? factory(materialId, data, props, graphics) : factory(materialId);
}

/**
 * A factory creating the variant `select` chooses by the data, props and context (e.g. GPU or CPU surfaces by a
 * unit's size); a visual is recreated once `select` chooses another variant (`mustRecreate`). `candidates`, the variants
 * the props and context allow before any data, name the factory's visuals for signatures
 */
export function selectVisual<D, P extends PD.Params, K extends string>(
  variants: { readonly [k in K]: (materialId: number) => Visual<D, P> },
  select: (data: D, props: PD.Values<P>, graphics?: GraphicsContext) => K,
  candidates: (props: PD.Values<P>, graphics?: GraphicsContext) => readonly K[],
): SelectingVisualFactory<D, P> {
  const factory = (materialId: number, data: D, props: PD.Values<P>, graphics?: GraphicsContext) => {
    const kind = select(data, props, graphics);
    const visual = variants[kind](materialId);
    const mustRecreate = visual.mustRecreate;
    visual.mustRecreate = (data, props, graphics) =>
      select(data, props, graphics) !== kind || !!mustRecreate?.(data, props, graphics);
    return visual;
  };
  return Object.assign(factory, {
    signatureVisuals: (props: PD.Values<P>, graphics?: GraphicsContext) =>
      candidates(props, graphics).map((k) => variants[k](-1)),
  });
}

/** As `selectVisual`, for factories choosing by the props and context only (e.g. impostors where supported) */
export function chooseVisual<D, P extends PD.Params, K extends string>(
  variants: { readonly [k in K]: (materialId: number) => Visual<D, P> },
  choose: (props: PD.Values<P>, graphics?: GraphicsContext) => K,
): SelectingVisualFactory<D, P> {
  return selectVisual(
    variants,
    (_data, props, graphics) => choose(props, graphics),
    (props, graphics) => [choose(props, graphics)],
  );
}

/** The signatures of the visuals `factory` may create with `props` in `ctx` */
export function getFactorySignatures<D, P extends PD.Params>(
  factory: VisualFactory<D, P>,
  props: PD.Values<P>,
  ctx: SignatureContext,
): ObjectSignature[] {
  const visuals = 'signatureVisuals' in factory ? factory.signatureVisuals(props, ctx.graphics) : [factory(-1)];
  return visuals.flatMap((v) => v.getSignatures?.(props, ctx) ?? []);
}

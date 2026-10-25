/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { ObjectSignature } from '@molstar/graphics/gl/precompile';
import type { SignatureContext } from '../signatures.js';
import type { Visual } from '../visual.js';
import {
  chooseVisual,
  createVisual,
  getFactorySignatures,
  selectVisual,
  type VisualFactory,
} from '../visual-factory.js';

type Data = { size: number };
type Params = PD.Params;

/** a visual of `kind`, recreated by itself while `recreate.value` is set */
function fakeVisual(kind: string, recreate = { value: false }) {
  return {
    kind,
    mustRecreate: () => recreate.value,
    getSignatures: () => [{ kind } as unknown as ObjectSignature],
  } as unknown as Visual<Data, Params> & { kind: string };
}

const kindOf = (visual: Visual<Data, Params>) => (visual as Visual<Data, Params> & { kind?: string }).kind;

/** the kinds of the visuals `factory` names for signatures with `props` */
function signatureKinds(factory: VisualFactory<Data, Params>, props: PD.Values<Params>, graphics?: GraphicsContext) {
  const ctx = { graphics } as SignatureContext;
  return getFactorySignatures(factory, props, ctx).map((s) => (s as unknown as { kind: string }).kind);
}

describe('selectVisual', () => {
  const recreate = { value: false };
  const factory = selectVisual(
    { large: () => fakeVisual('large', recreate), small: () => fakeVisual('small') },
    (data: Data) => (data.size > 10 ? 'large' : 'small'),
    (props: PD.Values<Params>) => (props.onlySmall ? ['small'] : ['large', 'small']),
  );

  it('creates the variant chosen by the data', () => {
    expect(kindOf(createVisual(factory, 0, { size: 20 }, {}))).toBe('large');
    expect(kindOf(createVisual(factory, 0, { size: 5 }, {}))).toBe('small');
  });

  it('recreates a visual once the data chooses another variant', () => {
    const visual = createVisual(factory, 0, { size: 20 }, {});
    expect(visual.mustRecreate?.({ size: 30 }, {})).toBe(false);
    expect(visual.mustRecreate?.({ size: 5 }, {})).toBe(true);
  });

  it('keeps the variant reasons to recreate', () => {
    const visual = createVisual(factory, 0, { size: 20 }, {});
    recreate.value = true;
    expect(visual.mustRecreate?.({ size: 30 }, {})).toBe(true);
    recreate.value = false;
  });

  it('names the candidate variants for signatures', () => {
    expect(signatureKinds(factory, {})).toEqual(['large', 'small']);
    expect(signatureKinds(factory, { onlySmall: true })).toEqual(['small']);
  });
});

describe('chooseVisual', () => {
  const factory = chooseVisual(
    { impostor: () => fakeVisual('impostor'), mesh: () => fakeVisual('mesh') },
    (props: PD.Values<Params>, graphics?: GraphicsContext) => (props.tryUseImpostor && graphics ? 'impostor' : 'mesh'),
  );
  const graphics = {} as GraphicsContext;

  it('creates the variant chosen by the props and context, whatever the data', () => {
    expect(kindOf(createVisual(factory, 0, { size: 1 }, { tryUseImpostor: true }, graphics))).toBe('impostor');
    expect(kindOf(createVisual(factory, 0, { size: 1 }, { tryUseImpostor: true }))).toBe('mesh');
  });

  it('recreates a visual once the props choose another variant, not otherwise', () => {
    const visual = createVisual(factory, 0, { size: 1 }, { tryUseImpostor: false }, graphics);
    expect(visual.mustRecreate?.({ size: 1 }, { tryUseImpostor: false }, graphics)).toBe(false);
    expect(visual.mustRecreate?.({ size: 1 }, { tryUseImpostor: true }, graphics)).toBe(true);
  });

  it('names only the chosen variant for signatures', () => {
    expect(signatureKinds(factory, { tryUseImpostor: true }, graphics)).toEqual(['impostor']);
    expect(signatureKinds(factory, { tryUseImpostor: true })).toEqual(['mesh']);
  });
});

describe('plain visual factories', () => {
  it('are called with the material id only, for visuals and for signatures', () => {
    const calls: number[] = [];
    const factory = (...args: unknown[]) => {
      calls.push(args.length);
      return fakeVisual('plain');
    };
    expect(kindOf(createVisual(factory, 7, { size: 1 }, {}))).toBe('plain');
    expect(signatureKinds(factory, {})).toEqual(['plain']);
    expect(calls).toEqual([1, 1]);
  });
});

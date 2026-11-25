/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DirectVolume } from '../direct-volume.js';

/** a transfer function of alpha 1 between `lo` and `hi`, 0 elsewhere */
function plateau(lo: number, hi: number) {
  return [Vec2.create(lo, 0), Vec2.create(lo + 0.01, 1), Vec2.create(hi - 0.01, 1), Vec2.create(hi, 0)];
}

function values(props: Partial<PD.Values<DirectVolume.Params>>) {
  return DirectVolume.Utils.createValuesSimple(DirectVolume.createEmpty(), props, ColorNames.orange, 1);
}

describe('direct-volume render mode', () => {
  it('is composite by default', () => {
    expect(values({}).dRenderMode.ref.value).toBe('composite');
  });

  it('starts the maximum where the transfer function rises above zero', () => {
    const v = values({ renderMode: 'maximum', controlPoints: plateau(0.2, 0.6) });
    expect(v.dRenderMode.ref.value).toBe('maximum');
    expect(v.uMaximumStart.ref.value).toBeCloseTo(0.2, 1);
  });

  it('starts the maximum above every value where the transfer function is zero throughout', () => {
    expect(values({ renderMode: 'maximum', controlPoints: [] }).uMaximumStart.ref.value).toBeGreaterThan(1);
  });

  it('skips the bricks below the window in the maximum mode, those the transfer function hides in the others', () => {
    const dv = DirectVolume.createEmpty();
    // two bricks: values 0-0.1 and 0.7-0.9 (times 255), the transfer function 0.2-0.6
    ValueCell.update(dv.brickRange, { array: new Uint8Array([0, 26, 178, 230]), dim: Vec3.create(2, 1, 1) });
    const visible = (renderMode: PD.Values<DirectVolume.Params>['renderMode']) =>
      Array.from(
        DirectVolume.Utils.createValuesSimple(
          dv,
          { renderMode, controlPoints: plateau(0.2, 0.6) },
          ColorNames.orange,
          1,
        ).tBrickTex.ref.value.array.slice(0, 2),
      );
    expect(visible('composite')).toEqual([0, 0]);
    expect(visible('additive')).toEqual([0, 0]);
    expect(visible('attenuated')).toEqual([0, 0]);
    expect(visible('maximum')).toEqual([0, 255]);
  });

  it('follows the transfer function and the mode on update', () => {
    const v = values({ renderMode: 'maximum', controlPoints: plateau(0.2, 0.6) });
    const props = { ...PD.getDefaultValues(DirectVolume.Params), renderMode: 'composite' as const };
    DirectVolume.Utils.updateValues(v, { ...props, controlPoints: plateau(0.5, 0.9) });
    expect(v.dRenderMode.ref.value).toBe('composite');
    expect(v.uMaximumStart.ref.value).toBeCloseTo(0.5, 1);
  });

  it('follows the attenuation on update', () => {
    const v = values({ renderMode: 'attenuated' });
    expect(v.uAttenuation.ref.value).toBe(0.15);
    DirectVolume.Utils.updateValues(v, {
      ...PD.getDefaultValues(DirectVolume.Params),
      renderMode: 'attenuated',
      attenuation: 1,
    });
    expect(v.uAttenuation.ref.value).toBe(1);
  });
});

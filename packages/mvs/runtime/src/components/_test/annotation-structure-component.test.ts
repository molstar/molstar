/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Structure } from '@molstar/model/model/structure';
import { createModel } from '../../helpers/_test/annotation-rows-fixture.js';
import { MVSAnnotationStructureComponent } from '../annotation-structure-component.js';

const dispose = (b: any, cache: unknown) => (MVSAnnotationStructureComponent as any).definition.dispose({ b, cache });

describe('MVSAnnotationStructureComponent dispose', async () => {
  const structure = Structure.ofModel(await createModel({ nChains: 2, nResidues: 5 }));
  const part = structure.subsetBuilder(true);
  part.setUnit(structure.units[0].id, structure.units[0].elements as any);
  const child = part.getStructure();

  it('disposes the properties of a component which is its own structure', () => {
    const spy = vi.spyOn(child.customPropertyDescriptors, 'dispose');
    dispose({ data: child }, { isSource: false });
    expect(spy).toHaveBeenCalledTimes(1);
  });
  it('does not dispose those of the source (everything is selected) or of the shared empty structure', () => {
    const source = vi.spyOn(structure.customPropertyDescriptors, 'dispose');
    dispose({ data: structure }, { isSource: true });
    expect(source).not.toHaveBeenCalled();
    const empty = vi.spyOn(Structure.Empty.customPropertyDescriptors, 'dispose');
    dispose({ data: Structure.Empty }, { isSource: false });
    dispose({ data: Structure.Empty }, undefined);
    expect(empty).not.toHaveBeenCalled();
  });
  it('has nothing to dispose without a component', () => {
    expect(() => dispose(undefined, undefined)).not.toThrow();
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { EveryLoci } from '@molstar/model/model/loci';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { createTemplateCopyLayers, getCopyLayerStates } from '../target/copy-layers.js';
import type { TemplateCopies } from '../target/expanded.js';
import { CopyLayerBlock, CopyLayerKinds, CopyLayerValue } from '@molstar/graphics/gl/instance-expansion';

function createList(fields: Partial<ParticleList> & Pick<ParticleList, 'count' | 'coordinates'>): ParticleList {
  return {
    targets: new Int32Array(fields.count),
    targetInfo: new Map([[0, {}]]),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    ...fields,
  };
}

/** three copies (top particles 0 to 2) of a template of four particles */
function createNested() {
  const template = createList({ count: 4, coordinates: new Float32Array(12) });
  const list = createList({
    count: 3,
    coordinates: new Float32Array(9),
    targetMapping: new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: template, center: Vec3.create(0, 0, 0) }],
    ]),
  });
  return { list, template, tops: [0, 1, 2] };
}

const layer = (loci: any, value: number) => ({ loci, value });
const states = (n: ReturnType<typeof createNested>, layers: ReturnType<typeof layer>[]) =>
  getCopyLayerStates(n.list, n.template, n.tops, 'transparency', { layers });
/** as `applyTransparencyValue` writes it */
const byte = (v: number) => Math.floor(v * 255);

describe('copy layers', () => {
  it('gives a copy reached whole a value, one in part a block, the others none', () => {
    const n = createNested();
    const { list } = n;
    const s = states(n, [
      layer(ParticleNesting.ofCopy(list, 1), 0.5),
      layer(ParticleNesting.ofInstance(list, 2, 3), 1),
    ]);
    expect([...s.keys()].sort()).toEqual([1, 2]);
    expect(Array.from(s.get(1)!.value!)).toEqual([byte(0.5)]);
    expect(Array.from(s.get(2)!.block!)).toEqual([0, 0, 0, byte(1)]);
  });

  it('leaves a template particle in every copy to the template, a copy with its own taking it in order', () => {
    const n = createNested();
    const { list, template } = n;
    const everyCopy = ParticleNesting.ofTemplateParticles(list, template, OrderedSet.ofSingleton(1));
    // a copy whole, then a template particle of all: the copy's value a block with it
    const s = states(n, [layer(ParticleNesting.ofCopy(list, 0), 0.2), layer(everyCopy, 1)]);
    expect([...s.keys()]).toEqual([0]);
    expect(Array.from(s.get(0)!.block!)).toEqual([byte(0.2), byte(1), byte(0.2), byte(0.2)]);
    // the template particle first, then one instance of a copy: its block from the template's
    const t = states(n, [layer(everyCopy, 1), layer(ParticleNesting.ofInstance(list, 2, 0), 0.4)]);
    expect([...t.keys()]).toEqual([2]);
    expect(Array.from(t.get(2)!.block!)).toEqual([byte(0.4), byte(1), 0, 0]);
  });

  it('gives none for every loci and all copies whole, a value for a copy whole after them', () => {
    const n = createNested();
    const { list } = n;
    expect(states(n, [layer(EveryLoci, 1)]).size).toBe(0);
    expect(states(n, [layer(Particle.Loci(list, OrderedSet.ofBounds(0, 3)), 1)]).size).toBe(0);
    const s = states(n, [layer(EveryLoci, 1), layer(ParticleNesting.ofCopy(list, 1), 0)]);
    expect(Array.from(s.get(1)!.value!)).toEqual([0]);
  });

  it('layers a copy in part whole past the room: no block, its value, the layers after it whole', () => {
    const n = createNested();
    const { list } = n;
    // a block four bytes (four particles' transparency): room for one
    const room = { bytes: 4, coarse: 0 };
    const s = getCopyLayerStates(
      list,
      n.template,
      n.tops,
      'transparency',
      {
        layers: [
          layer(ParticleNesting.ofInstance(list, 0, 3), 1),
          layer(ParticleNesting.ofInstance(list, 1, 2), 0.4),
          layer(ParticleNesting.ofInstance(list, 1, 0), 0.2),
        ],
      },
      room,
    );
    expect(Array.from(s.get(0)!.block!)).toEqual([0, 0, 0, byte(1)]);
    expect(s.get(1)!.block).toBeUndefined();
    expect(s.get(1)!.coarse).toBe(true);
    // its layers whole in order: the last's
    expect(Array.from(s.get(1)!.value!)).toEqual([byte(0.2)]);
    expect(room).toEqual({ bytes: 0, coarse: 1 });
    // a block's own cost given (the paged store's span)
    const r = { bytes: 10, coarse: 0 };
    const t = getCopyLayerStates(
      list,
      n.template,
      n.tops,
      'transparency',
      { layers: [layer(ParticleNesting.ofInstance(list, 0, 3), 1), layer(ParticleNesting.ofInstance(list, 1, 2), 1)] },
      r,
      8,
    );
    expect(t.get(0)!.block).toBeDefined();
    expect(t.get(1)!.coarse).toBe(true);
    expect(r).toEqual({ bytes: 2, coarse: 1 });
  });

  it("keeps a layout's layers within its budget across kinds, the kinds kept holding theirs, warned once", () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const n = createNested();
    const { list, template } = n;
    const copies = { tops: Int32Array.from(n.tops) } as unknown as TemplateCopies;
    // room for two transparency blocks (four bytes each), one overpaint (sixteen)
    const layers = createTemplateCopyLayers(() => 8);
    layers.layout(list, [[template, copies]]);
    const part = (top: number) => ParticleNesting.ofInstance(list, top, 1);
    layers.set(list, { transparency: { layers: [layer(part(0), 1), layer(part(1), 1), layer(part(2), 1)] } });
    expect(layers.has(template, 'transparency')).toBe(true);
    expect(warn).toHaveBeenCalledTimes(1);
    // overpaint alone given: the transparency's two blocks kept, none left: a copy's overpaint whole, warned no more
    layers.set(list, { overpaint: { layers: [{ loci: part(0), color: 0xff0000, clear: false }] } });
    expect(layers.has(template, 'overpaint')).toBe(true);
    const word = (c: number, k: number) =>
      new Uint32Array(layers.image.array.buffer, 0, n.tops.length * CopyLayerKinds)[c * CopyLayerKinds + k];
    // copy 0's overpaint (kind 0) a value, its transparency (kind 1) a block
    expect((word(0, 0) & CopyLayerValue) >>> 0).toBe(CopyLayerValue);
    expect(word(0, 0) & CopyLayerBlock).toBe(0);
    expect(word(0, 1) & CopyLayerBlock).toBe(CopyLayerBlock);
    expect(warn).toHaveBeenCalledTimes(1);
    // laid out again: warned again
    layers.layout(list, [[template, copies]]);
    layers.set(list, { transparency: { layers: [layer(part(0), 1), layer(part(1), 1), layer(part(2), 1)] } });
    expect(warn).toHaveBeenCalledTimes(2);
  });
});

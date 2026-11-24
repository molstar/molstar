/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList, ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { createTemplateCopyMarkers, type TemplateCopies } from '../target/expanded.js';

// WebGL2's copies' markers: a template's particles marked in every copy its objects' own markers, the rest by copy;
// a removal from some copies or from every copy taking the mark from both alike

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

/** two copies (top particles 0, 1) of a template of four particles of one target */
function createNested() {
  const template = createList({ count: 4, coordinates: new Float32Array(12) });
  const list = createList({
    count: 2,
    coordinates: new Float32Array(6),
    targetMapping: new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: template, center: Vec3.create(0, 0, 0) }],
    ]),
  });
  const markers = createTemplateCopyMarkers();
  markers.layout(list, [[template, { tops: [0, 1] } as unknown as TemplateCopies]]);
  const target = ParticleNesting.getLeafList(list, template).targets;
  const targetId = typeof target === 'number' ? target : target[0];
  return { list, template, markers, targetId };
}

describe("a list with templates' copy markers", () => {
  it('moves a mark of every copy into the copies where some of them lose it', () => {
    const { list, template, markers, targetId } = createNested();
    const all = OrderedSet.ofBounds(0, 2);
    // in every copy: its objects' own, no copy marked
    const select = markers.mark(list, ParticleNesting.ofTemplateParticles(list, template, all), MarkerAction.Select);
    expect(select.templateWide.get(template)).toEqual(all);
    expect(markers.isMarked(template, targetId)).toBe(false);
    // one of them taken from the first copy: the mark moved into both copies, the objects' own letting go of it
    const one = markers.mark(list, ParticleNesting.ofInstance(list, 0, 1), MarkerAction.Deselect);
    expect(one.changed).toBe(true);
    expect(one.templateWide.get(template)).toEqual(all);
    expect(markers.isMarked(template, targetId)).toBe(true);
    // taken from every copy then: none left in the copies
    const every = markers.mark(list, ParticleNesting.ofTemplateParticles(list, template, all), MarkerAction.Deselect);
    expect(every.templateWide.get(template)).toEqual(all);
    expect(markers.isMarked(template, targetId)).toBe(false);
  });

  it("takes a mark of every copy from a copy's own marks of the same particles", () => {
    const { list, template, markers, targetId } = createNested();
    markers.mark(list, ParticleNesting.ofInstance(list, 1, 0), MarkerAction.Highlight);
    expect(markers.isMarked(template, targetId)).toBe(true);
    const everyCopy = ParticleNesting.ofTemplateParticles(list, template, OrderedSet.ofSingleton(0));
    markers.mark(list, everyCopy, MarkerAction.RemoveHighlight);
    expect(markers.isMarked(template, targetId)).toBe(false);
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/model/model/location';
import type { ElementIndex, Unit } from '@molstar/model/model/structure';
import type { Features } from '@molstar/model/props/computed/interactions/features';
import type { Bridges, Interactions } from '@molstar/model/props/computed/interactions/interactions';

// Type-only imports and tag tests (as `Interactions.isLocation`): size themes are in every bundle, the interactions
// computation is not.

function isInteractionsLocation(location: Location): location is Interactions.Location {
  return location.kind === 'data-location' && location.tag === 'interactions';
}

function isBridgesLocation(location: Location): location is Bridges.Location {
  return location.kind === 'data-location' && location.tag === 'bridges';
}

function featureSize(
  unit: Unit,
  features: Features,
  index: Features.FeatureIndex,
  elementSize: (unit: Unit, element: ElementIndex) => number,
) {
  return elementSize(unit, unit.elements[features.members[features.offsets[index]]]);
}

/**
 * The size of an interaction (or bridge) location as the interactions mesh bakes it: the smaller `elementSize` of the
 * first member atoms of its features; for a bridge, of all four (its location does not tell the legs apart).
 * Undefined for other locations.
 */
export function getInteractionLocationSize(
  location: Location,
  elementSize: (unit: Unit, element: ElementIndex) => number,
): number | undefined {
  if (isInteractionsLocation(location)) {
    const { unitA, indexA, unitB, indexB } = location.element;
    const { unitsFeatures } = location.data.interactions;
    return Math.min(
      featureSize(unitA, unitsFeatures.get(unitA.id), indexA, elementSize),
      featureSize(unitB, unitsFeatures.get(unitB.id), indexB, elementSize),
    );
  } else if (isBridgesLocation(location)) {
    const { structure, bridges, unitsFeatures } = location.data;
    const b = bridges[location.element.bridgeIndex];
    const uA = structure.unitMap.get(b.unitA);
    const uM = structure.unitMap.get(b.unitM);
    const uB = structure.unitMap.get(b.unitB);
    return Math.min(
      featureSize(uA, unitsFeatures.get(b.unitA), b.indexA, elementSize),
      featureSize(uM, unitsFeatures.get(b.unitM), b.indexMA, elementSize),
      featureSize(uM, unitsFeatures.get(b.unitM), b.indexMB, elementSize),
      featureSize(uB, unitsFeatures.get(b.unitB), b.indexB, elementSize),
    );
  }
  return undefined;
}

/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 */

import type { Iterator } from '@molstar/core/data';
import { NullLocation, type Location } from '@molstar/model/model/location';
import { PositionLocation, isPositionLocation } from '@molstar/model/model/position-location';

export { PositionLocation, isPositionLocation };

export interface LocationValue {
  location: Location;
  location2: Location;
  index: number;
  groupIndex: number;
  instanceIndex: number;
  isSecondary: boolean;
}

export interface LocationIterator extends Iterator<LocationValue> {
  readonly hasNext: boolean;
  readonly isNextNewInstance: boolean;
  readonly groupCount: number;
  readonly instanceCount: number;
  readonly count: number;
  readonly stride: number;
  readonly nonInstanceable: boolean;
  readonly hasLocation2: boolean;
  move(): LocationValue;
  reset(): void;
  skipInstance(): void;
  voidInstances(): void;
  /**
   * The location of instance `instanceIndex`'s first group, as `move` gives it, without moving: per-instance data
   * written by index (`createInstanceColor`); none without groups (`groupCount` 0)
   */
  getInstanceLocation(instanceIndex: number): Location;
  /** Whether instance `instanceIndex`'s first group is secondary, see `getInstanceLocation` */
  isInstanceSecondary(instanceIndex: number): boolean;
}

type LocationGetter = (groupIndex: number, instanceIndex: number) => Location;
type IsSecondaryGetter = (groupIndex: number, instanceIndex: number) => boolean;

/**
 * A class rather than closures: an iterator is created per visual and update (a theme change creates one per unit of a
 * representation), the closures and accessors of an object literal cost more to create than the methods of a class
 */
class LocationIteratorImpl implements LocationIterator {
  private readonly value: LocationValue = {
    location: NullLocation as Location,
    location2: NullLocation as Location,
    index: 0,
    groupIndex: 0,
    instanceIndex: 0,
    isSecondary: false,
  };
  hasNext: boolean;
  isNextNewInstance = false;
  readonly count: number;
  readonly hasLocation2: boolean;
  private groupIndex = 0;
  private instanceIndex = 0;
  private _voidInstances = false;

  constructor(
    readonly groupCount: number,
    readonly instanceCount: number,
    readonly stride: number,
    private readonly getLocation: LocationGetter,
    readonly nonInstanceable: boolean,
    private readonly isSecondary: IsSecondaryGetter,
    private readonly getLocation2: LocationGetter | undefined,
  ) {
    this.count = groupCount * instanceCount;
    this.hasLocation2 = !!getLocation2;
    this.hasNext = 0 < groupCount;
  }

  move() {
    const { value } = this;
    if (this.hasNext) {
      const { groupIndex, instanceIndex, groupCount } = this;
      const instance = this._voidInstances ? -1 : instanceIndex;
      value.groupIndex = groupIndex;
      value.instanceIndex = instanceIndex;
      value.index = instanceIndex * groupCount + groupIndex;
      value.location = this.getLocation(groupIndex, instance);
      if (this.hasLocation2) value.location2 = this.getLocation2!(groupIndex, instance);
      value.isSecondary = this.isSecondary(groupIndex, instance);
      let nextGroup = groupIndex + this.stride;
      if (nextGroup === groupCount) {
        this.instanceIndex = instanceIndex + 1;
        this.isNextNewInstance = true;
        if (instanceIndex + 1 < this.instanceCount) nextGroup = 0;
      } else {
        this.isNextNewInstance = false;
      }
      this.groupIndex = nextGroup;
      this.hasNext = nextGroup < groupCount;
    }
    return value;
  }

  reset() {
    const { value } = this;
    value.location = NullLocation;
    value.location2 = NullLocation;
    value.index = 0;
    value.groupIndex = 0;
    value.instanceIndex = 0;
    value.isSecondary = false;

    this.hasNext = 0 < this.groupCount;
    this.isNextNewInstance = false;
    this.groupIndex = 0;
    this.instanceIndex = 0;
    this._voidInstances = false;
  }

  skipInstance() {
    if (this.hasNext && this.value.instanceIndex === this.instanceIndex) {
      ++this.instanceIndex;
      this.groupIndex = 0;
      this.hasNext = this.instanceIndex < this.instanceCount;
    }
  }

  voidInstances() {
    this._voidInstances = true;
  }

  getInstanceLocation(instanceIndex: number) {
    return this.getLocation(0, this._voidInstances ? -1 : instanceIndex);
  }

  isInstanceSecondary(instanceIndex: number) {
    return this.isSecondary(0, this._voidInstances ? -1 : instanceIndex);
  }
}

const NotSecondary: IsSecondaryGetter = () => false;

export function LocationIterator(
  groupCount: number,
  instanceCount: number,
  stride: number,
  getLocation: LocationGetter,
  nonInstanceable = false,
  isSecondary: IsSecondaryGetter = NotSecondary,
  getLocation2?: LocationGetter,
): LocationIterator {
  if (groupCount % stride !== 0) {
    throw new Error('incompatible groupCount and stride');
  }
  return new LocationIteratorImpl(
    groupCount,
    instanceCount,
    stride,
    getLocation,
    nonInstanceable,
    isSecondary,
    getLocation2,
  );
}

export const EmptyLocationIterator: LocationIterator = {
  get hasNext() {
    return false;
  },
  get isNextNewInstance() {
    return false;
  },
  groupCount: 0,
  instanceCount: 0,
  count: 0,
  stride: 0,
  nonInstanceable: false,
  hasLocation2: false,
  move() {
    return {
      location: NullLocation as Location,
      location2: NullLocation as Location,
      index: 0,
      groupIndex: 0,
      instanceIndex: 0,
      isSecondary: false,
    };
  },
  reset() {},
  getInstanceLocation() {
    return NullLocation as Location;
  },
  isInstanceSecondary() {
    return false;
  },
  skipInstance() {},
  voidInstances() {},
};

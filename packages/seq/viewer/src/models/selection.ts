/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The selection and the highlight channels (_spec/README.md "SelectionModel and HighlightModel"): session state, as
// addresses and Locations, never in the spec. This version holds what a command on an address sets (element keys and
// a Location, and what the selection pass draws of them); the gestures, the anchor, the cursor and the drawing of
// highlights come with the selection's own work.

import { KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import {
  type Elements,
  EmptyElements,
  EmptyLocation,
  elementsEqual,
  type Location,
  locationsEqual,
} from '@molstar/seqx-core/query';

/** What is selected or highlighted: element keys, positions, and the positions drawn (theirs and the elements' extents). */
export interface Emphasized {
  readonly elements: Elements;
  readonly location: Location;
  /** What is drawn: the location and the elements' extents. */
  readonly shown: Location;
  /** The address it was set from, if it was. */
  readonly address?: string;
}

export const Nothing: Emphasized = Object.freeze({
  elements: EmptyElements,
  location: EmptyLocation,
  shown: EmptyLocation,
});

const same = (a: Emphasized, b: Emphasized) =>
  a === b ||
  (a.address === b.address &&
    elementsEqual(a.elements, b.elements) &&
    locationsEqual(a.location, b.location) &&
    locationsEqual(a.shown, b.shown));

export function isNothing(e: Emphasized): boolean {
  return !e.elements.keys.length && !e.location.parts.length && !e.location.rects.length;
}

export class SelectionModel extends ReactiveModel {
  private readonly _value = new State<Emphasized>(Nothing, same);
  /** The selection: an equal value tells nothing. */
  readonly value: ReadonlyState<Emphasized> = this._value.asReadonly;

  set(value: Emphasized): void {
    this._value.set(value);
  }

  clear(): void {
    this._value.set(Nothing);
  }
}

export class HighlightModel extends ReactiveModel {
  private readonly _channels = new KeyedState<string, Emphasized>();
  private readonly _names = new State<readonly string[]>([]);
  /** The channels that hold something, in the order they were first set. */
  readonly names: ReadonlyState<readonly string[]> = this._names.asReadonly;

  /** A channel as a state: only its readers are told of its changes. */
  channel(name: string): ReadonlyState<Emphasized | undefined> {
    return this._channels.key(name);
  }

  get(name: string): Emphasized | undefined {
    return this._channels.get(name);
  }

  set(name: string, value: Emphasized): void {
    const old = this._channels.get(name);
    if (old && same(old, value)) return;
    this._channels.set(name, value);
    if (!this._names.value.includes(name)) this._names.set([...this._names.value, name]);
  }

  clear(name: string): void {
    if (!this._channels.get(name)) return;
    this._channels.delete(name);
    this._names.set(this._names.value.filter((n) => n !== name));
  }
}

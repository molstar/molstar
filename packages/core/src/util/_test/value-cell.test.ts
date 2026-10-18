/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueBox, ValueCell } from '../value-cell.js';

describe('value-cell observers', () => {
  it('flags observers on update and set, until unobserved', () => {
    const cell = ValueCell.create(1);
    const observer: ValueCell.Observer = { changed: false };
    ValueCell.observe(cell, observer);

    ValueCell.updateIfChanged(cell, 1);
    expect(observer.changed).toBe(false);

    ValueCell.update(cell, 2);
    expect(observer.changed).toBe(true);

    observer.changed = false;
    ValueCell.set(cell, ValueBox.create(3));
    expect(observer.changed).toBe(true);

    observer.changed = false;
    ValueCell.unobserve(cell, observer);
    ValueCell.update(cell, 4);
    expect(observer.changed).toBe(false);
  });
});

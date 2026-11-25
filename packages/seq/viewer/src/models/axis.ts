/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An axis's offsets (_spec/README.md "ViewModel"): a viewport is held in offsets along the axis in display order (0
// at the first position shown), so a range of the system, a cyclic range over the origin, a reversed axis and one
// with collapsed pieces read alike. Positions and offsets are fractional; the core's view converts between them
// through the pieces and the direction (one implementation, the renderer's too).

export {
  axesEqual,
  axisBreaks,
  axisKey,
  axisTotal,
  breakHolding,
  cellOffset,
  elementAt,
  hiddenLength,
  nearestOffsetOf,
  offsetOf,
  positionAt,
} from '@molstar/seqx-core/view';

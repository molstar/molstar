/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The square's chunk in the viewer (_spec/README.md "Composition", Navigators, The square: lazy): core's square and
// the set's gestures, loaded by `NavigatorModel` (`loadSquare`) the first time a view shows a square; one chunk.

export { SquareNavigator } from '@molstar/seqx-core/render/navigator-square';
export { SquareSetGesture, setCellReadout, setCursor, setReadout } from './square-set.js';

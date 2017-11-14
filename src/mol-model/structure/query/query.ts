/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../structure'
import Selection from './selection'

// TODO: Query { (s: Structure): Computation<Selection> }

interface Query { (s: Structure): Selection }
export default Query
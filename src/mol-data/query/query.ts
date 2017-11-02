/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../structure'
import Selection from './selection'

interface Query { (s: Structure): Selection }
export default Query
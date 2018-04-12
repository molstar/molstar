/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UI from './ui'
import State from './state'
import * as React from 'react'
import * as ReactDOM from 'react-dom'

const state = new State()
ReactDOM.render(<UI state={ state } />, document.getElementById('app'));
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import * as ReactDOM from 'react-dom'

import './index.html'

import { App } from './app';
import { AppComponent } from './component/app';
import { urlQueryParameter } from 'mol-util/url-query';

const elm = document.getElementById('app') as HTMLElement
if (!elm) throw new Error('Can not find element with id "app".')

const app = new App()
ReactDOM.render(React.createElement(AppComponent, { app }), elm);

const assemblyId = urlQueryParameter('assembly')
const pdbId = urlQueryParameter('pdb')
if (pdbId) app.loadPdbId(pdbId, assemblyId)
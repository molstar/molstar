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
if (pdbId) app.loadPdbIdOrMmcifUrl(pdbId, { assemblyId })

// app.loadPdbIdOrMmcifUrl('3pqr')
// app.loadCcp4Url('http://localhost:8091/ngl/data/3pqr-mode0.ccp4')

// app.loadPdbIdOrMmcifUrl('1lee')
// app.loadCcp4Url('http://localhost:8091/ngl/data/1lee.ccp4')

// app.loadPdbIdOrMmcifUrl('6DRV')
// app.loadCcp4Url('http://localhost:8091/ngl/data/betaGal.mrc')

// app.loadPdbIdOrMmcifUrl('3pqr')
// app.loadVolCifUrl('https://webchem.ncbr.muni.cz/DensityServer/x-ray/3pqr/cell?space=fractional', true)

// app.loadPdbIdOrMmcifUrl('5ire')
// app.loadVolcifUrl('https://webchem.ncbr.muni.cz/DensityServer/em/emd-8116/cell?space=cartesian&detail=6', true)

// app.loadPdbIdOrMmcifUrl('5gag')
// app.loadVolcifUrl('https://webchem.ncbr.muni.cz/DensityServer/em/emd-8003/cell?detail=3', true)

// app.loadPdbIdOrMmcifUrl('http://localhost:8091/test/pdb-dev/carb/1B5F-carb.cif')
// app.loadPdbIdOrMmcifUrl('http://localhost:8091/test/pdb-dev/carb/2HYV-carb.cif')
// app.loadPdbIdOrMmcifUrl('http://localhost:8091/test/pdb-dev/carb/2WMG-carb.cif')
// app.loadPdbIdOrMmcifUrl('http://localhost:8091/test/pdb-dev/carb/5KDS-carb.cif')
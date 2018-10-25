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

app.loadPdbIdOrMmcifUrl('6DRV')
app.loadCcp4Url('http://localhost:8091/ngl/data/betaGal.mrc')

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

const basisX = [
    h.xlen,
    0,
    0
  ]

  const basisY = [
    h.ylen * Math.cos(Math.PI / 180.0 * h.gamma),
    h.ylen * Math.sin(Math.PI / 180.0 * h.gamma),
    0
  ]

  const basisZ = [
    h.zlen * Math.cos(Math.PI / 180.0 * h.beta),
    h.zlen * (
      Math.cos(Math.PI / 180.0 * h.alpha) -
      Math.cos(Math.PI / 180.0 * h.gamma) *
      Math.cos(Math.PI / 180.0 * h.beta)
    ) / Math.sin(Math.PI / 180.0 * h.gamma),
    0
  ]
  basisZ[ 2 ] = Math.sqrt(
    h.zlen * h.zlen * Math.sin(Math.PI / 180.0 * h.beta) *
    Math.sin(Math.PI / 180.0 * h.beta) - basisZ[ 1 ] * basisZ[ 1 ]
  )

  const basis = [ 0, basisX, basisY, basisZ ]
  const nxyz = [ 0, h.MX, h.MY, h.MZ ]
  const mapcrs = [ 0, h.MAPC, h.MAPR, h.MAPS ]

  const matrix = new Matrix4()

  matrix.set(
    basis[ mapcrs[1] ][0] / nxyz[ mapcrs[1] ],
    basis[ mapcrs[2] ][0] / nxyz[ mapcrs[2] ],
    basis[ mapcrs[3] ][0] / nxyz[ mapcrs[3] ],
    0,
    basis[ mapcrs[1] ][1] / nxyz[ mapcrs[1] ],
    basis[ mapcrs[2] ][1] / nxyz[ mapcrs[2] ],
    basis[ mapcrs[3] ][1] / nxyz[ mapcrs[3] ],
    0,
    basis[ mapcrs[1] ][2] / nxyz[ mapcrs[1] ],
    basis[ mapcrs[2] ][2] / nxyz[ mapcrs[2] ],
    basis[ mapcrs[3] ][2] / nxyz[ mapcrs[3] ],
    0,
    0, 0, 0, 1
  )
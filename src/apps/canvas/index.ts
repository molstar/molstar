/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'

import Viewer from 'mol-view/viewer';
import { EveryLoci } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/util/marker-data';
import { labelFirst } from 'mol-view/label';
import { getCifFromUrl, getModelFromMmcif } from './util';
import { SymmetryView } from './assembly-symmetry';

const container = document.getElementById('container')
if (!container) throw new Error('Can not find element with id "container".')

const canvas = document.getElementById('canvas') as HTMLCanvasElement
if (!canvas) throw new Error('Can not find element with id "canvas".')

const info = document.getElementById('info')
if (!info) throw new Error('Can not find element with id "info".')

const viewer = Viewer.create(canvas, container)
viewer.animate()

viewer.input.resize.subscribe(() => {
    // do whatever appropriate
})

viewer.input.move.subscribe(({x, y, inside, buttons}) => {
    if (!inside || buttons) return
    const p = viewer.identify(x, y)
    const loci = viewer.getLoci(p)

    viewer.mark(EveryLoci, MarkerAction.RemoveHighlight)
    viewer.mark(loci, MarkerAction.Highlight)

    const label = labelFirst(loci)
    info.innerText = `${label}`
})

async function init() {
    const assembly = '1'

    const cif = await getCifFromUrl('https://files.rcsb.org/download/4hhb.cif')
    const model = await getModelFromMmcif(cif)
    
    const view = await SymmetryView(model, assembly)
    viewer.center(view.structure.boundary.sphere.center)
    viewer.add(view.cartoon)
    viewer.add(view.ballAndStick)
    viewer.add(view.axes)

    // ensure the added representations get rendered, i.e. without mouse input
    viewer.requestDraw()
}

init()
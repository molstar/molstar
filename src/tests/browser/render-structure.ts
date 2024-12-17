/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { Canvas3D, Canvas3DContext } from '../../mol-canvas3d/canvas3d';
import { CIF, CifFrame } from '../../mol-io/reader/cif';
import { Model, Structure } from '../../mol-model/structure';
import { ColorTheme } from '../../mol-theme/color';
import { SizeTheme } from '../../mol-theme/size';
import { CartoonRepresentationProvider } from '../../mol-repr/structure/representation/cartoon';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { MolecularSurfaceRepresentationProvider } from '../../mol-repr/structure/representation/molecular-surface';
import { BallAndStickRepresentationProvider } from '../../mol-repr/structure/representation/ball-and-stick';
import { GaussianSurfaceRepresentationProvider } from '../../mol-repr/structure/representation/gaussian-surface';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Representation } from '../../mol-repr/representation';
import { throttleTime } from 'rxjs/operators';
import { MarkerAction } from '../../mol-util/marker-action';
import { EveryLoci } from '../../mol-model/loci';
import { lociLabel } from '../../mol-theme/label';
import { InteractionsRepresentationProvider } from '../../mol-model-props/computed/representations/interactions';
import { InteractionsProvider } from '../../mol-model-props/computed/interactions';
import { SecondaryStructureProvider } from '../../mol-model-props/computed/secondary-structure';
import { SyncRuntimeContext } from '../../mol-task/execution/synchronous';
import { AssetManager } from '../../mol-util/assets';
import { MembraneOrientationProvider } from '../../extensions/anvil/prop';
import { MembraneOrientationRepresentationProvider } from '../../extensions/anvil/representation';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);

const assetManager = new AssetManager();

const canvas3dContext = Canvas3DContext.fromCanvas(canvas, assetManager);
const canvas3d = Canvas3D.create(canvas3dContext);
resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
canvas3dContext.syncPixelScale();
canvas3d.requestResize();
canvas3d.animate();

const info = document.createElement('div');
info.style.position = 'absolute';
info.style.fontFamily = 'sans-serif';
info.style.fontSize = '16pt';
info.style.bottom = '20px';
info.style.right = '20px';
info.style.color = 'white';
parent.appendChild(info);

let prevReprLoci = Representation.Loci.Empty;
canvas3d.input.move.pipe(throttleTime(100)).subscribe(({ x, y }) => {
    const pickingId = canvas3d.identify(x, y)?.id;
    let label = '';
    if (pickingId) {
        const reprLoci = canvas3d.getLoci(pickingId);
        label = lociLabel(reprLoci.loci);
        if (!Representation.Loci.areEqual(prevReprLoci, reprLoci)) {
            canvas3d.mark(prevReprLoci, MarkerAction.RemoveHighlight);
            canvas3d.mark(reprLoci, MarkerAction.Highlight);
            prevReprLoci = reprLoci;
        }
    } else {
        canvas3d.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
        prevReprLoci = Representation.Loci.Empty;
    }
    info.innerHTML = label;
});

canvas3d.input.resize.subscribe(() => {
    resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
    canvas3dContext.syncPixelScale();
    canvas3d.requestResize();
});

async function parseCif(data: string | Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}

async function downloadFromPdb(pdb: string) {
    const parsed = await downloadCif(`https://models.rcsb.org/${pdb}.bcif`, true);
    return parsed.blocks[0];
}

async function getModels(frame: CifFrame) {
    return await trajectoryFromMmCIF(frame).run();
}

async function getStructure(model: Model) {
    return Structure.ofModel(model);
}

const reprCtx = {
    webgl: canvas3d.webgl,
    colorThemeRegistry: ColorTheme.createRegistry(),
    sizeThemeRegistry: SizeTheme.createRegistry()
};
function getCartoonRepr() {
    return CartoonRepresentationProvider.factory(reprCtx, CartoonRepresentationProvider.getParams);
}

function getInteractionRepr() {
    return InteractionsRepresentationProvider.factory(reprCtx, InteractionsRepresentationProvider.getParams);
}

function getBallAndStickRepr() {
    return BallAndStickRepresentationProvider.factory(reprCtx, BallAndStickRepresentationProvider.getParams);
}

function getMolecularSurfaceRepr() {
    return MolecularSurfaceRepresentationProvider.factory(reprCtx, MolecularSurfaceRepresentationProvider.getParams);
}

function getGaussianSurfaceRepr() {
    return GaussianSurfaceRepresentationProvider.factory(reprCtx, GaussianSurfaceRepresentationProvider.getParams);
}

function getMembraneOrientationRepr() {
    return MembraneOrientationRepresentationProvider.factory(reprCtx, MembraneOrientationRepresentationProvider.getParams);
}

async function init() {
    const ctx = { runtime: SyncRuntimeContext, assetManager };

    const cif = await downloadFromPdb('3pqr');
    const models = await getModels(cif);
    const structure = await getStructure(models.representative);

    console.time('compute SecondaryStructure');
    await SecondaryStructureProvider.attach(ctx, structure);
    console.timeEnd('compute SecondaryStructure');

    console.time('compute Membrane Orientation');
    await MembraneOrientationProvider.attach(ctx, structure);
    console.timeEnd('compute Membrane Orientation');

    console.time('compute Interactions');
    await InteractionsProvider.attach(ctx, structure);
    console.timeEnd('compute Interactions');
    console.log(InteractionsProvider.get(structure).value);

    const show = {
        cartoon: true,
        interaction: true,
        ballAndStick: true,
        molecularSurface: false,
        gaussianSurface: false,
        membrane: true
    };

    const cartoonRepr = getCartoonRepr();
    const interactionRepr = getInteractionRepr();
    const ballAndStickRepr = getBallAndStickRepr();
    const molecularSurfaceRepr = getMolecularSurfaceRepr();
    const gaussianSurfaceRepr = getGaussianSurfaceRepr();
    const membraneOrientationRepr = getMembraneOrientationRepr();

    if (show.cartoon) {
        cartoonRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('element-symbol', { structure }),
            size: reprCtx.sizeThemeRegistry.create('uniform', { structure })
        });
        await cartoonRepr.createOrUpdate({ ...CartoonRepresentationProvider.defaultValues, quality: 'auto' }, structure).run();
    }

    if (show.interaction) {
        interactionRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('interaction-type', { structure }),
            size: reprCtx.sizeThemeRegistry.create('uniform', { structure })
        });
        await interactionRepr.createOrUpdate({ ...InteractionsRepresentationProvider.defaultValues, quality: 'auto' }, structure).run();
    }

    if (show.ballAndStick) {
        ballAndStickRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('element-symbol', { structure }),
            size: reprCtx.sizeThemeRegistry.create('uniform', { structure }, { value: 1 })
        });
        await ballAndStickRepr.createOrUpdate({ ...BallAndStickRepresentationProvider.defaultValues, quality: 'auto' }, structure).run();
    }

    if (show.molecularSurface) {
        molecularSurfaceRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('physical', { structure })
        });
        console.time('molecular surface');
        await molecularSurfaceRepr.createOrUpdate({ ...MolecularSurfaceRepresentationProvider.defaultValues, quality: 'custom', alpha: 0.5, flatShaded: true, doubleSided: true, resolution: 0.3 }, structure).run();
        console.timeEnd('molecular surface');
    }

    if (show.gaussianSurface) {
        gaussianSurfaceRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('physical', { structure })
        });
        console.time('gaussian surface');
        await gaussianSurfaceRepr.createOrUpdate({ ...GaussianSurfaceRepresentationProvider.defaultValues, quality: 'custom', alpha: 1.0, flatShaded: true, doubleSided: true, resolution: 0.3 }, structure).run();
        console.timeEnd('gaussian surface');
    }

    if (show.membrane) {
        await membraneOrientationRepr.createOrUpdate({ ...MembraneOrientationRepresentationProvider.defaultValues, quality: 'auto' }, structure).run();
    }

    if (show.cartoon) canvas3d.add(cartoonRepr);
    if (show.interaction) canvas3d.add(interactionRepr);
    if (show.ballAndStick) canvas3d.add(ballAndStickRepr);
    if (show.molecularSurface) canvas3d.add(molecularSurfaceRepr);
    if (show.gaussianSurface) canvas3d.add(gaussianSurfaceRepr);
    if (show.membrane) canvas3d.add(membraneOrientationRepr);
    canvas3d.requestCameraReset();
    // canvas3d.setProps({ trackball: { ...canvas3d.props.trackball, spin: true } })
}

init();
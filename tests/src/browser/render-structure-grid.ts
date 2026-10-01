/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Herman Bergwerf <post@hbergwerf.nl>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssetManager } from '@molstar/core/util/assets';
import { Canvas3D, type Canvas3DProps, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { resizeCanvas } from '@molstar/graphics/canvas3d/util';
import { parseSdf as _parseSdf, type SdfFile } from '@molstar/io/reader/sdf/parser';
import { Structure } from '@molstar/model/model/structure';
import { trajectoryFromSdf } from '@molstar/model/formats/structure/sdf';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import type { RepresentationContext } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { StructureRepresentationProvider } from '@molstar/graphics/repr/structure/representation';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

async function downloadPubChemSdf(cid: number) {
    const root = 'https://pubchem.ncbi.nlm.nih.gov/rest';
    const r = await fetch(`${root}/pug/compound/cid/${cid}/sdf?record_type=3d`);
    return await r.text();
}

async function parseSdf(sdf: string): Promise<SdfFile> {
    const parsed = await _parseSdf(sdf).run();
    if (parsed.isError) {
        throw parsed;
    } else {
        return parsed.result;
    }
}

function getRepr<P extends PD.Params>(provider: StructureRepresentationProvider<P>, reprCtx: RepresentationContext) {
    return provider.factory(reprCtx, provider.getParams);
}

function main() {
    const parent = document.body;
    const canvas = document.createElement('canvas');
    parent.appendChild(canvas);

    const canvas3ds: Canvas3D[] = [];

    const canvas3dContext = Canvas3DContext.fromCanvas(canvas, new AssetManager());
    resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
    canvas3dContext.syncPixelScale();

    canvas3dContext.input.resize.subscribe(() => {
        resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
        canvas3dContext.syncPixelScale();
        for (const canvas3d of canvas3ds) canvas3d.requestResize();
    });

    const cids = [
        6440397,
        2244,
        2519,
        241
    ];

    const reprCtx: RepresentationContext = {
        webgl: canvas3dContext.webgl,
        colorThemeRegistry: ColorTheme.createRegistry(),
        sizeThemeRegistry: SizeTheme.createRegistry()
    };

    for (let i = 0; i < cids.length; i++) {
        const ix = i % 2;
        const iy = Math.floor(i / 2);
        addViewer(canvas3dContext, reprCtx, cids[i], {
            name: 'relative-frame',
            params: {
                x: ix / 2,
                y: iy / 2,
                width: 0.5,
                height: 0.5,
            }
        }, {
            adjustCylinderLength: ix === 0
        }).then(canvas3d => {
            canvas3ds.push(canvas3d);
        });
    }
}

async function addViewer(
    ctx: Canvas3DContext,
    reprCtx: RepresentationContext,
    cid: number,
    viewport: Canvas3DProps['viewport'],
    reprValues: Partial<typeof BallAndStickRepresentationProvider.defaultValues>
) {
    const file = await parseSdf(await downloadPubChemSdf(cid));
    const models = await trajectoryFromSdf(file.compounds[0]).run();
    const structure = Structure.ofModel(models.representative);
    const canvas3d = Canvas3D.create(ctx, { viewport });
    canvas3d.requestResize();

    const repr = getRepr(BallAndStickRepresentationProvider, reprCtx);
    repr.setTheme({
        color: reprCtx.colorThemeRegistry.create('element-symbol', { structure }, {
            carbonColor: { name: 'element-symbol' }
        } as ColorTheme.BuiltInParams<'element-symbol'>),
        size: reprCtx.sizeThemeRegistry.create('physical', { structure })
    });

    await repr.createOrUpdate({
        ...BallAndStickRepresentationProvider.defaultValues,
        ...reprValues
    }, structure).run();

    canvas3d.add(repr);
    canvas3d.animate();

    return canvas3d;
}

main();

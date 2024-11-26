/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import express from 'express';
import { MembraneServerConfig } from './config';
import { swaggerUiAssetsHandler, swaggerUiIndexHandler } from '../common/swagger-ui';
import { getSchema, shortcutIconLink } from './web-schema';
import { MembraneOrientationProvider } from '../../extensions/anvil/prop';
import { SyncRuntimeContext } from '../../mol-task/execution/synchronous';
import { AssetManager } from '../../mol-util/assets';
import { CIF, CifFrame } from '../../mol-io/reader/cif';
import { Model, Structure, StructureSymmetry } from '../../mol-model/structure';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { ANVILParams, ANVILProps } from '../../extensions/anvil/algorithm';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ConsoleLogger } from '../../mol-util/console-logger';

const assetManager = new AssetManager();

export function initWebApi(app: express.Express) {
    function makePath(p: string) {
        return MembraneServerConfig.apiPrefix + '/' + p;
    }

    app.get(makePath('predict/:id/'), async (req, res) => predictMembraneOrientation(req, res));

    app.get(makePath('openapi.json'), (_, res) => {
        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With'
        });
        res.end(JSON.stringify(getSchema()));
    });

    app.use(makePath(''), swaggerUiAssetsHandler());
    app.get(makePath(''), swaggerUiIndexHandler({
        openapiJsonUrl: makePath('openapi.json'),
        apiPrefix: MembraneServerConfig.apiPrefix,
        title: 'MembraneServer API',
        shortcutIconLink
    }));
}

async function predictMembraneOrientation(req: express.Request, res: express.Response) {
    try {
        const ctx = { runtime: SyncRuntimeContext, assetManager };

        const entryId = req.params.id;
        const assemblyId = req.query.assemblyId as string ?? '1';
        const p = parseParams(req);
        ConsoleLogger.log('predictMembraneOrientation', `${entryId}-${assemblyId} with params: ${JSON.stringify(p)}`);

        const cif = await downloadFromPdb(entryId);
        const models = await getModels(cif);
        const structure = await getStructure(models.representative, assemblyId);

        await MembraneOrientationProvider.attach(ctx, structure, p);
        const data = MembraneOrientationProvider.get(structure).value;

        res.status(200).json(data);
    } catch (e) {
        const error = 'Failed to compute membrane orientation';
        ConsoleLogger.error(error, e);
        res.status(500).json({ error });
    }
}

const defaults = PD.getDefaultValues(ANVILParams);
function parseParams(req: express.Request): ANVILProps {
    const {
        numberOfSpherePoints = defaults.numberOfSpherePoints,
        stepSize = defaults.stepSize,
        minThickness = defaults.minThickness,
        maxThickness = defaults.maxThickness,
        asaCutoff = defaults.asaCutoff,
        adjust = defaults.adjust,
        tmdetDefinition = defaults.tmdetDefinition,
    } = req.query;
    return {
        numberOfSpherePoints: Number(numberOfSpherePoints),
        stepSize: Number(stepSize),
        minThickness: Number(minThickness),
        maxThickness: Number(maxThickness),
        asaCutoff: Number(asaCutoff),
        adjust: Number(adjust),
        tmdetDefinition: tmdetDefinition === 'true',
    };
}

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
    const parsed = await downloadCif(MembraneServerConfig.bcifSource(pdb), true);
    return parsed.blocks[0];
}

async function getModels(frame: CifFrame) {
    return await trajectoryFromMmCIF(frame).run();
}

async function getStructure(model: Model, assemblyId: string) {
    const modelStructure = Structure.ofModel(model);
    return await StructureSymmetry.buildAssembly(modelStructure, assemblyId).run();
}
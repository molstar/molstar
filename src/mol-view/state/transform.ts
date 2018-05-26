/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif'
import { FileEntity, DataEntity, UrlEntity, CifEntity, MmcifEntity, ModelEntity, StructureEntity, SpacefillEntity, AnyEntity, NullEntity } from './entity';
import { Model, Structure } from 'mol-model/structure';

import { StateContext } from './context';
import Spacefill, { SpacefillProps } from 'mol-geo/representation/structure/spacefill';
import { StructureRepresentation } from 'mol-geo/representation/structure';
import StructureSymmetry from 'mol-model/structure/structure/symmetry';

type transformer<I extends AnyEntity, O extends AnyEntity, P extends {}> = (ctx: StateContext, inputEntity: I, props?: P) => Promise<O>

export interface StateTransform<I extends AnyEntity, O extends AnyEntity, P extends {}> {
    inputKind: I['kind']
    outputKind: O['kind']
    kind: string
    apply: transformer<I, O, P>
}

export namespace StateTransform {
    export function create<I extends AnyEntity, O extends AnyEntity, P extends {}>(inputKind: I['kind'], outputKind: O['kind'], kind: string, transformer: transformer<I, O, P>) {
        return { inputKind, outputKind, kind, apply: transformer }
    }
}

export type AnyTransform = StateTransform<AnyEntity, AnyEntity, {}>

export type UrlToData = StateTransform<UrlEntity, DataEntity, {}>
export const UrlToData: UrlToData = StateTransform.create('url', 'data', 'url-to-data',
    async function (ctx: StateContext, urlEntity: UrlEntity) {
        return DataEntity.ofData(ctx, await urlEntity.value.getData(), urlEntity.value.type)
    })

export type FileToData = StateTransform<FileEntity, DataEntity, {}>
export const FileToData: FileToData = StateTransform.create('file', 'data', 'file-to-data',
    async function (ctx: StateContext, fileEntity: FileEntity) {
        return DataEntity.ofData(ctx, await fileEntity.value.getData(), fileEntity.value.type)
    })

export type DataToCif = StateTransform<DataEntity, CifEntity, {}>
export const DataToCif: DataToCif = StateTransform.create('data', 'cif', 'data-to-cif',
    async function (ctx: StateContext, dataEntity: DataEntity) {
        const comp = CIF.parse(dataEntity.value.data)
        const parsed = await comp.run(ctx.log)
        if (parsed.isError) throw parsed
        return CifEntity.ofCifFile(ctx, parsed.result)
    })

export type CifToMmcif = StateTransform<CifEntity, MmcifEntity, {}>
export const CifToMmcif: CifToMmcif = StateTransform.create('cif', 'mmcif', 'cif-to-mmcif',
    async function (ctx: StateContext, cifEntity: CifEntity) {
        return MmcifEntity.ofMmcifDb(ctx, CIF.schema.mmCIF(cifEntity.value.blocks[0]))
    })

export type MmcifToModel = StateTransform<MmcifEntity, ModelEntity, {}>
export const MmcifToModel: MmcifToModel = StateTransform.create('mmcif', 'model', 'mmcif-to-model',
    async function (ctx: StateContext, mmcifEntity: MmcifEntity) {
        const models = await Model.create({ kind: 'mmCIF', data: mmcifEntity.value }).run(ctx.log)
        return ModelEntity.ofModels(ctx, models)
    })

export interface StructureProps {
    assembly?: string
}

export type ModelToStructure = StateTransform<ModelEntity, StructureEntity, StructureProps>
export const ModelToStructure: ModelToStructure = StateTransform.create('model', 'structure', 'model-to-structure',
    async function (ctx: StateContext, modelEntity: ModelEntity, props: StructureProps = {}) {
        const model = modelEntity.value[0]
        const assembly = props.assembly
        let structure: Structure
        const assemblies = model.symmetry.assemblies
        if (assemblies.length) {
            structure = await StructureSymmetry.buildAssembly(Structure.ofModel(model), assembly || '1').run(ctx.log)
        } else {
            structure = Structure.ofModel(model)
        }
        return StructureEntity.ofStructure(ctx, structure)
    })

export type StructureCenter = StateTransform<StructureEntity, NullEntity, {}>
    export const StructureCenter: StructureCenter = StateTransform.create('structure', 'null', 'structure-center',
        async function (ctx: StateContext, structureEntity: StructureEntity) {
            ctx.viewer.center(structureEntity.value.boundary.sphere.center)
            return NullEntity
        })

export type StructureToSpacefill = StateTransform<StructureEntity, SpacefillEntity, SpacefillProps>
export const StructureToSpacefill: StructureToSpacefill = StateTransform.create('structure', 'spacefill', 'structure-to-spacefill',
    async function (ctx: StateContext, structureEntity: StructureEntity, props: SpacefillProps = {}) {
        const spacefillRepr = StructureRepresentation(Spacefill)
        await spacefillRepr.create(structureEntity.value, props).run(ctx.log)
        ctx.viewer.add(spacefillRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        // ctx.viewer.input.drag.subscribe(async () => {
        //     console.log('drag')
        //     console.time('spacefill update')
        //     await spacefillRepr.update(props).run(ctx.log)
        //     console.timeEnd('spacefill update')
        //     ctx.viewer.add(spacefillRepr)
        //     ctx.viewer.update()
        //     ctx.viewer.requestDraw()
        // })
        return SpacefillEntity.ofRepr(ctx, spacefillRepr)
    })

export type SpacefillUpdate = StateTransform<SpacefillEntity, NullEntity, SpacefillProps>
export const SpacefillUpdate: SpacefillUpdate = StateTransform.create('spacefill', 'null', 'spacefill-update',
    async function (ctx: StateContext, spacefillEntity: SpacefillEntity, props: SpacefillProps = {}) {
        const spacefillRepr = spacefillEntity.value
        await spacefillRepr.update(props).run(ctx.log)
        ctx.viewer.add(spacefillRepr)
        ctx.viewer.update()
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return NullEntity
    })

// composed

export type MmcifUrlToModel = StateTransform<UrlEntity, ModelEntity, {}>
export const MmcifUrlToModel: MmcifUrlToModel = StateTransform.create('url', 'model', 'url-to-model',
    async function (ctx: StateContext, urlEntity: UrlEntity) {
        const dataEntity = await UrlToData.apply(ctx, urlEntity)
        return DataToModel.apply(ctx, dataEntity)
    })

export type MmcifFileToModel = StateTransform<FileEntity, ModelEntity, {}>
export const MmcifFileToModel: MmcifFileToModel = StateTransform.create('file', 'model', 'file-to-model',
    async function (ctx: StateContext, fileEntity: FileEntity) {
        const dataEntity = await FileToData.apply(ctx, fileEntity)
        return DataToModel.apply(ctx, dataEntity)
    })

export type DataToModel = StateTransform<DataEntity, ModelEntity, {}>
export const DataToModel: DataToModel = StateTransform.create('data', 'model', 'data-to-model',
    async function getModelFromData(ctx: StateContext, dataEntity: DataEntity) {
        const cifEntity = await DataToCif.apply(ctx, dataEntity)
        const mmcifEntity = await CifToMmcif.apply(ctx, cifEntity)
        return MmcifToModel.apply(ctx, mmcifEntity)
    })

export type ModelToSpacefill = StateTransform<ModelEntity, SpacefillEntity, SpacefillProps>
export const ModelToSpacefill: ModelToSpacefill = StateTransform.create('model', 'spacefill', 'model-to-spacefill',
    async function (ctx: StateContext, modelEntity: ModelEntity, props: SpacefillProps = {}) {
        const structureEntity = await ModelToStructure.apply(ctx, modelEntity)
        return StructureToSpacefill.apply(ctx, structureEntity, props)
    })

export type MmcifUrlToSpacefill = StateTransform<UrlEntity, SpacefillEntity, SpacefillProps>
export const MmcifUrlToSpacefill: MmcifUrlToSpacefill = StateTransform.create('url', 'spacefill', 'url-to-spacefill',
    async function (ctx: StateContext, urlEntity: UrlEntity, props: SpacefillProps = {}) {
        const modelEntity = await MmcifUrlToModel.apply(ctx, urlEntity)
        return ModelToSpacefill.apply(ctx, modelEntity, props)
    })

export type MmcifFileToSpacefill = StateTransform<FileEntity, SpacefillEntity, SpacefillProps>
export const MmcifFileToSpacefill: MmcifFileToSpacefill = StateTransform.create('file', 'spacefill', 'file-to-spacefill',
    async function (ctx: StateContext, fileEntity: FileEntity, props: SpacefillProps = {}) {
        const modelEntity = await MmcifFileToModel.apply(ctx, fileEntity)
        return ModelToSpacefill.apply(ctx, modelEntity, props)
    })
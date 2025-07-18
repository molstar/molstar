/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { CIF, CifBlock, CifCategory, CifFile } from '../../../mol-io/reader/cif';
import { toTable } from '../../../mol-io/reader/cif/schema';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Model } from '../../../mol-model/structure';
import { Structure, StructureElement } from '../../../mol-model/structure/structure';
import { Asset } from '../../../mol-util/assets';
import { Jsonable, canonicalJsonString } from '../../../mol-util/json';
import { objectOfArraysToArrayOfObjects, pickObjectKeysWithRemapping, promiseAllObj } from '../../../mol-util/object';
import { Choice } from '../../../mol-util/param-choice';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { AtomRanges } from '../helpers/atom-ranges';
import { IndicesAndSortings } from '../helpers/indexing';
import { MaybeStringParamDefinition } from '../helpers/param-definition';
import { MVSAnnotationRow, MVSAnnotationSchema, getCifAnnotationSchema } from '../helpers/schemas';
import { getAtomRangesForRow } from '../helpers/selections';
import { Maybe, isDefined, safePromise } from '../helpers/utils';


/** Allowed values for the annotation format parameter */
const MVSAnnotationFormat = new Choice({ json: 'json', cif: 'cif', bcif: 'bcif' }, 'json');
type MVSAnnotationFormat = Choice.Values<typeof MVSAnnotationFormat>
const MVSAnnotationFormatTypes = { json: 'string', cif: 'string', bcif: 'binary' } as const satisfies { [format in MVSAnnotationFormat]: 'string' | 'binary' };

/** Parameter definition for custom model property "MVS Annotations" */
export type MVSAnnotationsParams = typeof MVSAnnotationsParams
export const MVSAnnotationsParams = {
    annotations: PD.ObjectList(
        {
            source: PD.MappedStatic('source-cif', {
                'source-cif': PD.EmptyGroup(),
                'url': PD.Group({
                    url: PD.Text(''),
                    format: MVSAnnotationFormat.PDSelect(),
                }),
            }),
            schema: MVSAnnotationSchema.PDSelect(),
            cifBlock: PD.MappedStatic('index', {
                index: PD.Group({ index: PD.Numeric(0, { min: 0, step: 1 }, { description: '0-based index of the block' }) }),
                header: PD.Group({ header: PD.Text(undefined, { description: 'Block header' }) }),
            }, { description: 'Specify which CIF block contains annotation data (only relevant when format=cif or format=bcif)' }),
            cifCategory: MaybeStringParamDefinition({ placeholder: 'Take first category', description: 'Specify which CIF category contains annotation data (only relevant when format=cif or format=bcif)' }),
            fieldRemapping: PD.ObjectList({
                standardName: PD.Text('', { placeholder: ' ', description: 'Standard name of the selector field (e.g. label_asym_id)' }),
                actualName: MaybeStringParamDefinition({ placeholder: 'Ignore field', description: 'Actual name of the field in the annotation data (e.g. spam_chain_id), null to ignore the field with standard name' }),
            }, e => `"${e.standardName}": ${e.actualName === null ? 'null' : `"${e.actualName}"`}`, { description: 'Optional remapping of annotation field names { standardName1: actualName1, ... }. Use { "label_asym_id": "X" } to load actual field "X" as "label_asym_id". Use { "label_asym_id": null } to ignore actual field "label_asym_id". Fields not mentioned here are mapped implicitely (i.e. actual name = standard name).' }),
            id: PD.Text('', { description: 'Arbitrary identifier that can be referenced by MVSAnnotationColorTheme' }),
        },
        obj => obj.id
    ),
};

/** Parameter values for custom model property "MVS Annotations" */
export type MVSAnnotationsProps = PD.Values<MVSAnnotationsParams>

/** Parameter values for a single annotation within custom model property "MVS Annotations" */
export type MVSAnnotationSpec = MVSAnnotationsProps['annotations'][number]

/** Describes the source of an annotation file */
type MVSAnnotationSource = { kind: 'url', url: string, format: MVSAnnotationFormat } | { kind: 'source-cif' }

/** Data file with one or more (in case of CIF) annotations */
type MVSAnnotationFile = { format: 'json', data: Jsonable } | { format: 'cif', data: CifFile }

/** Data for a single annotation */
type MVSAnnotationData = { format: 'json', data: Jsonable } | { format: 'cif', data: CifCategory }


/** Provider for custom model property "Annotations" */
export const MVSAnnotationsProvider: CustomModelProperty.Provider<MVSAnnotationsParams, MVSAnnotations> = CustomModelProperty.createProvider({
    label: 'MVS Annotations',
    descriptor: CustomPropertyDescriptor({
        name: 'mvs-annotations',
    }),
    type: 'static',
    defaultParams: MVSAnnotationsParams,
    getParams: (data: Model) => MVSAnnotationsParams,
    isApplicable: (data: Model) => true,
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<MVSAnnotationsProps>) => {
        props = { ...PD.getDefaultValues(MVSAnnotationsParams), ...props };
        const specs: MVSAnnotationSpec[] = props.annotations ?? [];
        const annots = await MVSAnnotations.fromSpecs(ctx, specs, data);
        return { value: annots } satisfies CustomProperty.Data<MVSAnnotations>;
    }
});


/** Represents multiple annotations retrievable by their ID */
export class MVSAnnotations {
    private constructor(private dict: { [id: string]: MVSAnnotation }) { }
    static async fromSpecs(ctx: CustomProperty.Context, specs: MVSAnnotationSpec[], model?: Model): Promise<MVSAnnotations> {
        const sources: MVSAnnotationSource[] = specs.map(annotationSourceFromSpec);
        const files = await getFilesFromSources(ctx, sources, model);
        const annots: { [id: string]: MVSAnnotation } = {};
        for (let i = 0; i < specs.length; i++) {
            const spec = specs[i];
            try {
                const file = files[i];
                if (!file.ok) throw file.error;
                annots[spec.id] = await MVSAnnotation.fromSpec(ctx, spec, file.value);
            } catch (err) {
                ctx.errorContext?.add('mvs', `Failed to obtain annotation (${err}).\nAnnotation specification source params: ${JSON.stringify(spec.source.params)}`);
                console.error(`Failed to obtain annotation (${err}).\nAnnotation specification:`, spec);
                annots[spec.id] = MVSAnnotation.createEmpty(spec.schema);
            }
        }
        return new MVSAnnotations(annots);
    }
    getAnnotation(id: string): MVSAnnotation | undefined {
        return this.dict[id];
    }
    getAllAnnotations(): MVSAnnotation[] {
        return Object.values(this.dict);
    }
}


/** Retrieve annotation with given `annotationId` from custom model property "MVS Annotations" and the model from which it comes */
export function getMVSAnnotationForStructure(structure: Structure, annotationId: string): { annotation: MVSAnnotation, model: Model } | { annotation: undefined, model: undefined } {
    const models = structure.isEmpty ? [] : structure.models;
    for (const model of models) {
        if (model.customProperties.has(MVSAnnotationsProvider.descriptor)) {
            const annots = MVSAnnotationsProvider.get(model).value;
            const annotation = annots?.getAnnotation(annotationId);
            if (annotation) {
                return { annotation, model };
            }
        }
    }
    return { annotation: undefined, model: undefined };
}

type FieldRemapping = Record<string, string | null>;

/** Mapping `ElementIndex` -> annotation row index for all elements in a `Model`.
 * `-1` means no row applies to the element.
 * `null` means no row applies to any element. */
type IndexedModel = number[] | null;

/** Main class for processing MVS annotation */
export class MVSAnnotation {

    /** Number of annotation rows. */
    public nRows: number;

    constructor(
        public data: MVSAnnotationData,
        public schema: MVSAnnotationSchema,
        public fieldRemapping: FieldRemapping,
    ) {
        this.nRows = getRowCount(data);
    }

    /** Create a new `MVSAnnotation` based on specification `spec`. Use `file` if provided, otherwise download the file.
     * Throw error if download fails or problem with data. */
    static async fromSpec(ctx: CustomProperty.Context, spec: MVSAnnotationSpec, file?: MVSAnnotationFile): Promise<MVSAnnotation> {
        file ??= await getFileFromSource(ctx, annotationSourceFromSpec(spec));

        let data: MVSAnnotationData;
        switch (file.format) {
            case 'json':
                data = file;
                break;
            case 'cif':
                if (file.data.blocks.length === 0) throw new Error('No block in CIF');
                const blockSpec = spec.cifBlock;
                let block: CifBlock;
                switch (blockSpec.name) {
                    case 'header':
                        const foundBlock = file.data.blocks.find(b => b.header === blockSpec.params.header);
                        if (!foundBlock) throw new Error(`CIF block with header "${blockSpec.params.header}" not found`);
                        block = foundBlock;
                        break;
                    case 'index':
                        block = file.data.blocks[blockSpec.params.index];
                        if (!block) throw new Error(`CIF block with index ${blockSpec.params.index} not found`);
                        break;
                }
                const categoryName = spec.cifCategory ?? Object.keys(block.categories)[0];
                if (!categoryName) throw new Error('There are no categories in CIF block');
                const category = block.categories[categoryName];
                if (!category) throw new Error(`CIF category "${categoryName}" not found`);
                data = { format: 'cif', data: category };
                break;
        }
        return new MVSAnnotation(data, spec.schema, Object.fromEntries(spec.fieldRemapping.map(e => [e.standardName, e.actualName])));
    }

    static createEmpty(schema: MVSAnnotationSchema): MVSAnnotation {
        return new MVSAnnotation({ format: 'json', data: [] }, schema, {});
    }

    /** Return value of field `fieldName` assigned to location `loc`, if any */
    getValueForLocation(loc: StructureElement.Location, fieldName: string): string | undefined {
        const indexedModel = this.getIndexedModel(loc.unit.model, loc.unit.conformation.operator.instanceId);
        const iRow = (indexedModel !== null) ? indexedModel[loc.element] : -1;
        return this.getValueForRow(iRow, fieldName);
    }
    /** Return value of field `fieldName` assigned to `i`-th annotation row, if any */
    getValueForRow(i: number, fieldName: string): string | undefined {
        if (i < 0) return undefined;
        switch (this.data.format) {
            case 'json':
                const value = getValueFromJson(i, fieldName, this.data.data);
                if (value === undefined || typeof value === 'string') return value;
                else return `${value}`;
            case 'cif':
                return getValueFromCif(i, fieldName, this.data.data);
        }
    }

    /** Return cached `ElementIndex` -> `MVSAnnotationRow` mapping for `Model` (or create it if not cached yet) */
    private getIndexedModel(model: Model, instanceId: string): IndexedModel {
        const key = this.hasInstanceIds() ? `${model.id}:${instanceId}` : model.id;
        if (!this._indexedModels.has(key)) {
            const result = this.getRowForEachAtom(model, instanceId);
            this._indexedModels.set(key, result);
        }
        return this._indexedModels.get(key)!;
    }
    /** Cached `IndexedModel` per `Model.id` (if annotation contains no instanceIds)
     * or per `Model.id:instanceId` combination (if at least one row contains instanceId). */
    private _indexedModels = new Map<string, IndexedModel>();

    /** Create `ElementIndex` -> `MVSAnnotationRow` mapping for `Model` */
    private getRowForEachAtom(model: Model, instanceId: string): IndexedModel {
        const indices = IndicesAndSortings.get(model);
        const nAtoms = model.atomicHierarchy.atoms._rowCount;
        let result: IndexedModel = null;
        const rows = this.getRows();
        for (let i = 0, nRows = rows.length; i < nRows; i++) {
            const row = rows[i];
            const atomRanges = getAtomRangesForRow(row, model, instanceId, indices);
            if (AtomRanges.count(atomRanges) === 0) continue;
            result ??= Array(nAtoms).fill(-1);
            AtomRanges.foreach(atomRanges, (from, to) => result!.fill(i, from, to));
        }
        return result;
    }

    /** Parse and return all annotation rows in this annotation, or return cached result if available */
    getRows(): readonly MVSAnnotationRow[] {
        return this._rows ??= this._getRows();
    }
    /** Cached annotation rows. Do not use directly, use `getRows` instead. */
    private _rows: MVSAnnotationRow[] | undefined = undefined;
    /** Parse and return all annotation rows in this annotation */
    private _getRows(): MVSAnnotationRow[] {
        switch (this.data.format) {
            case 'json':
                return getRowsFromJson(this.data.data, this.schema, this.fieldRemapping);
            case 'cif':
                return getRowsFromCif(this.data.data, this.schema, this.fieldRemapping);
        }
    }

    /** Return `true` if some rows in the annotation contain `instance_id` field. */
    private hasInstanceIds(): boolean {
        return this._hasInstanceIds ??= this.getRows().some(row => isDefined(row.instance_id));
    }
    private _hasInstanceIds?: boolean = undefined;

    /** Return list of all distinct values appearing in field `fieldName`, in order of first occurrence. Ignores special values `.` and `?`. If `caseInsensitive`, make all values uppercase. */
    getDistinctValuesInField(fieldName: string, caseInsensitive: boolean): string[] {
        const seen = new Set<string | undefined>();
        const out = [];
        for (let i = 0; i < this.nRows; i++) {
            let value = this.getValueForRow(i, fieldName);
            if (caseInsensitive) value = value?.toUpperCase();
            if (value !== undefined && !seen.has(value)) {
                seen.add(value);
                out.push(value);
            }
        }
        return out;
    }
}

function getValueFromJson<T>(rowIndex: number, fieldName: string, data: Jsonable): T | undefined {
    const js = data as any;
    if (Array.isArray(js)) {
        const row = js[rowIndex] ?? {};
        return row[fieldName];
    } else {
        const column = js[fieldName] ?? [];
        return column[rowIndex];
    }
}
function getValueFromCif(rowIndex: number, fieldName: string, data: CifCategory): string | undefined {
    const column = data.getField(fieldName);
    if (!column) return undefined;
    if (column.valueKind(rowIndex) !== Column.ValueKind.Present) return undefined;
    return column.str(rowIndex);
}

/** Return number of rows in this annotation (without parsing all the data) */
function getRowCount(data: MVSAnnotationData): number {
    switch (data.format) {
        case 'json':
            return getRowCountFromJson(data.data);
        case 'cif':
            return getRowCountFromCif(data.data);
    }
}
function getRowCountFromJson(data: Jsonable): number {
    const js = data as any;
    if (Array.isArray(js)) {
        // array of objects
        return js.length;
    } else {
        // object of arrays
        const keys = Object.keys(js);
        if (keys.length > 0) {
            return js[keys[0]].length;
        } else {
            return 0;
        }
    }
}
function getRowCountFromCif(data: CifCategory): number {
    return data.rowCount;
}

function getRowsFromJson(data: Jsonable, schema: MVSAnnotationSchema, fieldRemapping: FieldRemapping): MVSAnnotationRow[] {
    const js = data as any;
    const cifSchema = getCifAnnotationSchema(schema);
    const cifSchemaKeys = Object.keys(cifSchema);
    if (Array.isArray(js)) {
        // array of objects
        return js.map(row => pickObjectKeysWithRemapping(row, cifSchemaKeys, fieldRemapping));
    } else {
        // object of arrays
        const selectedFields: Record<string, any[]> = pickObjectKeysWithRemapping(js, cifSchemaKeys, fieldRemapping);
        return objectOfArraysToArrayOfObjects(selectedFields);
    }
}

function getRowsFromCif(data: CifCategory, schema: MVSAnnotationSchema, fieldRemapping: FieldRemapping): MVSAnnotationRow[] {
    const cifSchema = getCifAnnotationSchema(schema);
    const cifSchemaKeys = Object.keys(cifSchema) as (keyof typeof cifSchema)[];
    const columns: Partial<Record<keyof typeof cifSchema, any[]>> = {};
    for (const key of cifSchemaKeys) {
        let srcKey = fieldRemapping[key];
        if (srcKey === null) continue; // Ignore key
        if (srcKey === undefined) srcKey = key; // Implicit key mapping
        const columnArray = getArrayFromCifCategory(data, srcKey, cifSchema[key]); // Avoiding `column.toArray` as it replaces . and ? fields by 0 or ''
        if (columnArray) columns[key] = columnArray;
    }
    return objectOfArraysToArrayOfObjects(columns);
}

/** Load data from a specific column in a CIF category into an array. Load `.` and `?` as undefined. */
function getArrayFromCifCategory<T>(data: CifCategory, columnName: string, columnSchema: Column.Schema): (T | undefined)[] | undefined {
    if (data.getField(columnName) === undefined) return undefined;

    const table = toTable({ [columnName]: columnSchema }, data); // a bit dumb, I don't know how to make column directly
    const column = table[columnName];
    return getArrayFromCifColumn(column); // Avoiding `column.toArray` as it replaces . and ? fields by 0 or ''
}

/** Same as `column.toArray` but reads `.` and `?` as undefined (instead of using type defaults) */
function getArrayFromCifColumn<T>(column: Column<T>): (T | undefined)[] {
    const nRows = column.rowCount;
    const Present = Column.ValueKind.Present;
    const out: (T | undefined)[] = new Array(nRows);
    for (let iRow = 0; iRow < nRows; iRow++) {
        out[iRow] = column.valueKind(iRow) === Present ? column.value(iRow) : undefined;
    }
    return out;
}

async function getFileFromSource(ctx: CustomProperty.Context, source: MVSAnnotationSource, model?: Model): Promise<MVSAnnotationFile> {
    switch (source.kind) {
        case 'source-cif':
            return { format: 'cif', data: getSourceFileFromModel(model) };
        case 'url':
            const url = Asset.getUrlAsset(ctx.assetManager, source.url);
            const dataType = MVSAnnotationFormatTypes[source.format];
            const dataWrapper = await ctx.assetManager.resolve(url, dataType).runInContext(ctx.runtime);
            const rawData = dataWrapper.data;
            if (!rawData) throw new Error('Missing data');
            switch (source.format) {
                case 'json':
                    const json = JSON.parse(rawData as string) as Jsonable;
                    return { format: 'json', data: json };
                case 'cif':
                case 'bcif':
                    const parsed = await CIF.parse(rawData).run();
                    if (parsed.isError) throw new Error(`Failed to parse ${source.format}`);
                    return { format: 'cif', data: parsed.result };
            }
    }
}

/** Like `sources.map(s => safePromise(getFileFromSource(ctx, s)))`
 * but downloads a repeating source only once. */
async function getFilesFromSources(ctx: CustomProperty.Context, sources: MVSAnnotationSource[], model?: Model): Promise<Maybe<MVSAnnotationFile>[]> {
    const promises: { [key: string]: Promise<Maybe<MVSAnnotationFile>> } = {};
    for (const src of sources) {
        const key = canonicalJsonString(src);
        promises[key] ??= safePromise(getFileFromSource(ctx, src, model));
    }
    const files = await promiseAllObj(promises);
    return sources.map(src => files[canonicalJsonString(src)]);
}

function getSourceFileFromModel(model?: Model): CifFile {
    if (model && MmcifFormat.is(model.sourceData)) {
        if (model.sourceData.data.file) {
            return model.sourceData.data.file;
        } else {
            const frame = model.sourceData.data.frame;
            const block = CifBlock(Array.from(frame.categoryNames), frame.categories, frame.header);
            const file = CifFile([block]);
            return file;
        }
    } else {
        console.warn('Could not get CifFile from Model, returning empty CifFile');
        return CifFile([]);
    }
}

function annotationSourceFromSpec(s: MVSAnnotationSpec): MVSAnnotationSource {
    switch (s.source.name) {
        case 'url':
            return { kind: 'url', ...s.source.params };
        case 'source-cif':
            return { kind: 'source-cif' };
    }
}

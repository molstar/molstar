import { CIF, CifBlock } from '../../mol-io/reader/cif';
import { Box3D } from '../../mol-math/geometry';
import { Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { volumeFromDensityServerData } from '../../mol-model-formats/volume/density-server';
import { CustomProperties } from '../../mol-model/custom-property';
import { Grid, Volume } from '../../mol-model/volume';
import { Segment } from './cellstar-api/data';
import { lazyGetter } from './helpers';


export class LatticeSegmentation {
    private segments: number[];
    private sets: number[];
    /** Maps setId to a set of segmentIds*/
    private segmentMap: Map<number, Set<number>>; // computations with objects might be actually faster than with Maps and Sets?
    /** Maps segmentId to a set of setIds*/
    private inverseSegmentMap: Map<number, Set<number>>;
    private grid: Grid;

    private constructor(segmentationDataBlock: CifBlock, grid: Grid) {
        const segmentationValues = segmentationDataBlock!.categories['segmentation_data_3d'].getField('values')?.toIntArray()!;
        this.segmentMap = LatticeSegmentation.makeSegmentMap(segmentationDataBlock);
        this.inverseSegmentMap = LatticeSegmentation.invertMultimap(this.segmentMap);
        this.segments = Array.from(this.inverseSegmentMap.keys());
        this.sets = Array.from(this.segmentMap.keys());
        this.grid = grid;
        this.grid.cells.data = Tensor.Data1(segmentationValues);
    }

    public static async fromCifBlock(segmentationDataBlock: CifBlock) {
        const densityServerCif = CIF.schema.densityServer(segmentationDataBlock);
        // console.log('dscif', densityServerCif);
        const volume = await volumeFromDensityServerData(densityServerCif).run();
        // console.log('volume', volume);
        const grid = volume.grid;
        return new LatticeSegmentation(segmentationDataBlock, grid);
    }

    public createSegment_old(segId: number): Volume {
        // console.time('createSegment_old');
        const n = this.grid.cells.data.length;
        const newData = new Float32Array(n);

        for (let i = 0; i < n; i++) {
            newData[i] = this.segmentMap.get(this.grid.cells.data[i])?.has(segId) ? 1 : 0;
        }

        const result: Volume = {
            sourceData: { kind: 'custom', name: 'test', data: newData as any },
            customProperties: new CustomProperties(),
            _propertyData: {},
            grid: {
                ...this.grid,
                // stats: { min: 0, max: 1, mean: newMean, sigma: arrayRms(newData) },
                stats: { min: 0, max: 1, mean: 0, sigma: 1 },
                cells: {
                    ...this.grid.cells,
                    data: newData as any,
                }
            }
        };
        // console.timeEnd('createSegment_old');
        return result;
    }

    public hasSegment(segId: number): boolean {
        return this.inverseSegmentMap.has(segId);
    }
    public createSegment(seg: Segment, propertyData?: {[key: string]: any}): Volume {
        const { space, data }: Tensor = this.grid.cells;
        const [nx, ny, nz] = space.dimensions;
        const axisOrder = [...space.axisOrderSlowToFast];
        const get = space.get;
        const cell = Box.create(0, nx, 0, ny, 0, nz);

        const EXPAND_START = 2; // We need to add 2 layers of zeros, probably because of a bug in GPU marching cubes implementation
        const EXPAND_END = 1;
        let bbox = this.getSegmentBoundingBoxes()[seg.id];
        bbox = Box.expand(bbox, EXPAND_START, EXPAND_END);
        bbox = Box.confine(bbox, cell);
        const [ox, oy, oz] = Box.origin(bbox);
        const [mx, my, mz] = Box.size(bbox);
        // n, i refer to original box; m, j to the new box

        const newSpace = Tensor.Space([mx, my, mz], axisOrder, Uint8Array);
        const newTensor = Tensor.create(newSpace, newSpace.create());
        const newData = newTensor.data;
        const newSet = newSpace.set;

        const sets = this.inverseSegmentMap.get(seg.id);
        if (!sets) throw new Error(`This LatticeSegmentation does not contain segment ${seg.id}`);

        for (let jz = 0; jz < mz; jz++) {
            for (let jy = 0; jy < my; jy++) {
                for (let jx = 0; jx < mx; jx++) {
                    // Iterating in ZYX order is faster (probably fewer cache misses)
                    const setId = get(data, ox + jx, oy + jy, oz + jz);
                    const value = sets.has(setId) ? 1 : 0;
                    newSet(newData, jx, jy, jz, value);
                }
            }
        }

        const transform = this.grid.transform;
        let newTransform: Grid.Transform;
        if (transform.kind === 'matrix') {
            throw new Error('Not implemented for transform of kind "matrix"'); // TODO ask if this is really needed
        } else if (transform.kind === 'spacegroup') {
            const newFractionalBox = Box.toFractional(bbox, cell);
            const origFractSize = Vec3.sub(Vec3.zero(), transform.fractionalBox.max, transform.fractionalBox.min);
            Vec3.mul(newFractionalBox.min, newFractionalBox.min, origFractSize);
            Vec3.mul(newFractionalBox.max, newFractionalBox.max, origFractSize);
            Vec3.add(newFractionalBox.min, newFractionalBox.min, transform.fractionalBox.min);
            Vec3.add(newFractionalBox.max, newFractionalBox.max, transform.fractionalBox.min);
            newTransform = { ...transform, fractionalBox: newFractionalBox };
        } else {
            throw new Error(`Unknown transform kind: ${transform}`);
        }
        const result = {
            sourceData: { kind: 'custom', name: 'test', data: newTensor.data as any },
            label: seg.biological_annotation.name ?? `Segment ${seg.id}`,
            customProperties: new CustomProperties(),
            _propertyData: propertyData ?? {},
            grid: {
                stats: { min: 0, max: 1, mean: 0, sigma: 1 },
                cells: newTensor,
                transform: newTransform,
            }
        } as Volume;
        return result;
    }

    private static _getSegmentBoundingBoxes(self: LatticeSegmentation) {
        const { space, data }: Tensor = self.grid.cells;
        const [nx, ny, nz] = space.dimensions;
        const get = space.get;

        const setBoxes: { [setId: number]: Box } = {}; // with object this is faster than with Map
        self.sets.forEach(setId => setBoxes[setId] = Box.create(nx, -1, ny, -1, nz, -1));

        for (let iz = 0; iz < nz; iz++) {
            for (let iy = 0; iy < ny; iy++) {
                for (let ix = 0; ix < nx; ix++) {
                    // Iterating in ZYX order is faster (probably fewer cache misses)
                    const setId = get(data, ix, iy, iz);
                    Box.addPoint_InclusiveEnd(setBoxes[setId], ix, iy, iz);
                }
            }
        }

        const segmentBoxes: { [segmentId: number]: Box } = {};
        self.segments.forEach(segmentId => segmentBoxes[segmentId] = Box.create(nx, -1, ny, -1, nz, -1));
        self.inverseSegmentMap.forEach((setIds, segmentId) => {
            setIds.forEach(setId => {
                segmentBoxes[segmentId] = Box.cover(segmentBoxes[segmentId], setBoxes[setId]);
            });
        });

        for (const segmentId in segmentBoxes) {
            if (segmentBoxes[segmentId][5] === -1) { // segment's box left unchanged -> contains no voxels
                segmentBoxes[segmentId] = Box.create(0, 1, 0, 1, 0, 1);
            } else {
                segmentBoxes[segmentId] = Box.expand(segmentBoxes[segmentId], 0, 1); // inclusive end -> exclusive end
            }
        }
        return segmentBoxes;
    }
    private getSegmentBoundingBoxes = lazyGetter(() => LatticeSegmentation._getSegmentBoundingBoxes(this));

    private static invertMultimap<K, V>(map: Map<K, Set<V>>): Map<V, Set<K>> {
        const inverted = new Map<V, Set<K>>();
        map.forEach((values, key) => {
            values.forEach(value => {
                if (!inverted.has(value)) inverted.set(value, new Set<K>());
                inverted.get(value)?.add(key);
            });
        });
        return inverted;
    }

    private static makeSegmentMap(segmentationDataBlock: CifBlock): Map<number, Set<number>> {
        const setId = segmentationDataBlock.categories['segmentation_data_table'].getField('set_id')?.toIntArray()!;
        const segmentId = segmentationDataBlock.categories['segmentation_data_table'].getField('segment_id')?.toIntArray()!;
        const map = new Map<number, Set<number>>();
        for (let i = 0; i < segmentId.length; i++) {
            if (!map.has(setId[i])) {
                map.set(setId[i], new Set());
            }
            map.get(setId[i])!.add(segmentId[i]);
        }
        return map;
    }

    public benchmark(segment: Segment) {
        const N = 100;

        console.time(`createSegment ${segment.id} ${N}x`);
        for (let i = 0; i < N; i++) {
            this.getSegmentBoundingBoxes = lazyGetter(() => LatticeSegmentation._getSegmentBoundingBoxes(this));
            this.createSegment(segment);
        }
        console.timeEnd(`createSegment ${segment.id} ${N}x`);
    }
}


type Box = [number, number, number, number, number, number];

/** Represents a 3D box in integer coordinates. xFrom... is inclusive, xTo... is exclusive. */
namespace Box {
    export function create(xFrom: number, xTo: number, yFrom: number, yTo: number, zFrom: number, zTo: number): Box {
        return [xFrom, xTo, yFrom, yTo, zFrom, zTo];
    }
    export function expand(box: Box, expandFrom: number, expandTo: number): Box {
        const [xFrom, xTo, yFrom, yTo, zFrom, zTo] = box;
        return [xFrom - expandFrom, xTo + expandTo, yFrom - expandFrom, yTo + expandTo, zFrom - expandFrom, zTo + expandTo];
    }
    export function confine(box1: Box, box2: Box): Box {
        const [xFrom1, xTo1, yFrom1, yTo1, zFrom1, zTo1] = box1;
        const [xFrom2, xTo2, yFrom2, yTo2, zFrom2, zTo2] = box2;
        return [
            Math.max(xFrom1, xFrom2), Math.min(xTo1, xTo2),
            Math.max(yFrom1, yFrom2), Math.min(yTo1, yTo2),
            Math.max(zFrom1, zFrom2), Math.min(zTo1, zTo2)
        ];
    }
    export function cover(box1: Box, box2: Box): Box {
        const [xFrom1, xTo1, yFrom1, yTo1, zFrom1, zTo1] = box1;
        const [xFrom2, xTo2, yFrom2, yTo2, zFrom2, zTo2] = box2;
        return [
            Math.min(xFrom1, xFrom2), Math.max(xTo1, xTo2),
            Math.min(yFrom1, yFrom2), Math.max(yTo1, yTo2),
            Math.min(zFrom1, zFrom2), Math.max(zTo1, zTo2)
        ];
    }
    export function size(box: Box): [number, number, number] {
        const [xFrom, xTo, yFrom, yTo, zFrom, zTo] = box;
        return [xTo - xFrom, yTo - yFrom, zTo - zFrom];
    }
    export function origin(box: Box): [number, number, number] {
        const xFrom = box[0];
        const yFrom = box[2];
        const zFrom = box[4];
        return [xFrom, yFrom, zFrom];
    }
    export function log(name: string, box: Box): void {
        const [xFrom, xTo, yFrom, yTo, zFrom, zTo] = box;
        console.log(`Box ${name}: [${xFrom}:${xTo}, ${yFrom}:${yTo}, ${zFrom}:${zTo}], size: ${size(box)}`);
    }
    export function toFractional(box: Box, relativeTo: Box): Box3D {
        const [xFrom, xTo, yFrom, yTo, zFrom, zTo] = box;
        const [x0, y0, z0] = origin(relativeTo);
        const [sizeX, sizeY, sizeZ] = size(relativeTo);
        const min = Vec3.create((xFrom - x0) / sizeX, (yFrom - y0) / sizeY, (zFrom - z0) / sizeZ);
        const max = Vec3.create((xTo - x0) / sizeX, (yTo - y0) / sizeY, (zTo - z0) / sizeZ);
        return Box3D.create(min, max);
    }
    export function addPoint_InclusiveEnd(box: Box, x: number, y: number, z: number): void {
        if (x < box[0]) box[0] = x;
        if (x > box[1]) box[1] = x;
        if (y < box[2]) box[2] = y;
        if (y > box[3]) box[3] = y;
        if (z < box[4]) box[4] = z;
        if (z > box[5]) box[5] = z;
    }
    export function equal(box1: Box, box2: Box): boolean {
        return box1.every((value, i) => value === box2[i]);
    }
}


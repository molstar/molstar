/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { StateTransforms } from '../../mol-plugin-state/transforms';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { VolsegEntryData } from './entry-root';
import { CreateShapePrimitivesProvider, ShapePrimitivesData } from './shape_primitives';


const GEOMETRIC_SEGMENTATION_GROUP_TAG = 'geometric-segmentation-group';

export class VolsegGeometricSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadGeometricSegmentation() {
        const hasGeometricSegmentation = this.entryData.metadata.raw.grid.geometric_segmentation?.exists;
        if (hasGeometricSegmentation) {
            const url = this.entryData.api.shapePrimitivesUrl(this.entryData.source, this.entryData.entryId);

            let group = this.entryData.findNodesByTags(GEOMETRIC_SEGMENTATION_GROUP_TAG)[0]?.transform.ref;
            if (!group) {
                const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
                    { label: 'Segmentation', description: 'Geometric segmentation' }, { tags: [GEOMETRIC_SEGMENTATION_GROUP_TAG], state: { isCollapsed: true } }).commit();
                group = newGroupNode.ref;
            }

            const primitivesData = await this.entryData._resolveStringUrl(url);

            const parsedData: ShapePrimitivesData = JSON.parse(primitivesData);
            console.log('parsedData', parsedData);
            const geometricSegmentationNode = await this.entryData.newUpdate().to(group)
                .apply(CreateShapePrimitivesProvider, { data: parsedData })
                .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha: 0.5 })
                .commit();
        }
    }
}
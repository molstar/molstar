import { StructureElement, StructureProperties } from '../../../mol-model/structure';
import { LociLabel } from '../../../mol-plugin-state/manager/loci-label';
import { SbNcbrPartialChargesPropertyProvider, hasPartialChargesCategories } from './property';
import { Loci } from '../../../mol-model/loci';
import { PluginContext } from '../../../mol-plugin/context';
import { LociLabelProvider } from '../../../mol-plugin-state/manager/loci-label';

export function SbNcbrPartialChargesLociLabelProvider(ctx: PluginContext): LociLabelProvider {
    return {
        label: (loci: Loci) => {
            if (!StructureElement.Loci.is(loci)) return;

            const model = loci.structure.model;
            if (!hasPartialChargesCategories(model)) return;
            const data = SbNcbrPartialChargesPropertyProvider.get(model).value;
            if (!data) return;

            const loc = StructureElement.Loci.getFirstLocation(loci);
            if (!loc) return;

            const granularity = ctx.managers.interactivity.props.granularity;
            if (granularity !== 'element' && granularity !== 'residue') {
                return;
            }

            const atomId = StructureProperties.atom.id(loc);
            const { typeIdToAtomIdToCharge, typeIdToResidueToCharge } = data;

            const typeId = SbNcbrPartialChargesPropertyProvider.props(model).typeId;
            const showResidueCharge = granularity === 'residue';
            const charge = showResidueCharge
                ? typeIdToResidueToCharge.get(typeId)?.get(atomId)
                : typeIdToAtomIdToCharge.get(typeId)?.get(atomId);
            const label = granularity === 'residue' ? 'Residue charge' : 'Atom charge';

            return `<strong>${label}: ${charge?.toFixed(4) || 'undefined'}</strong>`;
        },
        group: (label: LociLabel): string => (label as string).toString().replace(/Model [0-9]+/g, 'Models'),
    };
}

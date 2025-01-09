import { Bond, StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Location } from '../../../mol-model/location';
import { SbNcbrPartialChargesPropertyProvider } from './property';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';

const Colors = {
    Bond: Color(0xffffff),
    Error: Color(0x00ff00),
    MissingCharge: Color(0x66ff00),

    Negative: Color(0xff0000),
    Zero: Color(0xffffff),
    Positive: Color(0x0000ff),

    getColor: (charge: number, maxCharge: number): Color => {
        if (charge === 0) return Colors.Zero;
        if (charge <= -maxCharge) return Colors.Negative;
        if (charge >= maxCharge) return Colors.Positive;

        const t = maxCharge !== 0 ? Math.abs(charge) / maxCharge : 1;
        const endColor = charge < 0 ? Colors.Negative : Colors.Positive;

        return Color.interpolate(Colors.Zero, endColor, t);
    },
};

export const PartialChargesThemeParams = {
    maxAbsoluteCharge: PD.Numeric(
        0,
        { min: 0 },
        {
            label: 'Charge Range',
        }
    ),
    absolute: PD.Boolean(false, { isHidden: false, label: 'Use Range' }),
    chargeType: PD.Select(
        'residue',
        [
            ['atom', 'Atom charges'],
            ['residue', 'Residue charges'],
        ],
        { isHidden: false }
    ),
};
export type PartialChargesThemeParams = typeof PartialChargesThemeParams;

export function getPartialChargesThemeParams() {
    return PD.clone(PartialChargesThemeParams);
}

export function PartialChargesColorTheme(
    ctx: ThemeDataContext,
    props: PD.Values<PartialChargesThemeParams>
): ColorTheme<PartialChargesThemeParams> {
    const model = ctx.structure?.models[0];
    if (!model) {
        throw new Error('No model found');
    }
    const data = SbNcbrPartialChargesPropertyProvider.get(model).value;
    if (!data) {
        throw new Error('No partial charges data found');
    }

    const { absolute, chargeType } = props;
    const { typeIdToAtomIdToCharge, typeIdToResidueToCharge, maxAbsoluteAtomCharges, maxAbsoluteResidueCharges } = data;
    const typeId = SbNcbrPartialChargesPropertyProvider.props(model).typeId;
    const atomToCharge = typeIdToAtomIdToCharge.get(typeId);
    const residueToCharge = typeIdToResidueToCharge.get(typeId);

    let maxCharge = 0;
    if (absolute) {
        maxCharge = props.maxAbsoluteCharge < 0 ? 0 : props.maxAbsoluteCharge;
    } else if (chargeType === 'atom') {
        maxCharge = maxAbsoluteAtomCharges.get(typeId) || 0;
    } else {
        maxCharge = maxAbsoluteResidueCharges.get(typeId) || 0;
    }

    // forces coloring updates
    const contextHash = SbNcbrPartialChargesPropertyProvider.get(model)?.version;

    const chargeMap = chargeType === 'atom' ? atomToCharge : residueToCharge;

    let color: LocationColor;

    if (!chargeMap) {
        color = (_: Location): Color => Colors.MissingCharge;
    } else {
        color = (location: Location): Color => {
            let id = -1;
            if (StructureElement.Location.is(location)) {
                if (Unit.isAtomic(location.unit)) {
                    id = StructureProperties.atom.id(location);
                }
            } else if (Bond.isLocation(location)) {
                if (Unit.isAtomic(location.aUnit)) {
                    const l = StructureElement.Location.create(ctx.structure?.root);
                    l.unit = location.aUnit;
                    l.element = location.aUnit.elements[location.aIndex];
                    id = StructureProperties.atom.id(l);
                }
            }

            const charge = chargeMap.get(id);

            if (charge === undefined) {
                console.warn('No charge found for id', id);
                return Colors.MissingCharge;
            }

            return Colors.getColor(charge, maxCharge);
        };
    }

    return {
        factory: PartialChargesColorTheme,
        granularity: 'group',
        color,
        props,
        description: 'Color atoms and residues based on their partial charge.',
        preferSmoothing: false,
        contextHash,
    };
}

export const SbNcbrPartialChargesColorThemeProvider: ColorTheme.Provider<
PartialChargesThemeParams,
'sb-ncbr-partial-charges'
> = {
    label: 'SB NCBR Partial Charges',
    name: 'sb-ncbr-partial-charges',
    category: ColorTheme.Category.Atom,
    factory: PartialChargesColorTheme,
    getParams: getPartialChargesThemeParams,
    defaultValues: PD.getDefaultValues(PartialChargesThemeParams),
    isApplicable: (ctx: ThemeDataContext) =>
        !!ctx.structure &&
        ctx.structure.models.some((model) => SbNcbrPartialChargesPropertyProvider.isApplicable(model)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) =>
            data.structure
                ? SbNcbrPartialChargesPropertyProvider.attach(ctx, data.structure.models[0], void 0, true)
                : Promise.resolve(),
        detach: (data) => data.structure && SbNcbrPartialChargesPropertyProvider.ref(data.structure.models[0], false),
    },
};

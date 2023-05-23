/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PresetStructureRepresentations } from "../../mol-plugin-state/builder/structure/representation-preset";
import { Color } from "../../mol-util/color";
import { CollapsableControls, PurePluginUIComponent } from "../base";
import { Button } from "../controls/common";
import { MagicWandSvg } from "../controls/icons";
import { ParamDefinition as PD } from "../../mol-util/param-definition";
import { PostprocessingParams } from "../../mol-canvas3d/passes/postprocessing";
import { PluginConfig } from "../../mol-plugin/config";
import { StructureComponentManager } from "../../mol-plugin-state/manager/structure/component";

export class StructureQuickStylesControls extends CollapsableControls {
  defaultState() {
    return {
      isCollapsed: false,
      header: "Quick Styles",
      brand: { accent: "gray" as const, svg: MagicWandSvg },
    };
  }

  renderControls() {
    return (
      <>
        <QuickStyles />
      </>
    );
  }
}

export class QuickStyles extends PurePluginUIComponent {
  state = {
    previousStyle: null,
  };
  leaveTimeout: number | null = null;

  // Define a method to revert to previous style
  async revertStyle() {
    // Code to apply style saved in this.state.previousStyle

    // Check if there's a previous style saved
    if (!this.state.previousStyle) {
      this.default();
    } else {
      // Extract the properties from the previous style
      const { structure, postprocessing } = this.state.previousStyle;

      // Apply the previous structure style
      const { structures } = this.plugin.managers.structure.hierarchy.selection;
      await this.plugin.managers.structure.component.applyPreset(
        structures,
        structure
      );

      // Apply the previous postprocessing style
      if (this.plugin.canvas3d && postprocessing) {
        this.plugin.canvas3d.setProps({ postprocessing });
      }

      // Reset previousStyle state
      this.setState({ previousStyle: null });
    }
  }

  // Add getCurrentStyle() method to return current style
  getCurrentStyle() {
    // Return the current style.
    // const structureState = this.plugin.managers.structure.component.state;
    // Similarly, you could get the current postprocessing state like this:
    // const postprocessingState = this.plugin.canvas3d ? this.plugin.canvas3d.props.postprocessing : null;
    // Now return an object that captures the current state
    // return {
    // structure: structureState,
    // postprocessing: postprocessingState,
    // };
  }
  handleMouseLeave = () => {
    console.log("Mouse Leave");
    // Clear the previous timeout if it exists
    if (this.leaveTimeout !== null) {
      clearTimeout(this.leaveTimeout);
      this.leaveTimeout = null;
    }

    // Set a timeout to delay the onMouseLeave action
    this.leaveTimeout = window.setTimeout(
      () => this.defaultAfterPreview(),
      500
    );

    // Start a timer, giving the onClick event a chance to interrupt if it fires
    // this.leaveTimeout = window.setTimeout(async () => {
    //   // Your existing code here...
    //   // Do not forget to clear the timer at the end
    //   this.defaultAfterPreview();
    //   this.leaveTimeout = null;
    // }, 200);

    // Delay the execution for 200ms
  };

  // Add a method to apply the default style
  async defaultAfterPreview() {
    // Clear the previous timeout if it exists
    if (this.leaveTimeout) {
      clearTimeout(this.leaveTimeout);
    }

    const { structures } = this.plugin.managers.structure.hierarchy.selection;
    const preset =
      this.plugin.config.get(
        PluginConfig.Structure.DefaultRepresentationPreset
      ) || PresetStructureRepresentations.auto.id;
    const provider =
      this.plugin.builders.structure.representation.resolveProvider(preset);
    await this.plugin.managers.structure.component.applyPreset(
      structures,
      provider
    );

    this.plugin.managers.structure.component.setOptions(
      PD.getDefaultValues(StructureComponentManager.OptionsParams)
    );

    if (this.plugin.canvas3d) {
      const p = PD.getDefaultValues(PostprocessingParams);
      this.plugin.canvas3d.setProps({
        postprocessing: { outline: p.outline, occlusion: p.occlusion },
      });
    }
  }

  async default() {
    console.log("Clicked default");
    if (this.leaveTimeout) {
      // If timer exists, clear it to prevent defaultAfterPreview from executing
      clearTimeout(this.leaveTimeout);
      this.leaveTimeout = null;
    }
    const { structures } = this.plugin.managers.structure.hierarchy.selection;
    const preset =
      this.plugin.config.get(
        PluginConfig.Structure.DefaultRepresentationPreset
      ) || PresetStructureRepresentations.auto.id;
    const provider =
      this.plugin.builders.structure.representation.resolveProvider(preset);
    await this.plugin.managers.structure.component.applyPreset(
      structures,
      provider
    );

    this.plugin.managers.structure.component.setOptions(
      PD.getDefaultValues(StructureComponentManager.OptionsParams)
    );

    if (this.plugin.canvas3d) {
      const p = PD.getDefaultValues(PostprocessingParams);
      this.plugin.canvas3d.setProps({
        postprocessing: { outline: p.outline, occlusion: p.occlusion },
      });
    }
  }

  async illustrative() {
    if (this.leaveTimeout) {
      // If timer exists, clear it to prevent defaultAfterPreview from executing
      clearTimeout(this.leaveTimeout);
      this.leaveTimeout = null;
    }
    console.log("Clicked illustrative");

    const { structures } = this.plugin.managers.structure.hierarchy.selection;
    await this.plugin.managers.structure.component.applyPreset(
      structures,
      PresetStructureRepresentations.illustrative
    );

    if (this.plugin.canvas3d) {
      this.plugin.canvas3d.setProps({
        postprocessing: {
          outline: {
            name: "on",
            params: {
              scale: 1,
              color: Color(0x000000),
              threshold: 0.25,
              includeTransparent: true,
            },
          },
          occlusion: {
            name: "on",
            params: {
              multiScale: { name: "off", params: {} },
              radius: 5,
              bias: 0.8,
              blurKernelSize: 15,
              samples: 32,
              resolutionScale: 1,
              color: Color(0x000000),
            },
          },
          shadow: { name: "off", params: {} },
        },
      });
    }
  }

  async stylized() {
    console.log("Clicked stylized");
    if (this.leaveTimeout) {
      // If timer exists, clear it to prevent defaultAfterPreview from executing
      clearTimeout(this.leaveTimeout);
      this.leaveTimeout = null;
    }

    this.plugin.managers.structure.component.setOptions({
      ...this.plugin.managers.structure.component.state.options,
      ignoreLight: true,
    });

    if (this.plugin.canvas3d) {
      const pp = this.plugin.canvas3d.props.postprocessing;
      this.plugin.canvas3d.setProps({
        postprocessing: {
          outline: {
            name: "on",
            params:
              pp.outline.name === "on"
                ? pp.outline.params
                : {
                    scale: 1,
                    color: Color(0x000000),
                    threshold: 0.33,
                    includeTransparent: true,
                  },
          },
          occlusion: {
            name: "on",
            params:
              pp.occlusion.name === "on"
                ? pp.occlusion.params
                : {
                    multiScale: { name: "off", params: {} },
                    radius: 5,
                    bias: 0.8,
                    blurKernelSize: 15,
                    samples: 32,
                    resolutionScale: 1,
                    color: Color(0x000000),
                  },
          },
          shadow: { name: "off", params: {} },
        },
      });
    }
  }

  render() {
    return (
      <div className="msp-flex-row">
        <Button
          noOverflow
          title="Applies default representation preset. Set outline and occlusion effects to defaults."
          onClick={() => this.default()}
          style={{ width: "auto" }}
        >
          Default
        </Button>
        <Button
          noOverflow
          title="Applies no representation preset. Enables outline and occlusion effects. Enables ignore-light representation parameter."
          onClick={() => console.log("Clicked stylized")}
          onMouseEnter={() => this.illustrative().catch(console.error)}
          onMouseLeave={() => this.handleMouseLeave()}
          style={{ width: "auto" }}
        >
          Stylized
        </Button>
        <Button
          noOverflow
          title="Applies illustrative representation preset. Enables outline and occlusion effects. Enables ignore-light parameter."
          onClick={() => this.illustrative().catch(console.error)}
          onMouseEnter={() => this.illustrative().catch(console.error)}
          onMouseLeave={() => this.handleMouseLeave()}
          style={{ width: "auto" }}
        >
          Illustrative
        </Button>
      </div>
    );
  }
}

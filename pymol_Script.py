##################################################                                                                                                          
#pyMoL Figures generator                                                                                                                              
##################################################                                                                                                                                                                                                                                                                                                                   
def pyMoLFig(pdbARM,cavityList):
    from chromophore import chromophoreName
    from lib import chainName
    from initial_Setup import linker_aa_ID, counterion_ID
    from counterions_Placement import TargetInner_pymol, TargetOuter_pymol

    cavity_pymol = '+'.join(cavityList)

    pyMoLScript = pdbARM[:-4]+".pml"
    with open(pyMoLScript, "w") as pyMoLScript:
        pyMoLScript.writelines(["from pymol import cmd,stored \n",
                                "bg_color white \n",
                                "load "+pdbARM+" \n",
                                "set auto_zoom, off \n",
                                "hide everything \n",
                                "hide cartoon, all \n",
                                "hide spheres, all \n",
                                "hide nb_spheres, all \n",
                                "hide sticks, all \n",
                                "#Cavity fpocket \n",
                                "sele resi "+cavity_pymol+" \n",
                                "create Cavity, sele \n",
                                "show lines, Cavity \n",
                                "show surface, Cavity \n",
                                "set transparency, 0.65, Cavity \n",
                                "color red, Cavity \n",
                                "# Chromophore \n",
                                "select resn "+chromophoreName+" \n",
                                "create CHROMOPHORE, sele \n",
                                "show sticks, CHROMOPHORE \n",
                                "color tv_green, CHROMOPHORE \n",
                                "#Linker amino acid \n",
                                "select resi "+linker_aa_ID+" \n",
                                "create linkerAA, sele \n",
                                "show sticks, linkerAA \n",
                                "color blue, linkerAA \n",
                                "#Main counterion \n",
                                "select resi "+counterion_ID+" \n",
                                "create mainCounter, sele \n",
                                "show sticks, mainCounter \n",
                                "color lightblue, mainCounter \n",
                                "hide surface, linkerAA \n",
                                "hide surface, mainCounter \n"
                                "#CL \n",
                                "select resn CL \n",
                                "create CL, sele \n",
                                "show spheres, CL \n",
                                "color red, CL \n",
                                "#NA \n",
                                "select resn NA \n",
                                "create NA, sele \n",
                                "show spheres, NA \n",
                                "color blue, NA \n",
                                "#Target residues \n",
                                "sele resi "+TargetInner_pymol+" \n",
                                "create TargetInner, sele \n",
                                "show lines, TargetInner \n",
                                "color gray, TargetInner \n",
                                "sele resi "+TargetOuter_pymol+" \n",
                                "create TargetOuter, sele \n",
                                "show lines, TargetOuter \n",
                                "color gray, TargetOuter \n",
                                "#Chain \n",
                                "select resn HOH \n",
                                "create WaterHOH, sele \n",
                                "show nb_spheres, WaterHOH \n"
                                "select all and not CHROMOPHORE and not CL and not NA and not WaterHOH and not mainCounter and not linkerAA and not Cavity and not TargetInner and not TargetOuter and not resn HOH \n",
                                "create mainChain, sele \n",
                                "show cartoon, mainChain \n"
                                "color gray, mainChain \n"
                                "set_bond stick_radius, 0.35, all \n",
                                "set sphere_scale, 0.5 \n",
                                "set sphere_quality, 2 \n",
                                "set cartoon_transparency, 0.6 \n",
                                "set ray_trace_mode, 3 \n",
                                "set ray_shadows, 0 \n",
                                "set antialias, 2 \n",
                                "rotate x, -90 \n",
                                "rotate y, -90 \n",
                                "#ray 1200,1200 \n",
                                "#png "+pdbARM[:-4]+".png \n",
                                ])

        print( str("The pyMoL script \x1b[0;33;49m"+str(pdbARM[:-4]+".pml")+"\x1b[0m has been generated. ").rjust(100, '.'))

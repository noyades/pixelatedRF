def createOpenAel(pathName: str,
                  libName: str,
                  gdsName: str,
                  ports: int,
                  portPosition: float,
                  aelName: str,
                  cellName: str):
  """
  This file creates an AEL file that will take a GDS file, input the file into ADS, add ports and prepare the design for simulation.
  The file needs to know the number of ports that will be simulated and the relative position of the ports. These are outputs when 
  the GDS is greated. 
  """
  aelPath = pathName + libName + '_wrk/' + aelName
  with open(aelPath, 'w') as f:
    f.write('/* Begin */\n')
    f.write('de_open_workspace("' + pathName + libName + '_wrk/");\n')
    #f.write('decl macroContext = de_create_new_layout_view("' + libName + '_lib", "cell_1", "layout");\n')
    #f.write('de_show_context_in_new_window(macroContext);\n')
    f.write('de_load_translators_plugin_if_not_loaded();\n')
    f.write('detransdlg_execute_more_options_ok_cb(api_get_current_window(), DE_GDSII_FILE, 1, "' + gdsName + '");\n')
    f.write('de_import_design(DE_GDSII_FILE, FALSE, "' + gdsName + '", "' + libName + '_lib", "", NULL);\n')
    f.write('decl macroContext = de_get_design_context_from_name("' + libName + '_lib:' + cellName + ':layout");\n')
    #f.write('de_close_design("cell_1");\n')
    #f.write('de_delete_cell("' + libName + '_lib", "cell_1");\n') #After import of design cell, delete dummy file that is opened
    f.write('// For an artwork macro: decl macroContext = de_get_current_design_context();\n')
    f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
    f.write('db_set_entry_layerid( de_get_current_design_context(), \
             db_find_layerid_by_name( de_get_current_design_context(), "l_top:drawing" ));\n')
    f.write('de_set_grid_snap_type(19);') 
    if ports == 2:
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[0]) + ',' + str(portPosition[1]) \
              + ',' + str(180) + ', db_layerid(39), 1, "P1", 2);\n') 
      # If all goes well, l_top should come in on layer 39 of ADS metal stack, so put the pins on this layer.
      # It may need to adjust if things change. Working to solve.
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[2]) + ',' + str(portPosition[3]) \
              + ',' + str(90) + ', db_layerid(39), 2, "P2", 2);\n')
    elif ports == 3:
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[0]) + ',' + str(portPosition[1]) \
              + ',' + str(180) + ', db_layerid(39), 1, "P1", 2);\n') 
      # If all goes well, l_top should come in on layer 39 of ADS metal stack, so put the pins on this layer.
      # It may need to adjust if things change. Working to solve.
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[2]) + ',' + str(portPosition[3]) \
              + ',' + str(90) + ', db_layerid(39), 2, "P2", 2);\n')
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[4]) + ',' + str(portPosition[5]) \
              + ',' + str(0) + ', db_layerid(39), 3, "P3", 2);\n')
    elif ports == 4:
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[0]) + ',' + str(portPosition[1]) \
              + ',' + str(180) + ', db_layerid(39), 1, "P1", 2);\n') 
      # If all goes well, l_top should come in on layer 39 of ADS metal stack, so put the pins on this layer.
      # It may need to adjust if things change. Working to solve.
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[2]) + ',' + str(portPosition[3]) \
              + ',' + str(90) + ', db_layerid(39), 2, "P2", 2);\n')
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[4]) + ',' + str(portPosition[5]) \
              + ',' + str(0) + ', db_layerid(39), 3, "P3", 2);\n')
      f.write('decl pin = db_create_pin(macroContext, ' + str(portPosition[6]) + ',' + str(portPosition[7]) \
              + ',' + str(-90) + ', db_layerid(39), 4, "P4", 2);\n')   
    f.write('de_open_substrate_window("' + libName + '_lib", "substrate1");\n')
    f.write('decl macroContext = de_get_design_context_from_name("' + libName + '_lib:' + cellName + ':layout");\n')
    f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
    f.write('de_save_oa_design("' + libName + '_lib:' + cellName + ':layout");\n')
    f.write('decl macroContext = de_get_design_context_from_name("' + libName + '_lib:' + cellName + ':layout");\n')
    f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
    f.write('de_select_all();\n')
    f.write('de_union();\n')
    f.write('decl macroContext = de_get_design_context_from_name("' + libName + '_lib:' + cellName + ':layout");\n')
    f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
    f.write('de_save_oa_design("' + libName + '_lib:' + cellName + ':layout");\n')
    f.write('de_close_all();\n')
    f.write('dex_em_writeSimulationFiles("' + libName + '_lib", "' + cellName + '", "emSetup", "simulation/' \
            + libName + '_lib/' + cellName + '/layout/emSetup_MoM");\n')
    f.write('de_close_workspace_without_prompting();\n')
    f.write('de_exit();\n')

def createCloseAel(pathName: str,
                   libName: str,
                   aelName: str,
                   cellName: str):
  """
  This file creates an AEL file that will clean up the directories after an EM simulation is run and prepare the 
  environment for the next simulation.
  """
  aelPath = pathName + libName + '_wrk/' + aelName
  with open(aelPath, 'w') as f:
    f.write('/* Begin */\n')
    f.write('de_open_workspace("' + pathName + libName + '_wrk/");\n')
    f.write('de_delete_cellview("' + libName + '_lib", "' + cellName + '", "layout");\n')
    f.write('de_close_workspace_without_prompting();\n')
    f.write('de_exit();\n')

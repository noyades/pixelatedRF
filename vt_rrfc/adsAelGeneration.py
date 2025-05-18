class AelTools:
  def __init__(self,
               pathName: str,
               libName: str,
               gdsName: str,
               ports: int,
               portPosition: float,
               aelName: str,
               cellName: str):
    self.pathName = pathName
    self.libName = libName
    self.gdsName = gdsName
    self.ports = ports
    self.portPosition = portPosition
    self.aelName = aelName
    self.cellName = cellName
    
  def create_open_ael(self):
    """
    This file creates an AEL file that will take a GDS file, input the file into ADS, add ports and prepare the design for simulation.
    The file needs to know the number of ports that will be simulated and the relative position of the ports. These are outputs when 
    the GDS is greated. 
    """
    aelPath = self.pathName + self.libName + '_wrk/' + self.aelName
    with open(aelPath, 'w') as f:
      f.write('/* Begin */\n')
      f.write('de_open_workspace("' + self.pathName + self.libName + '_wrk/");\n')
      #f.write('decl macroContext = de_create_new_layout_view("' + libName + '_lib", "cell_1", "layout");\n')
      #f.write('de_show_context_in_new_window(macroContext);\n')
      f.write('de_load_translators_plugin_if_not_loaded();\n')
      f.write('detransdlg_execute_more_options_ok_cb(api_get_current_window(), DE_GDSII_FILE, 1, "' + self.gdsName + '");\n')
      f.write('de_import_design(DE_GDSII_FILE, FALSE, "' + self.gdsName + '", "' + self.libName + '_lib", "", NULL);\n')
      f.write('decl macroContext = de_get_design_context_from_name("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      #f.write('de_close_design("cell_1");\n')
      #f.write('de_delete_cell("' + libName + '_lib", "cell_1");\n') #After import of design cell, delete dummy file that is opened
      f.write('// For an artwork macro: decl macroContext = de_get_current_design_context();\n')
      f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
      f.write('db_set_entry_layerid( de_get_current_design_context(), \
               db_find_layerid_by_name( de_get_current_design_context(), "l_top:drawing" ));\n')
      f.write('de_set_grid_snap_type(19);') 

      for i in range(self.ports):
        angle = [180, 90, 0, -90][i % 4] # Assign angles based on port index
        f.write(f'decl pin = db_create_pin(macroContext, {self.portPosition[2*i]}, {self.portPosition[2*i+1]}, {angle}, db_layerid(39), {i+1}, "P{i+1}", 2);\n')
        
      f.write('de_open_substrate_window("' + self.libName + '_lib", "substrate1");\n')
      f.write('decl macroContext = de_get_design_context_from_name("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
      f.write('de_save_oa_design("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      f.write('decl macroContext = de_get_design_context_from_name("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
      f.write('de_select_all();\n')
      f.write('de_union();\n')
      f.write('decl macroContext = de_get_design_context_from_name("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      f.write('de_bring_context_to_top_or_open_new_window(macroContext);\n')
      f.write('de_save_oa_design("' + self.libName + '_lib:' + self.cellName + ':layout");\n')
      f.write('de_close_all();\n')
      f.write('dex_em_writeSimulationFiles("' + self.libName + '_lib", "' + self.cellName + '", "emSetup", "simulation/' \
              + self.libName + '_lib/' + self.cellName + '/layout/emSetup_MoM");\n')
      f.write('de_close_workspace_without_prompting();\n')
      f.write('de_exit();\n')

  def create_close_ael(self):
    """
    This file creates an AEL file that will clean up the directories after an EM simulation is run and prepare the 
    environment for the next simulation.
    """
    aelPath = self.pathName + self.libName + '_wrk/' + self.aelName
    with open(aelPath, 'w') as f:
      f.write('/* Begin */\n')
      f.write('de_open_workspace("' + self.pathName + self.libName + '_wrk/");\n')
      f.write('de_delete_cellview("' + self.libName + '_lib", "' + self.cellName + '", "layout");\n')
      f.write('de_close_workspace_without_prompting();\n')
      f.write('de_exit();\n')
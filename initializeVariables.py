import tkinter as tk
from tkinter import filedialog
import tkinter.ttk as ttk
import xml.etree.ElementTree as ET
import platform
import os

os_type = platform.system()
# Initialization variables
init_vars = {
    "os_type": os_type,
    "fc": 8.0e9,
    "z0": 50,
    "EL": 90,
    "t": 1.4,
    "cond": 5.88e7,
    "h": 30,
    "t_air": 2*(1.4+30),
    "er": 3.66,
    "simulator": "ADS",
    "simulatorPath": r"C:\Program Files\Keysight\ADS2025_Update1\bin",
    "libName": "dbsEmSim",
    "sim": True,
    "view": False,
    "write": True,
    "ports": 2,
    "sides": 2,
    "shape": 1, 
    "corner": "normal",
    "scale": 1,
    "pixelSize": 8,
    "minPixel": 6,
    "layoutRes": 1,
    "layoutUnit": 25.4e-6/1,
    "pathName": r"H:/My Drive/dbs_optimization",
    "min_freq": 0.0,
    "max_freq": 20.0e9
}

# Descriptions for tooltips
var_descriptions = {
    "fc": "Operating center frequency for electrical length calculations (Hz)",
    "z0": "Desired characteristic impedance of launches (Ohms)",
    "EL": "Desired unit of electrical length in degrees",
    "t": "Thickness of conductor in mils",
    "cond": "Conductivity of the conductors (S/m)",
    "h": "Height of conductor above substrate (mils)",
    "t_air": "Thickness of the air above the conductor layer (mils)",
    "er": "Relative permittivity of the substrate material",
    "simulator": "Simulation tool to use ('ADS' or 'EMX')",
    "simulatorPath": "Path to the simulator executable (if not in PATH)",
    "libName": "Name of the ADS workspace",
    "sim": "Whether to run a simulation (True/False)",
    "view": "View the GDS after each creation (True/False)",
    "write": "Write output files (True/False)",
    "ports": "Number of ports (2, 3, or 4)",
    "sides": "Number of sides for ports (2, 3, or 4)",
    "shape": "Shape of the pixel map (1 for square, 2 for octagon, 3 for ...)",
    "corner": "Corner type: 'overlap' or 'normal'",
    "scale": "Scaling factor for layout",
    "pixelSize": "Size of the randomized pixel in mils",
    "minPixel": "Minimum possible pixel size in mils",
    "layoutRes": "Layout resolution (sub-pixel grid factor)",
    "layoutUnit": "Layout unit in meters",
    "pathName": "Base path for file creation",
    "min_freq": "Minimum frequency for simulation (Hz)",
    "max_freq": "Maximum frequency for simulation (Hz)"
}

# Add this near your var_descriptions
var_units = {
    "fc": "Hz",
    "z0": "Ω",
    "EL": "deg",
    "t": "mil",
    "cond": "S/m",
    "h": "mil",
    "t_air": "mil",
    "er": "",
    "simulator": "",
    "simulatorPath": "",
    "libName": "",
    "sim": "",
    "view": "",
    "write": "",
    "ports": "",
    "sides": "",
    "shape": "",
    "corner": "",
    "scale": "",
    "pixelSize": "mil",
    "minPixel": "mil",
    "layoutRes": "",
    "layoutUnit": "m",
    "pathName": "",
    "min_freq": "Hz",
    "max_freq": "Hz"
}

# Tooltip class
class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        widget.bind("<Enter>", self.show_tip)
        widget.bind("<Leave>", self.hide_tip)

    def show_tip(self, event=None):
        if self.tipwindow or not self.text:
            return
        x, y, _, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 30
        y = y + cy + self.widget.winfo_rooty() + 10
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(tw, text=self.text, justify='left',
                         background="#ffffe0", relief='solid', borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hide_tip(self, event=None):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def ensure_output_folders(pathName):
    """
    Ensures that pathName and its subfolders exist (cross-platform).
    """
    subfolders = ["pixelMaps", "spFiles", "gdsFiles"]
    # Expand user (~) and normalize path for cross-platform compatibility
    base_path = os.path.expanduser(pathName)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    for sub in subfolders:
        subfolder_path = os.path.join(base_path, sub)
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)

def browse_path(entry):
    path = filedialog.askdirectory(title="Select Directory")
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)
        
def save_to_xml(vars_dict, filename):
    root = ET.Element("InitializationVariables")
    for k, v in vars_dict.items():
        child = ET.SubElement(root, k)
        child.text = str(v)
    tree = ET.ElementTree(root)
    tree.write(filename, encoding="utf-8", xml_declaration=True)

def load_from_xml():
    file_path = filedialog.askopenfilename(defaultextension=".xml", filetypes=[("XML files", "*.xml")])
    if not file_path:
        return
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        for child in root:
            key = child.tag
            val = child.text
            if key in entries:
                #if key == "simulator":
                #    entries[key].set(val)
                #else:
                #    entries[key].delete(0, tk.END)
                #    entries[key].insert(0, val)
                # Update init_vars as well
                entry = entries[key]
                if isinstance(entry, tk.StringVar):
                    entry.set(val)
                elif isinstance(entry, tuple):
                    var, custom_entry = entry
                    if val in ["0", "1", "2", "3", "4"]:
                        var.set(val)
                    else:
                        var.set("Other")
                        custom_entry.delete(0, tk.END)
                        custom_entry.insert(0, val)
                else:
                    entry.delete(0, tk.END)
                    entry.insert(0, val)

                if val is not None:
                    if val.lower() == "true":
                        init_vars[key] = True
                    elif val.lower() == "false":
                        init_vars[key] = False
                    else:
                        try:
                            if '.' in val or 'e' in val.lower():
                                init_vars[key] = float(val)
                            else:
                                init_vars[key] = int(val)
                        except ValueError:
                            init_vars[key] = val
        status_label.config(text=f"Loaded from {file_path}")
    except Exception as e:
        status_label.config(text=f"Failed to load: {e}")

def save_vars():
    for key, entry in entries.items():
        if key == "simulator":
            val = entry.get()
        elif key in ["ports", "sides"]:
            var, custom_entry = entry
            val = var.get()
            if val == "Other":
                val = custom_entry.get()
            try:
                init_vars[key] = int(val)
            except ValueError:
                init_vars[key] = val
            continue  # Skip the rest of the loop for these keys
        elif key in ["sim", "view", "write"]:
            val = entry.get()
            if val.lower() == "true":
                init_vars[key] = True
            elif val.lower() == "false":
                init_vars[key] = False
            else:
                init_vars[key] = val
            continue
        else:
            val = entry.get()
            if val.lower() == "true":
                init_vars[key] = True
            elif val.lower() == "false":
                init_vars[key] = False
            else:
                try:
                    if '.' in val or 'e' in val.lower():
                        init_vars[key] = float(val)
                    else:
                        init_vars[key] = int(val)
                except ValueError:
                    init_vars[key] = val
    file_path = filedialog.asksaveasfilename(defaultextension=".xml", filetypes=[("XML files", "*.xml")])
    if file_path:
        save_to_xml(init_vars, file_path)
        status_label.config(text=f"Saved to {file_path}")

def show_custom_entry(var, entry_widget, option_menu_widget):
    if var.get() == "Other":
        entry_widget.grid()
    else:
        entry_widget.grid_remove()

def edit_variables(parent):
    global frame, entries, status_label, libname_label, libname_entry

    window = tk.Toplevel(parent)
    window.title("Edit and Save Initialization Variables")
    frame = tk.Frame(window, padx=10, pady=10)
    frame.pack()

    tk.Label(frame, text="Edit variables and save to XML").pack(pady=5)


    # Define groups: group name -> list of variable keys
    var_groups = {
        "Physical Parameters": ["fc", "z0", "EL", "t", "cond", "h", "t_air", "er"],
        "Layout Parameters": ["scale", "pixelSize", "minPixel", "shape", "ports", "sides", "corner", "layoutRes", "layoutUnit", "pathName"],
        "Simulation Settings": ["simulator", "simulatorPath", "libName", "sim", "view", "write", "min_freq", "max_freq"],
    }

    simulator_options = ["ADS", "EMX", "MEEP", "Other"]
    corner_options = ["normal", "noverlap", "overlap"]
    ports_sides_options = ["0", "1", "2", "3", "4", "Other"]

    entries = {}
    
    def update_libname_visibility(*args):
        if entries["simulator"].get() == "ADS":
            libname_label.grid()
            libname_entry.grid()
        else:
            libname_label.grid_remove()
            libname_entry.grid_remove()

    entries = {}
    #vars_frame = tk.Frame(frame)
    #vars_frame.pack(pady=5)
    #row = 0
    # Create notebook (tabs)
    notebook = ttk.Notebook(frame)
    notebook.pack(pady=5, fill="both", expand=True)

    tab_frames = {}

    for group_name, keys in var_groups.items():
        #header = tk.Label(vars_frame, text=group_name, font=("Arial", 10, "bold"), anchor="w")
        #header.grid(row=row, column=0, columnspan=2, sticky="w", pady=(10,2))
        #row += 1
        tab = tk.Frame(notebook)
        notebook.add(tab, text=group_name)
        tab_frames[group_name] = tab

        row = 0
        for key in keys:
            if key == "libName":
                continue  # We'll insert libName right after simulator
            #label = tk.Label(vars_frame, text=key)
            #label.grid(row=row, column=0, sticky="w")
            label = tk.Label(tab, text=key)
            label.grid(row=row, column=0, sticky="w")
            if key == "simulator":
                var = tk.StringVar(value=init_vars[key])
                entry = tk.OptionMenu(tab, var, *simulator_options)
                entry.config(width=27)
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = var
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                # Attach callback for showing/hiding libName
                var.trace_add("write", update_libname_visibility)
                row += 1
                # Insert libName field immediately after simulator
                libname_label = tk.Label(tab, text="libName")
                libname_entry = tk.Entry(tab, width=30)
                libname_entry.insert(0, str(init_vars["libName"]))
                entries["libName"] = libname_entry
                desc = var_descriptions.get("libName", "")
                if desc:
                    ToolTip(libname_label, desc)
                libname_label.grid(row=row, column=0, sticky="w")
                libname_entry.grid(row=row, column=1, padx=5, pady=2)
                row += 1
                continue  # Already incremented row for libName
            elif key == "simulatorPath":
                entry = tk.Entry(tab, width=30)
                entry.insert(0, str(init_vars[key]))
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                browse_btn = tk.Button(tab, text="Browse", command=lambda e=entry: browse_path(e))
                browse_btn.grid(row=row, column=3, padx=2)
                entries[key] = entry
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1    
            elif key == "pathName":
                entry = tk.Entry(tab, width=30)
                entry.insert(0, str(init_vars[key]))
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                browse_btn = tk.Button(tab, text="Browse", command=lambda e=entry: browse_path(e))
                browse_btn.grid(row=row, column=3, padx=2)
                entries[key] = entry
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1                        
            elif key == "corner":
                var = tk.StringVar(value=init_vars[key])
                entry = tk.OptionMenu(tab, var, *corner_options)
                entry.config(width=27)
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = var
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1
            elif key in ["sim", "view", "write"]:
                var = tk.StringVar(value=str(init_vars[key]))
                entry = tk.OptionMenu(tab, var, "True", "False")
                entry.config(width=27)
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = var
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1
            elif key in ["ports", "sides"]:
                var = tk.StringVar(value=str(init_vars[key]))
                entry = tk.OptionMenu(tab, var, *ports_sides_options)
                entry.config(width=27)
                entry.grid(row=row, column=1, padx=5, pady=2)
                # Custom entry for "Other"
                custom_entry = tk.Entry(tab, width=10)
                custom_entry.grid(row=row, column=3, padx=2)
                custom_entry.grid_remove()  # Hide initially
                # Set initial value if not in options
                if str(init_vars[key]) not in ports_sides_options:
                    var.set("Other")
                    custom_entry.insert(0, str(init_vars[key]))
                    custom_entry.grid()
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = (var, custom_entry)
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                # Show/hide custom entry on change
                var.trace_add("write", lambda *args, v=var, e=custom_entry, o=entry: show_custom_entry(v, e, o))
                row += 1
            elif key == "shape":
                shape_options = [str(i) for i in range(1, 27)]
                var = tk.StringVar(value=str(init_vars[key]))
                entry = tk.OptionMenu(tab, var, *shape_options)
                entry.config(width=27)
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = var
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1
            else:
                entry = tk.Entry(tab, width=30)
                entry.insert(0, str(init_vars[key]))
                entry.grid(row=row, column=1, padx=5, pady=2)
                unit = var_units.get(key, "")
                unit_label = tk.Label(tab, text=unit)
                unit_label.grid(row=row, column=2, sticky="w")
                entries[key] = entry
                desc = var_descriptions.get(key, "")
                if desc:
                    ToolTip(label, desc)
                row += 1

    # Set initial visibility
    update_libname_visibility()

    def on_ok():
        # Save all variables from GUI to init_vars
        for key, entry in entries.items():
            if key == "simulator":
                val = entry.get()
            elif key in ["ports", "sides"]:
                var, custom_entry = entry
                val = var.get()
                if val == "Other":
                    val = custom_entry.get()
                try:
                    init_vars[key] = int(val)
                except ValueError:
                    init_vars[key] = val
                continue
            elif key in ["sim", "view", "write"]:
                val = entry.get()
                if val.lower() == "true":
                    init_vars[key] = True
                elif val.lower() == "false":
                    init_vars[key] = False
                else:
                    init_vars[key] = val
                continue
            else:
                val = entry.get()
                if val.lower() == "true":
                    init_vars[key] = True
                elif val.lower() == "false":
                    init_vars[key] = False
                else:
                    try:
                        if '.' in val or 'e' in val.lower():
                            init_vars[key] = float(val)
                        else:
                            init_vars[key] = int(val)
                    except ValueError:
                        init_vars[key] = val
        # Ensure output folders exist
        ensure_output_folders(init_vars["pathName"])
        window.destroy()

    def on_cancel():
        # Set a flag or return None
        global init_vars
        init_vars = None
        window.destroy()

    button_frame = tk.Frame(frame)
    button_frame.pack(pady=5)
    tk.Button(button_frame, text="OK", command=on_ok).pack(side="left", padx=5)
    tk.Button(button_frame, text="Cancel", command=on_cancel).pack(side="left", padx=5)
    tk.Button(button_frame, text="Save to XML", command=save_vars).pack(side="left", padx=5)
    tk.Button(button_frame, text="Load from XML", command=load_from_xml).pack(side="left", padx=5)

    status_label = tk.Label(frame, text="")
    status_label.pack(pady=5)

    window.grab_set()
    window.wait_window()
    return init_vars

def get_init_vars(parent=None):
    if parent is None:
        root = tk.Tk()
        root.withdraw()
        result = edit_variables(root)
        root.destroy()
        return result
    else:
        return edit_variables(parent)

if __name__ == "__main__":
    edit_variables()  # or get_init_vars()
import time, sys
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp1d
import tkinter as tk
from tkinter import ttk, messagebox
from initializeVariables import get_init_vars
from vt_rrfc import *
from PIL import Image, ImageTk
import threading

class MainApp:
    def __init__(self, master):
        self.master = master
        self.list_counter = 0
        self.xProto = 0
        self.yProto = 0
        self.rows = 0
        self.cols = 0
        self.launch_l_pixels = 0
        self.reSizePixMap = False
        self.initPixMap = None
        self.flatInitPixMap = None
        self.base_name = ""
        self.pixelSize = 8  # Default pixel size
        master.title("Pixelated RF Filter Design")

        self.notebook = ttk.Notebook(master)
        self.notebook.pack(fill='both', expand=True)

        # Setup Tab
        self.setup_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.setup_frame, text='Setup')
        self.init_vars = None
        self.setup_button = ttk.Button(self.setup_frame, text="Edit Variables", command=self.edit_vars)
        self.setup_button.pack(pady=20)
        self.vars_label = tk.Label(self.setup_frame, text="No variables loaded.")
        self.vars_label.pack(pady=10)
        self.vars_text = tk.Text(self.setup_frame, height=18, width=60, state='disabled')
        self.vars_text.pack(pady=5)

        # Geometry Tab
        self.geometry_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.geometry_frame, text='Geometry')

        self.x_pixels_var = tk.IntVar(value=30)
        self.y_pixels_var = tk.IntVar(value=30)
        self.el_ports_var = tk.DoubleVar(value=30)

        ttk.Label(self.geometry_frame, text="Number of Pixels (X):").grid(row=0, column=0, sticky="e")
        self.x_pixels_entry = ttk.Entry(self.geometry_frame, textvariable=self.x_pixels_var, width=10)
        self.x_pixels_entry.grid(row=0, column=1)

        ttk.Label(self.geometry_frame, text="Number of Pixels (Y):").grid(row=1, column=0, sticky="e")
        self.y_pixels_entry = ttk.Entry(self.geometry_frame, textvariable=self.y_pixels_var, width=10)
        self.y_pixels_entry.grid(row=1, column=1)

        ttk.Label(self.geometry_frame, text="Electrical Length (deg):").grid(row=2, column=0, sticky="e")
        self.el_ports_entry = ttk.Entry(self.geometry_frame, textvariable=self.el_ports_var, width=10)
        self.el_ports_entry.grid(row=2, column=1)

        self.result_label = ttk.Label(self.geometry_frame, text="Port Width: --\nPort Length: --")
        self.result_label.grid(row=3, column=0, columnspan=2, pady=10)

        self.init_pixelmap_button = ttk.Button(
            self.geometry_frame, text="Initialize Pixel Map", command=self.initialize_pixelmap
        )
        self.init_pixelmap_button.grid(row=12, column=0, columnspan=3, pady=10)

        # DC Connection Section
        self.dc_frame = ttk.LabelFrame(self.geometry_frame, text="Check Boxes for Which Ports You Would Like a DC Connections Between")
        self.dc_frame.grid(row=4, column=0, columnspan=2, pady=10, sticky="ew")
        self.dc_vars = {}  # {(i, j): tk.BooleanVar}
        
        # Update geometry when any entry changes
        self.x_pixels_var.trace_add("write", self.update_geometry)
        self.y_pixels_var.trace_add("write", self.update_geometry)
        self.el_ports_var.trace_add("write", self.update_geometry)

        # Also update DC checkboxes when variables are loaded/changed
        self.num_ports = 2  # Default
        self.update_dc_checkboxes()
        self.dc_dontcare_var = tk.BooleanVar(value=False)
        self.dc_dontcare_cb = ttk.Checkbutton(
            self.geometry_frame, text="Don't Care (no DC connections enforced)", variable=self.dc_dontcare_var, command=self.update_dc_checkboxes
        )
        self.dc_dontcare_cb.grid(row=4, column=1, sticky="w", padx=5)

        # Symmetry dropdown
        self.symmetry_var = tk.StringVar(value="none")
        ttk.Label(self.geometry_frame, text="Symmetry:").grid(row=5, column=0, sticky="e", pady=(5,0))
        self.symmetry_menu = ttk.OptionMenu(
            self.geometry_frame, self.symmetry_var, "none", "xy-axis", "x-axis", "y-axis", "none"
        )
        self.symmetry_menu.grid(row=5, column=1, sticky="w", padx=0, pady=(5,0))

        # Custom name question
        self.custom_name_var = tk.StringVar(value="no")
        ttk.Label(self.geometry_frame, text="Do you want to give a custom name for your files?").grid(row=6, column=0, sticky="w", pady=(10,0))
        self.custom_name_yes = ttk.Radiobutton(self.geometry_frame, text="Yes", variable=self.custom_name_var, value="yes", command=self.update_custom_name_entry)
        self.custom_name_no = ttk.Radiobutton(self.geometry_frame, text="No", variable=self.custom_name_var, value="no", command=self.update_custom_name_entry)
        self.custom_name_yes.grid(row=7, column=0, sticky="w")
        self.custom_name_no.grid(row=7, column=1, sticky="w")

        # Entry for custom name
        self.custom_name_entry_var = tk.StringVar()
        self.custom_name_entry = ttk.Entry(self.geometry_frame, textvariable=self.custom_name_entry_var, width=40)
        self.custom_name_entry.grid(row=8, column=0, columnspan=2, sticky="w", padx=(0,5))
        self.custom_name_entry.grid_remove()

        # Pixel map existence question
        self.pixelmap_var = tk.StringVar(value="no")
        ttk.Label(self.geometry_frame, text="Does a pixel map already exist?").grid(row=9, column=0, sticky="w", pady=(10,0))
        self.pixelmap_yes = ttk.Radiobutton(self.geometry_frame, text="Yes", variable=self.pixelmap_var, value="yes", command=self.update_pixelmap_entry)
        self.pixelmap_no = ttk.Radiobutton(self.geometry_frame, text="No", variable=self.pixelmap_var, value="no", command=self.update_pixelmap_entry)
        self.pixelmap_yes.grid(row=10, column=0, sticky="w")
        self.pixelmap_no.grid(row=10, column=1, sticky="w")

        # Entry and browse for pixel map path
        self.pixelmap_path_var = tk.StringVar()
        self.pixelmap_entry = ttk.Entry(self.geometry_frame, textvariable=self.pixelmap_path_var, width=40)
        self.pixelmap_browse = ttk.Button(self.geometry_frame, text="Browse", command=self.browse_pixelmap)
        # Hide initially
        self.pixelmap_entry.grid(row=11, column=0, columnspan=2, sticky="w", padx=(0,5))
        self.pixelmap_browse.grid(row=11, column=2, sticky="w")
        self.pixelmap_entry.grid_remove()
        self.pixelmap_browse.grid_remove()

        self.pixelmap_img_label = tk.Label(self.geometry_frame)
        self.pixelmap_img_label.grid(row=0, column=2, rowspan=3, columnspan=3, pady=10)

        # Optimization Tab
        self.opt_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.opt_frame, text='Optimization')

        # Create a canvas and a vertical scrollbar for the Optimization tab
        self.opt_canvas = tk.Canvas(self.opt_frame, height=350)
        self.opt_scrollbar = ttk.Scrollbar(self.opt_frame, orient="vertical", command=self.opt_canvas.yview)
        self.opt_canvas.configure(yscrollcommand=self.opt_scrollbar.set)

        self.opt_canvas.grid(row=0, column=0, sticky="nsew")
        self.opt_scrollbar.grid(row=0, column=1, sticky="ns")

        # Make the canvas expandable
        self.opt_frame.grid_rowconfigure(0, weight=1)
        self.opt_frame.grid_columnconfigure(0, weight=1)

        # Create a frame inside the canvas to hold all widgets
        self.opt_inner_frame = ttk.Frame(self.opt_canvas)
        self.opt_inner_frame.bind("<Configure>", self._on_opt_frame_configure)
        self.opt_canvas.create_window((0, 0), window=self.opt_inner_frame, anchor="nw")

        # Number of Iterations
        ttk.Label(self.opt_inner_frame, text="Number of Iterations:").grid(row=0, column=0, sticky="e", pady=(10,0))
        self.num_iterations_var = tk.StringVar(value="3")
        self.num_iterations_entry = ttk.Entry(self.opt_inner_frame, textvariable=self.num_iterations_var, width=10)
        self.num_iterations_entry.grid(row=0, column=1, sticky="w", pady=(10,0))

        # Simultaneous Positions
        ttk.Label(self.opt_inner_frame, text="Simultaneous Positions:").grid(row=1, column=0, sticky="e", pady=(10,0))
        self.simul_positions_var = tk.StringVar(value="1")
        self.simul_positions_entry = ttk.Entry(self.opt_inner_frame, textvariable=self.simul_positions_var, width=10)
        self.simul_positions_entry.grid(row=1, column=1, sticky="w", pady=(10,0))

        # Dropdown for number of specs (now inside opt_inner_frame)
        specs_row = 2
        ttk.Label(self.opt_inner_frame, text="Number of Specifications:").grid(row=specs_row, column=0, sticky="e", pady=(10,0))
        self.num_specs_var = tk.StringVar(value="1")
        spec_options = [str(i) for i in range(1, 11)] + ["Other"]
        self.num_specs_menu = ttk.OptionMenu(
            self.opt_inner_frame, self.num_specs_var, "1", *spec_options, command=self.on_num_specs_change
        )
        self.num_specs_menu.grid(row=specs_row, column=1, sticky="w", pady=(10,0))

        self.other_specs_var = tk.StringVar()
        self.other_specs_entry = ttk.Entry(self.opt_inner_frame, textvariable=self.other_specs_var, width=5)
        self.other_specs_entry.grid(row=specs_row, column=2, sticky="w", padx=(5,0), pady=(10,0))
        self.other_specs_entry.grid_remove()
        self.other_specs_var.trace_add("write", lambda *args: self.update_spec_entries())

        # Frame to hold spec entries (inside opt_inner_frame)
        self.spec_entries_frame = ttk.Frame(self.opt_inner_frame)
        self.spec_entries_frame.grid(row=specs_row+1, column=0, columnspan=3, pady=(10,0), sticky="w")

        self.spec_entries = []
        self.update_spec_entries(1)

        # Run Tab
        self.run_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.run_frame, text='Run')
        self.run_button = ttk.Button(self.run_frame, text="Run Optimization", command=self.run_optimization)
        self.run_button.pack(pady=20)
        self.status_text = tk.Text(self.run_frame, height=15, width=80)
        self.status_text.pack(pady=10)

    # Bind the frame to update the scrollregion
    def _on_opt_frame_configure(self, event):
        self.opt_canvas.configure(scrollregion=self.opt_canvas.bbox("all"))

    def get_connect_map(self):
        """
        Returns a list of 6 entries corresponding to DC connections:
        [1-2, 1-3, 1-4, 2-3, 2-4, 3-4]
        If 'Don't Care' is checked, returns an empty list.
        """
        if hasattr(self, "dc_dontcare_var") and self.dc_dontcare_var.get():
            return []
        pairs = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
        connect_map = []
        for pair in pairs:
            var = self.dc_vars.get(pair)
            connect_map.append(1 if var and var.get() else 0)
        return connect_map

    def on_num_specs_change(self, *args):
        if self.num_specs_var.get() == "Other":
            self.other_specs_entry.grid()
        else:
            self.other_specs_entry.grid_remove()
        self.update_spec_entries()

    def update_spec_entries(self, *args):
        # Remove old widgets
        for widgets in self.spec_entries:
            for w in widgets:
                w.destroy()
        self.spec_entries = []

        # Determine number of specs
        try:
            if self.num_specs_var.get() == "Other":
                num_specs = int(self.other_specs_var.get())
            else:
                num_specs = int(self.num_specs_var.get())
        except Exception:
            num_specs = 0

        for i in range(num_specs):
            row = i
            label = ttk.Label(self.spec_entries_frame, text=f"Spec {i+1}:")
            start_entry = ttk.Entry(self.spec_entries_frame, width=12)
            stop_entry = ttk.Entry(self.spec_entries_frame, width=12)
            attn_entry = ttk.Entry(self.spec_entries_frame, width=12)
            maxmin_var = tk.StringVar(value="max")
            maxmin_menu = ttk.OptionMenu(self.spec_entries_frame, maxmin_var, "max", "max", "min")

            label.grid(row=row, column=0, sticky="w")
            ttk.Label(self.spec_entries_frame, text="Start Freq (Hz):").grid(row=row, column=1, sticky="e")
            start_entry.grid(row=row, column=2)
            ttk.Label(self.spec_entries_frame, text="Stop Freq (Hz):").grid(row=row, column=3, sticky="e")
            stop_entry.grid(row=row, column=4)
            ttk.Label(self.spec_entries_frame, text="Attenuation (dB):").grid(row=row, column=5, sticky="e")
            attn_entry.grid(row=row, column=6)
            ttk.Label(self.spec_entries_frame, text="Max/Min:").grid(row=row, column=7, sticky="e")
            maxmin_menu.grid(row=row, column=8)
            
            self.spec_entries.append((label, start_entry, stop_entry, attn_entry, maxmin_var, maxmin_menu))

    def update_custom_name_entry(self):
        if self.custom_name_var.get() == "yes":
            self.custom_name_entry.grid()
        else:
            self.custom_name_entry.grid_remove()
            self.custom_name_entry_var.set("")

    def update_pixelmap_entry(self):
        if self.pixelmap_var.get() == "yes":
            self.pixelmap_entry.grid()
            self.pixelmap_browse.grid()
        else:
            self.pixelmap_entry.grid_remove()
            self.pixelmap_browse.grid_remove()
            self.pixelmap_path_var.set("")

    def browse_pixelmap(self):
        from tkinter import filedialog
        path = filedialog.askopenfilename(title="Select Pixel Map CSV", filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if path:
            self.pixelmap_path_var.set(path)

    def get_csv_file(self):
        """Returns the value for csv_file based on user input in the geometry tab."""
        if self.pixelmap_var.get() == "yes" and self.pixelmap_path_var.get():
            return self.pixelmap_path_var.get()
        else:
            return 0

    def initialize_pixelmap(self):
        try:
            # Gather variables from GUI and init_vars
            if not self.init_vars:
                messagebox.showerror("Error", "Please set up variables first in the Setup tab.")
                return

            # --- Gather variables from Setup tab ---
            layoutUnit = self.init_vars.get("layoutUnit", 25.4e-6)
            ports = self.init_vars.get("ports", 2)
            sides = self.init_vars.get("sides", 2)
            shape = self.init_vars.get("shape", 1)
            corner = self.init_vars.get("corner", "normal")
            minPixel = self.init_vars.get("minPixel", 6)
            pixelSize = self.init_vars.get("pixelSize", 8)
            layoutRes = self.init_vars.get("layoutRes", 1)
            scale = self.init_vars.get("scale", 1)
            simulator = self.init_vars.get("simulator", "ADS")
            view = self.init_vars.get("view", False)
            write = self.init_vars.get("write", True)
            outFile = self.init_vars.get("pathName", "./") + "/pixelmap_init"
            connectMap = self.get_connect_map()
            sym = self.symmetry_var.get()
            portPosition = ''
            
            # Microstrip substrate and z0
            t = self.init_vars.get("t", 1.4)
            cond = self.init_vars.get("cond", 5.88e7)
            h = self.init_vars.get("h", 30)
            er = self.init_vars.get("er", 3.66)
            fc = self.init_vars.get("fc", 8.0e9)
            z0 = self.init_vars.get("z0", 50)
            sub1 = MicrostripSub(t, cond, h, er, fc)

            self.layoutUnit = layoutUnit
            self.ports = ports
            self.sides = sides
            self.corner = corner
            self.connectMap = connectMap
            self.minPixel = minPixel
            self.pixelSize = pixelSize
            self.layoutRes = layoutRes
            self.scale = scale
            self.simulator = simulator
            self.sym = sym
            self.z0 = z0
            self.sub1 = sub1
            self.pathName = self.init_vars.get("pathName", "./")
            
            simulator = self.init_vars.get("simulator", "ADS")
            simulatorPath = self.init_vars.get("simulatorPath", "")
            libName = self.init_vars.get("libName", "")

            self.simulator = simulator
            self.simulatorPath = simulatorPath
            self.libName = libName

            # --- Gather variables from Geometry tab ---
            x_pixels = self.x_pixels_var.get()
            y_pixels = self.y_pixels_var.get()
            el_ports = self.el_ports_var.get()
            w_l, l_l = self.w_l, self.l_l  # Port width and length from geometry calculation
            self.el_ports = el_ports

            # Pixel geometry
            self.pixelSize = pixelSize
            refPixelSize = pixelSize  # or set as needed
            launch_pixels = round(l_l/refPixelSize)        # or set as needed

            launch_l_pixels = int(launch_pixels*refPixelSize/pixelSize)
            xProto = x_pixels * int(refPixelSize / pixelSize) #- 2 * launch_l_pixels
            yProto = y_pixels * int(refPixelSize / pixelSize)
            self.xProto = xProto
            self.yProto = yProto
            
            rrfc1 = RandomComponent(
                        unit=layoutUnit,
                        ports=ports,
                        sides=sides,
                        corner=corner,
                        connect=connectMap,
                        minPix=minPixel,
                        pixelSize=pixelSize,
                        layoutRes=layoutRes,
                        scale=scale,
                        launchLen=el_ports,
                        seed=1,
                        sim=simulator,
                        view=view,
                        write=write,
                        outF=outFile,
                        sym=sym,
                        shape=shape,
                        portPosition=portPosition
                    )
            
            if self.custom_name_var.get() == "yes" and self.custom_name_entry_var.get():
                base_name = self.custom_name_entry_var.get()
            else:
                base_name = f"{ports}Port_{x_pixels}x{y_pixels}_pixelsize={pixelSize}_init"

            self.base_name = base_name
            outFile = self.pathName + "/" + base_name
            
            csv_file = self.get_csv_file()
            base = csv_file
            self.gds_file = self.pathName + 'data/' + base_name + ".gds"
            print(self.gds_file)
            if csv_file and os.path.isfile(csv_file):
                if base.lower().endswith('.csv'):
                    base = csv_file[:-4]  # Remove .csv extension for output file naming
                else:
                    base = csv_file
                refPixMap = np.loadtxt(csv_file, delimiter=',')
                if refPixelSize >= pixelSize:
                    reSizePixMap = np.repeat(refPixMap, int(refPixelSize/pixelSize), axis=0)
                    reSizePixMap = np.repeat(reSizePixMap, int(refPixelSize/pixelSize), axis=1)
                    xProto = np.size(reSizePixMap, 0)
                    yProto = np.size(reSizePixMap, 1)
                    launch_l_pixels = int(launch_pixels*refPixelSize/pixelSize)
                    xProto = xProto - 2 * launch_l_pixels
                    portPosition, _, _, _, _, _, _ = rrfc1.random_gds_dim(
                        sub1, xProto * pixelSize, yProto * pixelSize, z0)
                else:
                    a = np.zeros((int(np.size(refPixMap,0)/int(pixelSize/refPixelSize)),int(np.size(refPixMap,1)/int(pixelSize/refPixelSize)),int(pixelSize/refPixelSize)),dtype=int)
                    for p in range(1,int(pixelSize/refPixelSize)):
                        a[:,:,p] = refPixMap[p::int(pixelSize/refPixelSize),p::int(pixelSize/refPixelSize)]
                    reSizePixMap = np.mean(a, axis=2)
                    reSizePixMap = np.where(reSizePixMap >= 0.5, 1, 0)
                    yProto = np.size(reSizePixMap,0)
                    xProto = np.size(reSizePixMap,1)
                    launch_l_pixels = int(launch_pixels*refPixelSize/pixelSize)
                    xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
                    portPosition, _, _, _, _, _, _ = rrfc1.random_gds_dim(
                        sub1, xProto*pixelSize, yProto*pixelSize, z0)
                self.portPosition = portPosition
                # Optionally, display a message or update the GUI
                messagebox.showinfo("Pixel Map Initialized", f"Pixel map read!\nCSV file: {csv_file}")

                # --- Render PNG using print_image ---
                # Create a RandomUstrip object with the same parameters as rrfc1
                            
                rrfc1.print_image(sub1,z0)
                # If print_image saves a PNG, get the path
                png_path = outFile + ".png"
                # If print_image does not save, you can use matplotlib to save:
                # plt.savefig(png_path)

                # --- Display PNG in geometry tab ---
                img = Image.open(png_path)
                img = img.resize((256, 256))  # Resize as needed
                self.pixelmap_img = ImageTk.PhotoImage(img)
                self.pixelmap_img_label.config(image=self.pixelmap_img)
                self.pixelmap_img_label.image = self.pixelmap_img  # Prevent garbage collection
            else: #No existing pixel map, generate a new one
                # --- Generate pixel map ---
                reSizePixMap = False
                max_attempts = 1000
                attempt = 0
                csv_file = ""
                while attempt < max_attempts:
                    self.x = datetime.now()
                    # --- Create RandomComponent object ---
                    rrfc2 = RandomComponent(
                        unit=layoutUnit,
                        ports=ports,
                        sides=sides,
                        corner=corner,
                        connect=connectMap,
                        minPix=minPixel,
                        pixelSize=pixelSize,
                        layoutRes=layoutRes,
                        scale=scale,
                        launchLen=el_ports,
                        seed=self.x,
                        sim=simulator,
                        view=view,
                        write=write,
                        outF=outFile,
                        sym=sym,
                        shape=shape,
                        portPosition=portPosition
                    )
                    
                    portPosition, xBoard, yBoard, csv_file, gds_file, cell, _ = rrfc2.random_gds_dim(
                        sub1, xProto * pixelSize, yProto * pixelSize, z0
                    )
                    self.portPosition = portPosition
                    # Check if csv_file is a non-empty file
                    if csv_file and os.path.isfile(csv_file) and os.path.getsize(csv_file) > 0:
                        break
                    attempt += 1

                if not csv_file or not os.path.isfile(csv_file) or os.path.getsize(csv_file) == 0:
                    messagebox.showerror(
                        "Error",
                        "No design satisfied constraints, likely due to connectivity constraint"
                    )
                    return

                ## Define optimization parameters
                rows = int(yProto)
                cols = int(xProto)+2*launch_l_pixels
                self.rows = rows
                self.cols = cols
                self.launch_l_pixels = launch_l_pixels
                self.reSizePixMap = reSizePixMap
                pixels = rows * (cols-2*launch_l_pixels)
                list_counter = 0

                # Optionally, display a message or update the GUI
                messagebox.showinfo("Pixel Map Initialized", f"Pixel map created!\nCSV file: {csv_file}")

                print('Rows=',rows,' Cols=',cols)
                # Initialize pixel map and reshape to flat array
                if np.any(reSizePixMap) != 0: 
                    flatInitPixMap = reSizePixMap[:,launch_l_pixels:cols-launch_l_pixels]
                    flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))
                else:
                    print(csv_file)
                    initPixMap = np.loadtxt(csv_file, delimiter=',')
                    print('Pixmap Shape=',initPixMap.shape)
                    flatInitPixMap = initPixMap[:,launch_l_pixels:cols-launch_l_pixels]
                    print('Pixmap ReShape=',flatInitPixMap.shape)
                    flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))
                    self.initPixMap = initPixMap
                    self.flatInitPixMap = flatInitPixMap

                # --- Render PNG using print_image ---
                # Create a RandomUstrip object with the same parameters as rrfc1
                            
                rrfc2.print_image(sub1,z0)
                # If print_image saves a PNG, get the path
                png_path = outFile + ".png"
                # If print_image does not save, you can use matplotlib to save:
                # plt.savefig(png_path)

                # --- Display PNG in geometry tab ---
                img = Image.open(png_path)
                img = img.resize((256, 256))  # Resize as needed
                self.pixelmap_img = ImageTk.PhotoImage(img)
                self.pixelmap_img_label.config(image=self.pixelmap_img)
                self.pixelmap_img_label.image = self.pixelmap_img  # Prevent garbage collection

        except Exception as e:
            messagebox.showerror("Error", f"Failed to initialize pixel map:\n{e}")
                    
    def update_vars_display(self):
        self.vars_text.config(state='normal')
        self.vars_text.delete(1.0, tk.END)
        if self.init_vars:
            for k, v in self.init_vars.items():
                self.vars_text.insert(tk.END, f"{k}: {v}\n")
        else:
            self.vars_text.insert(tk.END, "No variables loaded.")
        self.vars_text.config(state='disabled')
    def run_optimization(self):
        if not self.init_vars:
            messagebox.showerror("Error", "Please set up variables first in the Setup tab.")
            return
        # Run optimization in a separate thread to keep GUI responsive
        threading.Thread(target=self.optimization_task, daemon=True).start()

    def edit_vars(self):
        self.init_vars = get_init_vars(self.master)
        self.vars_label.config(text="Variables loaded. Ready to run.")
        self.update_vars_display()
        self.update_geometry()  # Update geometry display after editing variables

    def update_geometry(self, *args):
        if not self.init_vars:
            self.result_label.config(text="Port Width: --\nPort Length: --")
            self.w_l = None
            self.l_l = None
            return
        try:
            t = self.init_vars.get("t", 1.4)
            cond = self.init_vars.get("cond", 5.88e7)
            h = self.init_vars.get("h", 30)
            er = self.init_vars.get("er", 3.66)
            fc = self.init_vars.get("fc", 8.0e9)
            z0 = self.init_vars.get("z0", 50)
            el = self.el_ports_var.get()
            sub1 = MicrostripSub(t, cond, h, er, fc)
            w, l = MicrostripCalc.synth_microstrip(sub1, z0, el)
            self.w_l = w
            self.l_l = l
            self.result_label.config(
                text=f"Port Width: {w:.2f} mil\nPort Length: {l:.2f} mil"
            )
        except Exception as e:
            self.result_label.config(text=f"Error: {e}")
            self.w_l = None
            self.l_l = None

        # Update DC checkboxes if number of ports changed
        ports = self.init_vars.get("ports", 2)
        if ports != self.num_ports:
            self.num_ports = ports
            self.update_dc_checkboxes()

    def update_dc_checkboxes(self):
        # Clear previous checkboxes
        for widget in self.dc_frame.winfo_children():
            widget.destroy()
        self.dc_vars = {}

        ports = self.init_vars.get("ports", 2) if self.init_vars else 2
        row = 0
        for i in range(1, ports):
            for j in range(i+1, ports+1):
                var = tk.BooleanVar()
                cb = ttk.Checkbutton(self.dc_frame, text=f"Ports {i} to {j}", variable=var)
                cb.grid(row=row, column=0, sticky="w")
                self.dc_vars[(i, j)] = var
                row += 1

        # Disable checkboxes if Don't Care is checked
        if hasattr(self, "dc_dontcare_var") and self.dc_dontcare_var.get():
            for cb in self.dc_frame.winfo_children():
                cb.configure(state="disabled")
        else:
            for cb in self.dc_frame.winfo_children():
                cb.configure(state="normal")        

    def run_optimization(self):
        if not self.init_vars:
            messagebox.showerror("Error", "Please set up variables first in the Setup tab.")
            return
        # Run optimization in a separate thread to keep GUI responsive
        threading.Thread(target=self.optimization_task, daemon=True).start()

    def optimization_task(self):
        if self.xProto == 0 or self.yProto == 0:
            messagebox.showerror("Error", "Pixel map not initialized. Please initialize pixel map first.")
            return
        pixels = self.xProto * self.yProto
        rows = self.rows
        cols = self.cols
        flatInitPixMap = self.flatInitPixMap if hasattr(self, 'flatInitPixMap') else None
        if flatInitPixMap is None:
            messagebox.showerror("Error", "Initial pixel map not set. Please initialize pixel map first.")
            return
        simulPositions = int(self.simul_positions_var.get())
        max_iterations = int(self.num_iterations_var.get())
        min_freq = float(self.init_vars.get("min_freq", 0))
        max_freq = float(self.init_vars.get("max_freq", 20e9))
        self.spIntFreq = np.arange(min_freq, max_freq + 2e6, 2e6)
        specs = []
        for widgets in self.spec_entries:
            _, start_entry, stop_entry, attn_entry, maxmin_var, _ = widgets
            try:
                start_freq = float(start_entry.get())
                stop_freq = float(stop_entry.get())
                attenuation = float(attn_entry.get())
                maxmin = maxmin_var.get()
                specs.append({
                    "start_freq": start_freq, 
                    "stop_freq": stop_freq, 
                    "attenuation": attenuation, 
                    "maxmin": maxmin
                })
            except ValueError:
                continue  # Skip invalid entries
        spec_indices = []    
        for spec in specs:
            start_idx = np.searchsorted(self.spIntFreq, spec["start_freq"])
            stop_idx = np.searchsorted(self.spIntFreq, spec["stop_freq"])
            spec_indices.append((start_idx, stop_idx, spec["attenuation"], spec["maxmin"]))
        print(spec)    
        # Import and run your optimization code here, using self.init_vars
        # For demonstration, just print to the status_text
        self.status_text.insert(tk.END, "Starting optimization...\n")
        self.status_text.see(tk.END)
        # ... Place your optimization code here ...
        # Example:
        # result = run_your_optimization(self.init_vars)
        # self.status_text.insert(tk.END, f"Result: {result}\n")

        self.DBS = dbsAlgo(pixels, self.cost_function, rows, 'none', 'max', simulPositions, max_iterations, self.call_back, initial_solution=flatInitPixMap)
        self.status_text.insert(tk.END, "Optimization complete!\n")
        self.status_text.see(tk.END)

    def cost_function(self, pixMap):
        reSizePixMap = self.reSizePixMap
        initPixMap = self.initPixMap
        pixelSize = self.pixelSize
        rows = self.rows
        cols = self.cols
        launch_l_pixels = self.launch_l_pixels
        pOps = pixOps(pixelSize = pixelSize)
        newPixMap = pixOps.updatePixels(pOps,pixMap)
        dbsPixMap = newPixMap.reshape(rows,cols-2*launch_l_pixels)
        if np.any(reSizePixMap) != 0:
            dbsPixMap = np.concatenate((reSizePixMap[:,0:launch_l_pixels],dbsPixMap,reSizePixMap[:,cols-launch_l_pixels:cols]),1)
        else:
            dbsPixMap = np.concatenate((initPixMap[:,0:launch_l_pixels],dbsPixMap,initPixMap[:,cols-launch_l_pixels:cols]),1)

        outFile = self.pathName + 'data/' + self.base_name[:-5] + '_' + str(self.list_counter)
        csv_file2 = outFile + ".csv"
        # Export Pixel Map file
        np.savetxt(csv_file2, dbsPixMap, fmt = '%d', delimiter = ",")  

        #rfc2 = rfc(unit=layoutUnit,pixelSize=pixelSize,sim=simulator,\
        #      view=False,write=True,outF=csv_file2)
        rfc2 = RandomUstrip(
            unit=self.layoutUnit,
            ports=self.ports,
            sides=self.sides,
            corner=self.corner,
            connect=self.connectMap,
            minPix=self.minPixel,
            pixelSize=self.pixelSize,
            layoutRes=self.layoutRes,
            scale=self.scale,
            launchLen=self.el_ports,
            seed=self.x,
            sim=self.simulator,
            view=False,
            write=True,
            outF=csv_file2,
            sym=self.sym,
            portPosition='')
        cell = rfc2.create_gds_from_pixmap_file(self.sub1,self.z0,None)
        
        gds_file2 = self.gds_file.replace('_init.gds', '_' + str(self.list_counter) + '.gds')

        print(gds_file2)
        em1 = EmSim(
            workingPath = self.pathName, 
            adsLibName = self.libName, 
            gdsFile = gds_file2,
            csvFile = csv_file2,
            numPorts = self.ports,
            portPositions = self.portPosition,
            gdsCellName = cell,
            dataFile = outFile,
            simPath = self.simulatorPath)
        
        em1.mom_run()

        # Read the init s-parameter data from the file and pixelMap
        dataPath = self.pathName + 'data/spfiles/afs/' + self.base_name[:-5] + '_' + str(self.list_counter) + '.afs'
        initFreq, initS = readCiti(dataPath)
        reFreq = np.arange(0,20.002e9,2e6)
        s11interp = interp1d(initFreq,initS[:,0])
        s12interp = interp1d(initFreq,initS[:,1])
        s21interp = interp1d(initFreq,initS[:,2])
        s22interp = interp1d(initFreq,initS[:,3])
        reS11 = s11interp(spIntFreq)
        reS12 = s12interp(spIntFreq)
        reS21 = s21interp(spIntFreq)
        reS22 = s22interp(spIntFreq)

        aS11 = 20*np.log10(abs(reS11)) 
        aS21 = 20*np.log10(abs(reS21))
        aS22 = 20*np.log10(abs(reS22))
        
        #AS11_pass = np.log10(EmSim.sigmoid(abs(aS11) - 10))
        #AS22_pass = np.log10(EmSim.sigmoid(abs(aS22) - 10))
        #AS21_pass = np.log10(EmSim.sigmoid(abs(aS21) - 0))
        #AS11_stop = np.log10(EmSim.sigmoid(abs(aS11) - 0))
        #AS22_stop = np.log10(EmSim.sigmoid(abs(aS22) - 0))
        #AS21_stop = np.log10(EmSim.sigmoid(abs(aS21) - 20))

        fs1_s11 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS11[fs1s_index:fs1e_index])/20))**2)))
        fs1_s21 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS21[fs1s_index:fs1e_index])/20))**2)))
        fs1_s22 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS22[fs1s_index:fs1e_index])/20))**2)))

        fp1_s11 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS11[fp1s_index:fp1e_index])/20))**2)))
        fp1_s21 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS21[fp1s_index:fp1e_index])/20))**2)))
        fp1_s22 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS22[fp1s_index:fp1e_index])/20))**2)))

        fs2_s11 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS11[fs2s_index:fs2e_index])/20))**2)))
        fs2_s21 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS21[fs2s_index:fs2e_index])/20))**2)))
        fs2_s22 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS22[fs2s_index:fs2e_index])/20))**2)))

        print('RMS S-params:\nS11(<7.124G)=' + str(fs1_s11) + \
                        '\nS21(<7.124G)=' + str(fs1_s21) + \
                        '\nS22(<7.124G)=' + str(fs1_s22) + \
                        '\nS11(7.124-9.124G)=' + str(fp1_s11) + \
                        '\nS21(7.124-9.124G)=' + str(fp1_s21) + \
                        '\nS22(7.124-9.124G)=' + str(fp1_s22) + \
                        '\nS11(>9.124G)=' + str(fs2_s11) + \
                        '\nS21(>9.124G)=' + str(fs2_s21) + \
                        '\nS22(>9.124G)=' + str(fs2_s22))

        # Passband rewards. Goal is to maximize S21 in the passband while minimizing S11 and S22. Added non-linear functions to 
        # reduce overfitting to a single parameter. The gain (S21) is weighted more heavily than the return loss (factor of 5 here)
        rw1 = abs(fp1_s11) / (1 + np.exp(-2 + abs(fp1_s21))) + abs(fp1_s22) / (1 + np.exp(-2 + abs(fp1_s21))) + 5*(10**((fp1_s21 + 2)/20) / (1 + np.exp(15-abs(fp1_s11)) + np.exp(15-abs(fp1_s22))))

        # Stopband rewards. Goal is to minimize S21 in the stopband while maximizing S11 and S22. Added non-linear functions to 
        # reduce overfitting to a single parameter. The gain (S21) is weighted more heavily than the return loss (factor of 5 here)
        rw2 = 10**((-35-(fs1_s21))/20) #/ (1 + np.exp(3-abs(aS11[fs1_index])) + np.exp(3-abs(aS22[fs1_index])))
        rw3 = 10**((-35-(fs2_s21))/20) #/ (1 + np.exp(3-abs(aS11[fs1_index])) + np.exp(3-abs(aS22[fs1_index])))


        rw1_list[self.list_counter] = rw1
        rw2_list[self.list_counter] = rw2
        rw3_list[self.list_counter] = rw3

        self.list_counter += 1
        #fomUpdate = "Current RMSE Errors:" + "\nIL (ISM):" + str(rw1) + "\nIL (UNII):" + str(rw3) + \
        #                                     "\nATT (LOW):" + str(rw9) + "\nATT (MID):" + str(rw2) + \
        #                                     "\nATT (HIGH):" + str(rw8) + "\nRL (ISM):" + str(rw4+rw6) + \
        #                                     "\nRL (UNII):" + str(rw5+rw7) + "\nTotal Error:" + \
        #                                     str(rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8 + rw9) + "\n"
        fomUpdate = "Current RMSE Errors:" + "\nfp1:" + str(rw1) + "\nfs1:" + str(rw2) + "\nfs2:" + str(rw3) + "\nTotal Error:" + str(1 / (1/rw1 + 1/rw2 + 1/rw3)) + "\n"
        print(fomUpdate)
        return 1 / (1/rw1 + 1/rw2 + 1/rw3)
        #return rw1 + rw2 + rw3 + rw8 + rw9
        #return max(rw1, rw2, rw3, rw4, rw5, rw6, rw7, rw8, rw9) #minimax


    def call_back():
        print("Size of remained:", DBS.get_remained_size())
        print("Number of iteration:", DBS.get_iteration_number())
        print("The minimum cost:", DBS.get_cost())
        print("Best Solution:\n", DBS.get_best_solution())
        ## save results into files
        np.savetxt(self.pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW1", np.array(rw1_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW2", np.array(rw2_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW3", np.array(rw3_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        csvFinal = self.pathName + 'data/' + self.base_name[:-5] + '_finalDesign.csv'
        if np.any(reSizePixMap) != 0:
            np.savetxt(csvFinal, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
            reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
        else:
            np.savetxt(csvFinal, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
            initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

        DBS = dbsAlgo(pixels, cost_function, rows, 'none', 'max', simulPositions, max_iterations, call_back, initial_solution=flatInitPixMap)

        ## save initial solution into the file
        np.save("init_solution", DBS.best_solution)

        ## save initial solution into the file
        csvInit = self.pathName + 'data/' + self.base_name[:-5] + '_initialDesign.csv'
        if np.any(reSizePixMap) != 0:
            np.savetxt(csvInit, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
            reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
        else:
            np.savetxt(csvInit, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
            initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

        ## start DBS optimization
        DBS.run()
        
        ## save results into files
        np.save("solution", DBS.best_solution.reshape(rows,cols-2*launch_l_pixels))
        np.savetxt(self.pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW1", np.array(rw1_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW2", np.array(rw2_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        np.savetxt(self.pathName + "data/RW3", np.array(rw3_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

        ## plot for the cost
        x = range(0, DBS.cg_curve.shape[0])
        y = DBS.cg_curve
        xlabel = "Operation Times"
        ylabel = "Cost (Lower Better)"
        plt.figure()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x, y)
        plt.savefig("PBS", dpi=600)
        plt.show()
        plt.close()
        plt.clf()

if __name__ == "__main__":
    root = tk.Tk()
    app = MainApp(root)
    root.mainloop()
import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfilename
from typing import Any

CHECKBOX_PROPERTIES = ('LogP', 'LogD', 'LogS', 'tPSA', 'Toxicity', 'EGG', 'Calculate docking score')
DROPDOWN_PROPERTIES = ('All', 'LogP', 'LogD', 'LogS', 'tPSA', 'Toxicity', 'EGG')
FRAME_X_PADDING = 200
FRAME_Y_PADDING = 50

def get_filename() -> None:
    filename = askopenfilename()

    csv_filename.set(filename)

def make_resizable(tk_object: Any, nrows: int, ncols: int) -> None:
    for row in range(nrows):
        tk_object.rowconfigure(row, weight=1)
        for col in range(ncols):
            tk_object.columnconfigure(col, weight=1)

if __name__ == '__main__':
    window = tk.Tk()
    window.title('Flower Pot')
    make_resizable(window, 8, 4)

    logo_frame = tk.Frame()
    logo_frame.grid(row=0, column=1, rowspan=2, columnspan=2)
    logo = ttk.Label(logo_frame, text='Flower Pot')
    logo.pack()

    smiles_frame = tk.Frame()
    smiles_frame.grid(row=2, column=0, rowspan=6, columnspan=2, padx=FRAME_X_PADDING, pady=FRAME_Y_PADDING)
    smiles_label = ttk.Label(smiles_frame, text='SMILES Input:')
    smiles_label.grid(row=0, column=0, sticky='w')
    smiles_text = tk.Text(smiles_frame)
    smiles_text.grid(row=1, column=0, rowspan=10, sticky='nsew')
    go_button = ttk.Button(smiles_frame, text='Go', width=10)
    go_button.grid(row=1, column=1, sticky='nw')
    dropdown_value = tk.StringVar()
    smiles_dropdown = ttk.OptionMenu(smiles_frame, dropdown_value, DROPDOWN_PROPERTIES[0], *DROPDOWN_PROPERTIES)
    smiles_dropdown.config(width=10)
    smiles_dropdown.grid(row=2, column=1, sticky='nw')
    docking_box_value = tk.BooleanVar()
    smiles_docking_box = ttk.Checkbutton(smiles_frame, variable=docking_box_value, text='Calculate docking score')
    smiles_docking_box.grid(row=3, column=1, sticky='nw')
    make_resizable(smiles_frame, 11, 2)

    batch_frame = tk.Frame()
    batch_frame.grid(row=2, column=2, rowspan=6, columnspan=2, padx=FRAME_X_PADDING, pady=FRAME_Y_PADDING)
    csv_filename = tk.StringVar()
    csv_entry = ttk.Entry(batch_frame, textvariable=csv_filename)
    csv_entry.grid(row=0, column=0, sticky='nw')
    browse_button = ttk.Button(batch_frame, text='Browse', command=get_filename)
    browse_button.grid(row=0, column=1, sticky='nw')
    checkbox_indices = [(i, j) for i in range(3, 7) for j in range(2)]
    del checkbox_indices[-1]
    checkbox_values = {}
    for idx, prop in zip(checkbox_indices, CHECKBOX_PROPERTIES):
        i, j = idx
        checkbox_values[prop] = tk.BooleanVar()
        checkbox = ttk.Checkbutton(batch_frame, variable=checkbox_values[prop], text=prop)
        checkbox.grid(row=i, column=j, sticky='w')
    go_button = ttk.Button(batch_frame, text='Go')
    go_button.grid(row=6, column=1, sticky='nw')
    make_resizable(batch_frame, 5, 2)

    window.mainloop()

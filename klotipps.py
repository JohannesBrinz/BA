#simple GUI app (https://www.youtube.com/watch?v=jE-SpRI3K5g&list=RDCMUClb90NQQcskPUGDIXsQEz5Q&start_radio=1&rv=jE-SpRI3K5g&t=60)


import tkinter as tk
from tkinter import filedialog, Text
import os, sys, subprocess

root = tk.Tk()
apps = []

#Anlegen der Klo_Dokumente
klos = []


#Definition Variablen

neues_klo = "Musterklo"

#Funktionen

def anlegen():
    neues_klo = input("Gib den Namen des Klos ein, das du anlegen möchtest:    ")
    klos.append(neues_klo)
    print("Liste der Klos:  ", klos)
    for widget in frame.winfo_children():
        widget.destroy()
    for klo in klos:
        label = tk.Button(frame, text=klo, bg = "gray")
        label.pack()

#Rest
canvas = tk.Canvas(root, height=700, width=700, bg="brown")
canvas.pack()

frame = tk.Frame(root, bg = "white")
frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)


NeuesKlo = tk.Button(root, text="Neues Klo anlegen", padx=15, pady=5, fg="black", command = anlegen)
NeuesKlo.pack()




'''
def open_file(filename):
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        subprocess.call([opener, filename])

def addApp():
    for widget in frame.winfo_children():
        widget.destroy()

    filename = filedialog.askopenfilename(initialdir="/", title = "Select File",
                                        filetypes=(("executables", "*.exe"), ("all files", "*.*")))
    apps.append(filename)
    print(filename)
    for app in apps:
        label = tk.Label(frame, text=app, bg = "gray")
        label.pack()

def runApps():
    for app in apps:
        open_file(app)


canvas = tk.Canvas(root, height=700, width=700, bg="brown")
canvas.pack()

frame = tk.Frame(root, bg = "white")
frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

openFile = tk.Button(root, text="Open File", padx=10, pady=5, fg="white", command=addApp)
openFile.pack()

runApps = tk.Button(root, text="Run Apps", padx=15, pady=5, fg="white", command = runApps)
runApps.pack()'''

root.mainloop()

"""Small window shown while the program starts: the heavy imports (pandas,
matplotlib, scipy, the backend) run in a thread behind a spinner."""

import threading
import tkinter as tk

from .gui_texts import TEXTS
from .ui_helpers import FONT_SIZE, FONT_UI, PALETTE, Spinner


def load_with_splash(load):
    """Run load() in a thread while a splash window turns; return its result or
    raise its exception."""
    root = tk.Tk()
    root.overrideredirect(True)
    w, h = 340, 170
    root.geometry(f"{w}x{h}+{(root.winfo_screenwidth() - w) // 2}+{(root.winfo_screenheight() - h) // 2}")
    root.configure(bg=PALETTE['bg_panel'])
    tk.Label(root, text="EasyPAML", bg=PALETTE['bg_panel'], fg=PALETTE['text_primary'],
             font=(FONT_UI, -FONT_SIZE['xl'], 'bold')).pack(pady=(26, 10))
    spinner = Spinner(root, size=36, bg=PALETTE['bg_panel'])
    spinner.canvas.pack()
    tk.Label(root, text=TEXTS["loading_app"], bg=PALETTE['bg_panel'], fg=PALETTE['text_secondary'],
             font=(FONT_UI, -FONT_SIZE['sm'])).pack(pady=(10, 0))

    out = {}

    def work():
        try:
            out['value'] = load()
        except BaseException as exc:      # re-raised in the main thread
            out['error'] = exc
    thread = threading.Thread(target=work, daemon=True)
    thread.start()

    def poll():
        if thread.is_alive():
            root.after(40, poll)
        else:
            spinner.stop()
            root.destroy()
    root.after(40, poll)
    root.mainloop()
    if 'error' in out:
        raise out['error']
    return out.get('value')

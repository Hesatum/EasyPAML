"""Interface pieces shared by the main window and the results panel: palettes and
theme, dialogs, file and folder picker, window sizing."""

import fnmatch
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import customtkinter as ctk

from .gui_texts import TEXTS, get_language

_ON_LINUX = platform.system() == "Linux"
_ON_WIN = platform.system() == "Windows"
FONT_UI = "DejaVu Sans" if _ON_LINUX else "Roboto"
FONT_MONO = "Cascadia Code" if _ON_WIN else "DejaVu Sans Mono"

# Text colours have a contrast of at least 4.5:1 (WCAG AA) on every background
# of their palette (tests/test_gui_texts.py).
DARK_PALETTE = {
    'bg_dark': '#0d0d11',
    'bg_card': '#17171f',
    'bg_card_hover': '#20202c',
    'text_primary': '#eeeef2',
    'text_secondary': '#b8b8ca',
    'text_tertiary': '#a2a2b6',
    'text_muted': '#9c9cb2',
    'accent_text': '#818cf8',
    'accent_fill': '#4f46e5',
    'success_text': '#22c55e',
    'success_fill': '#15803d',
    'warning_text': '#f59e0b',
    'warning_fill': '#b45309',
    'danger_text': '#f87171',
    'danger_fill': '#b91c1c',
    'info_text': '#22d3ee',
    'info_fill': '#0e7490',
    'border': '#262632',
    'border_hover': '#3a3a4e',
    'bg_window': '#0a0b0e',          # fundo da janela
    'bg_panel': '#111318',           # painel lateral
    'bg_surface': '#181a21',
    'bg_elevated': '#20232c',
    'bg_elevated_hover': '#292c37',
    'bg_inset': '#0e1014',           # campos de entrada
    'divider': '#262a34',            # linha de 1 px entre grupos
    'control_border': '#2c303c',
    'control_border_hover': '#3b4050',
    'neutral_fill': '#343846',
    'switch_knob': '#e6e7ee',
    'switch_track': '#3b4050',
    'success_subtle': '#0f2a1d',
    'success_fg': '#86efac',         # texto sobre success_subtle (10,9:1)
    'warning_subtle': '#2e2108',
    'warning_fg': '#fcd34d',         # texto sobre warning_subtle (10,9:1)
    'danger_subtle': '#341719',
    'danger_fg': '#fca5a5',          # texto sobre danger_subtle (8,6:1)
    'row_alt': '#15171d',            # zebra suave das tabelas (entre bg_panel e bg_surface)
    'accent_blue': '#6366f1',
    'accent_cyan': '#22d3ee',
    'accent_purple': '#a78bfa',
    'success_light': '#86efac',
    'plot_bg': '#111115',
    'plot_fg': '#c8c8d4',
    'plot_muted': '#9898a6',
    'plot_line': '#5a5a6a',
    'plot_node': '#2a2a36',
    'canvas_bg': '#14141a',
    'canvas_line': '#666677',
    'canvas_text': '#e5e5e5',
}

LIGHT_PALETTE = {
    'bg_dark': '#eef0f3',
    'bg_card': '#ffffff',
    'bg_card_hover': '#f3f4f7',
    'text_primary': '#111827',
    'text_secondary': '#363c47',
    'text_tertiary': '#474e5b',
    'text_muted': '#474e5b',
    'accent_text': '#4338ca',
    'accent_fill': '#4f46e5',
    'success_text': '#167038',
    'success_fill': '#15803d',
    'warning_text': '#9a4d07',
    'warning_fill': '#b45309',
    'danger_text': '#b42318',
    'danger_fill': '#b91c1c',
    'info_text': '#0b6378',
    'info_fill': '#0e7490',
    'border': '#d5d9e0',
    'border_hover': '#b9bfca',
    'bg_window': '#eef0f3',
    'bg_panel': '#f6f7f9',
    'bg_surface': '#ffffff',
    'bg_elevated': '#f3f4f7',
    'bg_elevated_hover': '#e7e9ee',
    'bg_inset': '#ffffff',
    'divider': '#e1e4e9',
    'control_border': '#b4bbc6',
    'control_border_hover': '#aab1bd',
    'neutral_fill': '#4b5563',
    'switch_knob': '#4b5563',
    'switch_track': '#d5d9e0',
    'success_subtle': '#dcf5e4',
    'success_fg': '#14532d',
    'warning_subtle': '#fdf0c8',
    'warning_fg': '#7a3e06',
    'danger_subtle': '#fde4e2',
    'danger_fg': '#8f1d16',
    'row_alt': '#f7f8fa',
    'accent_blue': '#4f46e5',
    'accent_cyan': '#0e7490',
    'accent_purple': '#6d28d9',
    'success_light': '#167038',
    'plot_bg': '#ffffff',
    'plot_fg': '#1f2937',
    'plot_muted': '#4b5563',
    'plot_line': '#9ca3af',
    'plot_node': '#e5e7eb',
    'canvas_bg': '#ffffff',
    'canvas_line': '#6b7280',
    'canvas_text': '#111827',
}

# Active palette, filled by apply_theme() before any window is built.
PALETTE: dict = {}
THEME_CHOICES = ('system', 'light', 'dark')
_THEME_PREF = Path.home() / '.easypaml_theme'
CURRENT_THEME = {'choice': 'system', 'mode': 'dark'}


def load_theme_pref() -> str:
    try:
        choice = _THEME_PREF.read_text(encoding='utf-8').strip()
    except OSError:
        return 'system'
    return choice if choice in THEME_CHOICES else 'system'


def save_theme_pref(choice: str) -> None:
    try:
        _THEME_PREF.write_text(choice, encoding='utf-8')
    except OSError:
        pass


def system_theme() -> str:
    """'light' or 'dark' from the system (darkdetect); 'light' if unknown."""
    try:
        import darkdetect
        return 'dark' if (darkdetect.theme() or '').lower() == 'dark' else 'light'
    except Exception:
        return 'light'


def apply_theme(choice: str = None) -> str:
    """Define PALETTE e o modo do customtkinter. Retorna 'light' ou 'dark'."""
    choice = choice if choice in THEME_CHOICES else load_theme_pref()
    mode = system_theme() if choice == 'system' else choice
    PALETTE.clear()
    PALETTE.update(LIGHT_PALETTE if mode == 'light' else DARK_PALETTE)
    CURRENT_THEME.update(choice=choice, mode=mode)
    ctk.set_appearance_mode(mode)
    return mode


apply_theme()

# spacing (4 px grid) and corner radii
SPACE = {'xs': 4, 'sm': 8, 'md': 12, 'lg': 16, 'xl': 24, 'xxl': 32}
RADIUS = {'field': 6, 'card': 8, 'panel': 12}
# font sizes (pt), none below 11
FONT_SIZE = {'xs': 13, 'sm': 13, 'md': 14, 'lg': 16, 'xl': 18, 'xxl': 21}


def mix(color_a: str, color_b: str, t: float) -> str:
    """Mistura duas cores hex: t = 0 -> a, t = 1 -> b."""
    a = [int(color_a.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    b = [int(color_b.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    return '#' + ''.join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))


def hover_tint(accent: str, background: str = '#17171f') -> str:
    """Hover background of an outlined button: a dark shade of the accent, so
    accent-coloured text stays legible."""
    return mix(background, accent, 0.22)


def fit_to_screen(win, width: int, height: int, min_w: int = 1024, min_h: int = 640) -> None:
    """Open the window at the requested size, never larger than the screen."""
    win.update_idletasks()
    sw, sh = win.winfo_screenwidth(), win.winfo_screenheight()
    w = min(width, sw - 40)
    h = min(height, sh - 80)
    x = max(0, (sw - w) // 2)
    y = max(0, (sh - h) // 3)
    win.geometry(f"{w}x{h}+{x}+{y}")
    win.minsize(min(min_w, w), min(min_h, h))


def disable_mouse_wheel(widget) -> None:
    """Keep the mouse wheel from changing a CTkSlider while the panel scrolls."""
    canvas = getattr(widget, '_canvas', None)
    if canvas is None:
        return
    for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        try:
            canvas.unbind(seq)
        except Exception:
            pass


def open_folder(path, parent=None) -> bool:
    """Open a folder in the file manager. With parent, a failure (also one that
    xdg-open or open reports a moment later) shows the path, copied to the clipboard."""
    path = str(Path(path))

    def _failed():
        if parent is None:
            return
        try:
            parent.clipboard_clear()
            parent.clipboard_append(path)
        except Exception:
            pass
        show_message(parent, "EasyPAML", TEXTS["msg_open_folder_failed"].format(path=path), 'warning')

    try:
        if _ON_WIN:
            os.startfile(path)  # type: ignore[attr-defined]
            return True
        cmd = ['open', path] if sys.platform == 'darwin' else ['xdg-open', path]
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        _failed()
        return False
    if parent is not None:
        def _check(tries=20):
            code = proc.poll()
            if code is None and tries:
                parent.after(150, lambda: _check(tries - 1))
            elif code:
                _failed()
        parent.after(150, _check)
    return True


class _Modal(ctk.CTkToplevel):
    def __init__(self, parent, title: str, width: int, height: int):
        super().__init__(parent)
        self.title(title)
        self.configure(fg_color=PALETTE['bg_window'])
        fit_to_screen(self, width, height, min_w=420, min_h=240)
        self.transient(parent)
        self.bind("<Escape>", lambda e: self._close(None))
        self.protocol("WM_DELETE_WINDOW", lambda: self._close(None))
        self.result = None
        self.after(50, self._grab)

    def _grab(self):
        try:
            self.grab_set()
            self.focus_force()
        except Exception:
            pass

    def _close(self, value):
        self.result = value
        try:
            self.grab_release()
        except Exception:
            pass
        self.destroy()

    def show(self):
        self.wait_window()
        return self.result


def _button(parent, text, command, fill, **kw):
    return ctk.CTkButton(parent, text=text, command=command, fg_color=fill,
                         hover_color=mix(fill, '#000000', 0.2), text_color='#ffffff',
                         text_color_disabled=mix('#ffffff', fill, 0.45),
                         font=(FONT_UI, FONT_SIZE['md'], 'bold'), height=36,
                         corner_radius=RADIUS['card'], **kw)


def ask_yes_no(parent, title: str, message: str, yes: str = None, no: str = None) -> bool:
    """Yes/no dialog with translated buttons (Tk's messagebox uses the system's)."""
    dlg = _Modal(parent, title, 460, 220)
    ctk.CTkLabel(dlg, text=message, font=(FONT_UI, FONT_SIZE['md']), wraplength=400, justify='left',
                 text_color=PALETTE['text_primary']).pack(padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['lg']),
                                                          anchor='w')
    row = ctk.CTkFrame(dlg, fg_color='transparent')
    row.pack(fill='x', padx=SPACE['xl'], pady=(0, SPACE['xl']), side='bottom')
    _button(row, no or TEXTS['btn_no'], lambda: dlg._close(False), PALETTE['neutral_fill']).pack(side='right')
    _button(row, yes or TEXTS['btn_yes'], lambda: dlg._close(True), PALETTE['accent_fill']).pack(
        side='right', padx=(0, SPACE['sm']))
    return bool(dlg.show())


def ask_choice(parent, title: str, message: str, options) -> Optional[str]:
    """Dialog with one button per (label, value); returns the value, or None when closed."""
    dlg = _Modal(parent, title, 560, 260)
    ctk.CTkLabel(dlg, text=message, font=(FONT_UI, FONT_SIZE['md']), wraplength=500, justify='left',
                 text_color=PALETTE['text_primary']).pack(padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['lg']),
                                                          anchor='w')
    row = ctk.CTkFrame(dlg, fg_color='transparent')
    row.pack(fill='x', padx=SPACE['xl'], pady=(0, SPACE['xl']), side='bottom')
    for i, (label, value) in enumerate(reversed(list(options))):
        _button(row, label, lambda v=value: dlg._close(v),
                PALETTE['accent_fill'] if i == len(options) - 1 else PALETTE['neutral_fill']).pack(
            side='right', padx=(0, SPACE['sm'] if i else 0))
    return dlg.show()


def show_message(parent, title: str, message: str, kind: str = 'info') -> None:
    dlg = _Modal(parent, title, 560, 300)
    color = {'error': PALETTE['danger_text'], 'warning': PALETTE['warning_text']}.get(
        kind, PALETTE['text_primary'])
    box = ctk.CTkTextbox(dlg, font=(FONT_UI, FONT_SIZE['md']), fg_color=PALETTE['bg_surface'],
                         text_color=color, wrap='word', corner_radius=RADIUS['card'],
                         border_width=1 if kind in ('error', 'warning') else 0,
                         border_color=mix(PALETTE['bg_surface'], color, 0.45))
    box.pack(fill='both', expand=True, padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['md']))
    box.insert('end', message)
    box.configure(state='disabled')
    _button(dlg, "OK", lambda: dlg._close(True), PALETTE['accent_fill'], width=100).pack(
        pady=(0, SPACE['lg']))
    dlg.show()


class PreflightDialog(_Modal):
    """Problems found before running. Returns 'continue' or 'fix'."""

    def __init__(self, parent, report):
        super().__init__(parent, TEXTS['preflight_title'], 900, 620)
        lang = get_language()
        n_err = sum(1 for i in report.issues if i.severity == 'error')
        n_warn = sum(1 for i in report.issues if i.severity == 'warning')
        head = ctk.CTkFrame(self, fg_color='transparent')
        head.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['lg'], SPACE['sm']))
        ctk.CTkLabel(head, text=TEXTS['preflight_heading'].format(
                         genes=len(report.genes), errors=n_err, warnings=n_warn),
                     font=(FONT_UI, FONT_SIZE['lg'], 'bold'), text_color=PALETTE['text_primary'],
                     justify='left', wraplength=840).pack(anchor='w')
        ctk.CTkLabel(head, text=TEXTS['preflight_explain'], font=(FONT_UI, FONT_SIZE['sm']),
                     text_color=PALETTE['text_secondary'], justify='left',
                     wraplength=840).pack(anchor='w', pady=(SPACE['xs'], 0))

        box = ctk.CTkTextbox(self, font=(FONT_UI, FONT_SIZE['sm']), fg_color=PALETTE['bg_surface'],
                             text_color=PALETTE['text_primary'], wrap='word',
                             corner_radius=RADIUS['card'], border_width=1,
                             border_color=mix(PALETTE['bg_surface'], PALETTE['warning_text'], 0.35))
        box.pack(fill='both', expand=True, padx=SPACE['xl'], pady=SPACE['sm'])
        box.tag_config('error', foreground=PALETTE['danger_text'])
        box.tag_config('warning', foreground=PALETTE['warning_text'])
        box.tag_config('info', foreground=PALETTE['text_tertiary'])
        box.tag_config('gene', foreground=PALETTE['accent_text'])
        tags = {'error': TEXTS['preflight_tag_error'], 'warning': TEXTS['preflight_tag_warning'],
                'info': TEXTS['preflight_tag_info']}
        for gene, issues in report.by_gene().items():
            box.insert('end', (gene or TEXTS['preflight_general']) + "\n", 'gene')
            for issue in sorted(issues, key=lambda i: {'error': 0, 'warning': 1, 'info': 2}[i.severity]):
                box.insert('end', f"  [{tags[issue.severity]}] ", issue.severity)
                box.insert('end', issue.message(lang) + "\n")
            box.insert('end', "\n")
        box.configure(state='disabled')

        row = ctk.CTkFrame(self, fg_color='transparent')
        row.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['sm'], SPACE['lg']))
        _button(row, TEXTS['preflight_btn_continue'], lambda: self._close('continue'),
                PALETTE['warning_fill']).pack(side='right')
        _button(row, TEXTS['preflight_btn_fix'], lambda: self._close('fix'),
                PALETTE['accent_fill']).pack(side='right', padx=(0, SPACE['sm']))
        ctk.CTkLabel(row, text=TEXTS['preflight_continue_hint'], font=(FONT_UI, FONT_SIZE['sm']),
                     text_color=PALETTE['text_tertiary'], wraplength=480,
                     justify='left').pack(side='left')


def show_about(parent) -> None:
    from backend.codeml_backend import codeml_version, find_codeml
    from backend.version import version_string
    path = find_codeml()
    ver = codeml_version(path) if path else None
    text = TEXTS['about_text'].format(
        version=version_string(), codeml=path or TEXTS['about_codeml_missing'],
        codeml_version=ver or '?', python=platform.python_version(),
        platform=platform.platform())
    show_message(parent, TEXTS['about_title'], text)


def _select_all_binding(entry) -> None:
    """Ctrl+A selects the whole field (Tk's default moves to the line start)."""
    def _sel(event):
        event.widget.select_range(0, 'end')
        event.widget.icursor('end')
        return 'break'
    inner = getattr(entry, '_entry', entry)
    inner.bind('<Control-a>', _sel)
    inner.bind('<Control-A>', _sel)


def os_error_text(exc: OSError) -> str:
    """A file-system error in words, instead of '[Errno 13] Permission denied: …'."""
    path = getattr(exc, 'filename', None) or ''
    if isinstance(exc, PermissionError):
        return TEXTS['err_permission'].format(path=path)
    if isinstance(exc, FileNotFoundError):
        return TEXTS['err_not_found'].format(path=path)
    if getattr(exc, 'errno', None) == 28:
        return TEXTS['err_disk_full'].format(path=path)
    return str(exc)


class Spinner:
    """A turning arc on a tk.Canvas: shows that the program is working."""

    def __init__(self, parent, size: int = 40, bg: str = None, color: str = None, width: int = 4):
        import tkinter as tk
        self.canvas = tk.Canvas(parent, width=size, height=size, highlightthickness=0,
                                bg=bg or PALETTE['bg_window'])
        pad = width + 1
        self.canvas.create_oval(pad, pad, size - pad, size - pad, outline=mix(bg or PALETTE['bg_window'],
                                PALETTE['text_tertiary'], 0.35), width=width)
        self.arc = self.canvas.create_arc(pad, pad, size - pad, size - pad, start=90, extent=100,
                                          style='arc', outline=color or PALETTE['accent_fill'], width=width)
        self.angle = 90
        self.running = True
        self._job = None
        self._tick()

    def _tick(self):
        if not self.running:
            return
        self.angle = (self.angle - 12) % 360
        try:
            self.canvas.itemconfigure(self.arc, start=self.angle)
            self._job = self.canvas.after(30, self._tick)
        except Exception:
            self.running = False

    def stop(self):
        self.running = False
        if self._job is not None:
            try:
                self.canvas.after_cancel(self._job)
            except Exception:
                pass
            self._job = None


class LoadingOverlay:
    """Covers a window with a spinner and a line of text until close() is called."""

    def __init__(self, parent, text: str):
        self.frame = ctk.CTkFrame(parent, fg_color=PALETTE['bg_window'], corner_radius=0)
        self.frame.place(relx=0, rely=0, relwidth=1, relheight=1)
        box = ctk.CTkFrame(self.frame, fg_color='transparent')
        box.place(relx=0.5, rely=0.45, anchor='center')
        self.spinner = Spinner(box, size=44, bg=PALETTE['bg_window'])
        self.spinner.canvas.pack()
        ctk.CTkLabel(box, text=text, font=(FONT_UI, FONT_SIZE['md']),
                     text_color=PALETTE['text_secondary']).pack(pady=(SPACE['md'], 0))
        self.frame.lift()

    def close(self):
        self.spinner.stop()
        try:
            self.frame.destroy()
        except Exception:
            pass


def add_tooltip(widget, text: str, wraplength: int = 320) -> None:
    """Small box with text while the mouse is over widget."""
    tip = {}

    def _show(event):
        if tip or not text:
            return
        w = ctk.CTkToplevel(widget)
        w.overrideredirect(True)
        w.attributes('-topmost', True)
        ctk.CTkLabel(w, text=text, font=(FONT_UI, FONT_SIZE['sm']), wraplength=wraplength,
                     justify='left', fg_color=PALETTE['bg_elevated'], corner_radius=6,
                     text_color=PALETTE['text_primary'], padx=SPACE['sm'], pady=SPACE['xs']).pack()
        w.geometry(f"+{event.x_root + 12}+{event.y_root + 16}")
        tip['w'] = w

    def _hide(_event=None):
        w = tip.pop('w', None)
        if w is not None:
            w.destroy()

    widget.bind('<Enter>', _show, add='+')
    widget.bind('<Leave>', _hide, add='+')


def ask_string(parent, title: str, prompt: str, initial: str = '') -> Optional[str]:
    """One-line text input in the program's theme."""
    dlg = _Modal(parent, title, 440, 200)
    ctk.CTkLabel(dlg, text=prompt, font=(FONT_UI, FONT_SIZE['md']), wraplength=390, justify='left',
                 text_color=PALETTE['text_primary']).pack(padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['sm']),
                                                          anchor='w')
    var = ctk.StringVar(value=initial)
    entry = ctk.CTkEntry(dlg, textvariable=var, height=32, fg_color=PALETTE['bg_inset'],
                         border_color=PALETTE['control_border'], text_color=PALETTE['text_primary'],
                         corner_radius=RADIUS['field'], font=(FONT_UI, FONT_SIZE['md']))
    entry.pack(fill='x', padx=SPACE['xl'])
    _select_all_binding(entry)
    row = ctk.CTkFrame(dlg, fg_color='transparent')
    row.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['md'], SPACE['xl']), side='bottom')
    _button(row, "OK", lambda: dlg._close(var.get()), PALETTE['accent_fill'], width=90).pack(side='right')
    _button(row, TEXTS['picker_cancel'], lambda: dlg._close(None), PALETTE['neutral_fill']).pack(
        side='right', padx=(0, SPACE['sm']))
    entry.bind('<Return>', lambda e: dlg._close(var.get()))
    dlg.after(80, entry.focus_set)
    return dlg.show()


def _patterns(filetypes) -> List[str]:
    pats = []
    for _, spec in (filetypes or []):
        pats += spec.split()
    return [p for p in pats if p not in ('*', '*.*')] if pats else []


class FilePicker(_Modal):
    """Themed folder and file chooser, used on Linux instead of Tk's dialogs.
    mode 'dir' picks a folder, 'open' an existing file, 'save' a new file name."""

    def __init__(self, parent, title: str, initialdir=None, mode: str = 'dir',
                 allow_new: bool = False, must_exist: bool = True, filetypes=None,
                 initialfile: str = '', defaultextension: str = ''):
        super().__init__(parent, title, 720, 540)
        self.mode, self.allow_new, self.must_exist = mode, allow_new, must_exist
        self.patterns = _patterns(filetypes)
        self.defaultextension = defaultextension
        start = Path(initialdir).expanduser() if initialdir else Path.home()
        while not start.is_dir() and start != start.parent:
            start = start.parent
        self.cwd = start.resolve()
        self._entries = []
        pad = SPACE['lg']

        ctk.CTkLabel(self, text=title, font=(FONT_UI, FONT_SIZE['md'], 'bold'), anchor='w',
                     justify='left', wraplength=660, text_color=PALETTE['text_primary']
                     ).pack(fill='x', padx=pad, pady=(pad, SPACE['sm']))
        bar = ctk.CTkFrame(self, fg_color='transparent')
        bar.pack(fill='x', padx=pad)
        back = ctk.CTkButton(bar, text="←", width=36, height=30, command=self._up,
                             fg_color=PALETTE['neutral_fill'],
                             hover_color=mix(PALETTE['neutral_fill'], '#000000', 0.2),
                             text_color='#ffffff', corner_radius=RADIUS['field'],
                             font=(FONT_UI, FONT_SIZE['lg'], 'bold'))
        back.pack(side='left')
        add_tooltip(back, TEXTS['picker_up'])
        self.path_var = ctk.StringVar(value=str(self.cwd))
        self.path_entry = ctk.CTkEntry(bar, textvariable=self.path_var, height=30,
                                       fg_color=PALETTE['bg_inset'], border_color=PALETTE['control_border'],
                                       text_color=PALETTE['text_primary'], corner_radius=RADIUS['field'],
                                       font=(FONT_UI, FONT_SIZE['sm']))
        self.path_entry.pack(side='left', fill='x', expand=True, padx=(SPACE['sm'], 0))
        self.path_entry.bind('<Return>', lambda e: self._typed_path())
        _select_all_binding(self.path_entry)
        self.path_error = ctk.CTkLabel(self, text='', font=(FONT_UI, FONT_SIZE['xs']), anchor='w',
                                       text_color=PALETTE['danger_text'], height=0)

        import tkinter as tk
        frame = ctk.CTkFrame(self, fg_color=PALETTE['bg_inset'], corner_radius=RADIUS['card'])
        frame.pack(fill='both', expand=True, padx=pad, pady=SPACE['sm'])
        self.listbox = tk.Listbox(frame, activestyle='none', borderwidth=0, highlightthickness=0,
                                  bg=PALETTE['bg_inset'], fg=PALETTE['text_primary'],
                                  selectbackground=PALETTE['accent_fill'], selectforeground='#ffffff',
                                  font=(FONT_UI, FONT_SIZE['sm']), exportselection=False)
        sb = ctk.CTkScrollbar(frame, command=self.listbox.yview)
        self.listbox.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y', pady=SPACE['xs'])
        self.listbox.pack(side='left', fill='both', expand=True, padx=SPACE['sm'], pady=SPACE['sm'])
        self.listbox.bind('<Double-Button-1>', lambda e: self._open_selected())
        self.listbox.bind('<Return>', lambda e: self._open_selected())
        self.listbox.bind('<BackSpace>', lambda e: self._up())
        self.listbox.bind('<<ListboxSelect>>', lambda e: self._on_select())

        self.name_var = ctk.StringVar(value=initialfile)
        if mode == 'save':
            name_row = ctk.CTkFrame(self, fg_color='transparent')
            name_row.pack(fill='x', padx=pad, pady=(0, SPACE['xs']))
            ctk.CTkLabel(name_row, text=TEXTS['picker_file_name'], font=(FONT_UI, FONT_SIZE['sm'], 'bold'),
                         text_color=PALETTE['text_primary']).pack(side='left', padx=(0, SPACE['sm']))
            self.name_entry = ctk.CTkEntry(name_row, textvariable=self.name_var, height=30,
                                           fg_color=PALETTE['bg_inset'], border_color=PALETTE['control_border'],
                                           text_color=PALETTE['text_primary'], corner_radius=RADIUS['field'],
                                           font=(FONT_UI, FONT_SIZE['sm']))
            self.name_entry.pack(side='left', fill='x', expand=True)
            self.name_entry.bind('<Return>', lambda e: self._choose())
            _select_all_binding(self.name_entry)
        hint = {'dir': 'picker_hint', 'open': 'picker_hint_open', 'save': 'picker_hint_save'}[mode]
        ctk.CTkLabel(self, text=TEXTS[hint], font=(FONT_UI, FONT_SIZE['xs']), anchor='w',
                     justify='left', wraplength=660, text_color=PALETTE['text_secondary']
                     ).pack(fill='x', padx=pad)
        row = ctk.CTkFrame(self, fg_color='transparent')
        row.pack(fill='x', padx=pad, pady=(SPACE['sm'], pad))
        if allow_new or mode == 'save':
            _button(row, TEXTS['picker_new_folder'], self._new_folder, PALETTE['neutral_fill']).pack(side='left')
        self.choose_btn = _button(row, TEXTS['picker_choose'], self._choose, PALETTE['accent_fill'])
        self.choose_btn.pack(side='right')
        _button(row, TEXTS['picker_cancel'], lambda: self._close(None), PALETTE['neutral_fill']).pack(
            side='right', padx=(0, SPACE['sm']))
        self._load()

    def _matches(self, p: Path) -> bool:
        if not self.patterns:
            return True
        return any(fnmatch.fnmatch(p.name.lower(), pat.lower()) for pat in self.patterns)

    def _load(self):
        self.path_var.set(str(self.cwd))
        self.path_entry.configure(border_color=PALETTE['control_border'])
        self.path_error.pack_forget()
        self.listbox.delete(0, 'end')
        self._entries = []
        try:
            items = sorted(self.cwd.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower()))
        except OSError:
            items = []
        for p in items:
            if p.name.startswith('.'):
                continue
            is_dir = p.is_dir()
            pickable = is_dir if self.mode == 'dir' else (is_dir or self._matches(p))
            self.listbox.insert('end', (p.name + '/') if is_dir else '    ' + p.name)
            if not pickable:
                self.listbox.itemconfig('end', fg=PALETTE['text_tertiary'],
                                        selectbackground=PALETTE['bg_inset'],
                                        selectforeground=PALETTE['text_tertiary'])
            self._entries.append((p, is_dir, pickable))
            if len(self._entries) >= 2000:
                break
        self._refresh_choose()

    def _selected(self):
        sel = self.listbox.curselection()
        return self._entries[sel[0]] if sel else None

    def _selected_dir(self):
        e = self._selected()
        return e[0] if e and e[1] else None

    def _selected_file(self):
        e = self._selected()
        return e[0] if e and not e[1] and e[2] else None

    def _on_select(self):
        f = self._selected_file()
        if f is not None and self.mode == 'save':
            self.name_var.set(f.name)
        self._refresh_choose()

    def _refresh_choose(self):
        if self.mode == 'dir':
            target = self._selected_dir() or self.cwd
            name = target.name or str(target)
            self.choose_btn.configure(text=TEXTS['picker_choose_named'].format(name=name))
        elif self.mode == 'open':
            f = self._selected_file()
            self.choose_btn.configure(text=TEXTS['picker_choose_named'].format(name=f.name) if f
                                      else TEXTS['picker_open'], state='normal' if f else 'disabled')
        else:
            self.choose_btn.configure(text=TEXTS['picker_save'])

    def _open_selected(self):
        e = self._selected()
        if e is None:
            return
        if e[1]:
            self.cwd = e[0].resolve()
            self._load()
        elif self.mode == 'open' and e[2]:
            self._close(str(e[0].resolve()))
        elif self.mode == 'save':
            self._choose()

    def _up(self):
        if self.cwd.parent != self.cwd:
            self.cwd = self.cwd.parent
            self._load()

    def _show_path_error(self, text: str):
        self.path_entry.configure(border_color=PALETTE['danger_fg'])
        self.path_error.configure(text=text)
        self.path_error.pack(fill='x', padx=SPACE['lg'], after=self.path_entry.master)

    def _typed_path(self):
        p = Path(self.path_var.get().strip()).expanduser()
        if p.is_dir():
            self.cwd = p.resolve()
            self._load()
        elif self.mode == 'open' and p.is_file():
            self._close(str(p.resolve()))
        elif self.mode == 'dir' and not self.must_exist and p.is_absolute() and p.parent.is_dir():
            self._close(str(p))
        else:
            self._show_path_error(TEXTS['picker_path_missing'].format(path=p))

    def _new_folder(self):
        name = (ask_string(self, TEXTS['picker_new_folder'], TEXTS['picker_new_folder_prompt']) or '').strip()
        if not name or '/' in name or name in ('.', '..'):
            return
        new = self.cwd / name
        try:
            new.mkdir(exist_ok=True)
        except OSError as exc:
            show_message(self, TEXTS['picker_new_folder'], os_error_text(exc), 'error')
            return
        self.cwd = new.resolve()
        self._load()

    def _choose(self):
        if self.mode == 'open':
            f = self._selected_file()
            if f is not None:
                self._close(str(f.resolve()))
            return
        if self.mode == 'save':
            name = self.name_var.get().strip()
            if not name or '/' in name:
                return
            if self.defaultextension and not Path(name).suffix:
                name += self.defaultextension
            target = self.cwd / name
            if target.exists() and not ask_yes_no(self, TEXTS['picker_save'],
                                                  TEXTS['picker_overwrite'].format(name=name)):
                return
            self._close(str(target))
            return
        typed = Path(self.path_var.get().strip()).expanduser()
        if typed != self.cwd and typed.is_absolute() and not self._selected_dir():
            if typed.is_dir() or (not self.must_exist and typed.parent.is_dir()):
                self._close(str(typed))
                return
        self._close(str((self._selected_dir() or self.cwd).resolve()))


def ask_directory(parent, title: str, initialdir=None, allow_new: bool = False,
                  must_exist: bool = True):
    """Absolute path of the chosen folder, or None. FilePicker on Linux, the
    native dialog on Windows and macOS."""
    if platform.system() == 'Linux':
        return FilePicker(parent, title, initialdir, 'dir', allow_new, must_exist).show()
    from tkinter import filedialog
    path = filedialog.askdirectory(parent=parent, initialdir=str(initialdir or Path.home()),
                                   title=title, mustexist=must_exist)
    return str(Path(path).resolve()) if path else None


def ask_open_file(parent, title: str, initialdir=None, filetypes=None):
    if platform.system() == 'Linux':
        return FilePicker(parent, title, initialdir, 'open', filetypes=filetypes).show()
    from tkinter import filedialog
    path = filedialog.askopenfilename(parent=parent, initialdir=str(initialdir or Path.home()),
                                      title=title, filetypes=filetypes or [('*', '*.*')])
    return path or None


def ask_save_file(parent, title: str, initialdir=None, initialfile: str = '',
                  defaultextension: str = '', filetypes=None):
    if platform.system() == 'Linux':
        return FilePicker(parent, title, initialdir, 'save', filetypes=filetypes,
                          initialfile=initialfile, defaultextension=defaultextension).show()
    from tkinter import filedialog
    path = filedialog.asksaveasfilename(parent=parent, initialdir=str(initialdir or Path.home()),
                                        title=title, initialfile=initialfile,
                                        defaultextension=defaultextension,
                                        filetypes=filetypes or [('*', '*.*')])
    return path or None
